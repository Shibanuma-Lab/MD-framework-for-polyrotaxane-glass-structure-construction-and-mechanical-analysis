"""
Script: Automated multi-stage annealing and equilibration of polymer systems
----------------------------------------------------------------------------
This script performs a multi-stage annealing and equilibration procedure
for polymeric systems (e.g., PRg). The simulation proceeds through a series
of NPT stages, each automatically terminated based on a density-based
convergence criterion.

Convergence criterion:
    1. Collect instantaneous density over a sliding window of length `history_length`
    2. Fit a linear regression to (time, density) data within the window
    3. Compute the slope per picosecond (slope_per_ps)
    4. If |slope_per_ps| < slope_tolerance holds for `consecutive_needed` checks,
       the stage is considered converged and the simulation proceeds to the next stage.
"""

from mpi4py import MPI
from lammps import lammps
import numpy as np
from collections import deque

# ========================== User Parameters ==========================
input_data = "prg_system.data"          # input LAMMPS data file
log_file = "anneal_prg_system.lammps"
max_total_steps = 1000000                   # global maximum steps
check_interval = 50                        # steps per run chunk
consecutive_needed = 200                   # consecutive satisfied checks required

DEFAULT_STAGE_SETTINGS = {
    "timestep": 0.5,                       # fs
    "slope_tolerance": 1.0e-3,             # density/ps
    "history_length": 401,                 # sliding window length
}
# =====================================================================


# ---------------------- MPI and Initialization -----------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
lmp = lammps(comm=comm)

if rank == 0:
    print("Initializing LAMMPS with GPU support...")

lmp.command("package gpu 1 neigh no")
lmp.command("suffix gpu")
lmp.command("units real")
lmp.command("atom_style full")
lmp.command("bond_style hybrid harmonic")
lmp.command("angle_style hybrid harmonic")
lmp.command("dihedral_style hybrid fourier")
lmp.command("improper_style hybrid cvff")
lmp.command("pair_style hybrid lj/charmm/coul/long 9.0 10.0 10.0")
lmp.command("kspace_style pppm 0.0001")
lmp.command("pair_modify mix arithmetic")
lmp.command("special_bonds amber")

lmp.command(f"read_data {input_data}")
lmp.command("include system.in.settings")
lmp.command(f"log {log_file}")
lmp.command("thermo 50")
lmp.command("thermo_style custom step time temp press pe ke etotal ebond vol density")

# Pressure variables
lmp.command("variable Pmax equal 28000")
lmp.command("variable 001P equal 0.01*${Pmax}")
lmp.command("variable 01P  equal 0.1*${Pmax}")
lmp.command("variable 05P  equal 0.5*${Pmax}")

lmp.command("minimize 0.0 0.0 1000 100000")
lmp.command("reset_timestep 0")
lmp.command("velocity all create 300 32479 dist gaussian loop geom")
# ---------------------------------------------------------------------


# -------------------------- Stage Settings ---------------------------
stages = [
    {"name": "Step1",  "fix": "fix 1 all npt temp 300.0 300.0 $(100*dt) iso 1.0 1.0 $(1000*dt)"},
    {"name": "Step2",  "fix": "fix 1 all npt temp 1000.0 1000.0 $(100*dt) iso ${001P} ${001P} $(1000*dt)"},
    {"name": "Step3",  "fix": "fix 1 all npt temp 1000.0 1000.0 $(100*dt) iso ${01P} ${01P} $(1000*dt)"},
    {"name": "Step4",  "fix": "fix 1 all npt temp 1000.0 1000.0 $(100*dt) iso ${05P} ${05P} $(1000*dt)"},
    {"name": "Step5",  "fix": "fix 1 all npt temp 1000.0 1000.0 $(100*dt) iso ${Pmax} ${Pmax} $(1000*dt)"},
    {"name": "Step6",  "fix": "fix 1 all npt temp 300.0 300.0 $(100*dt) iso ${Pmax} ${Pmax} $(1000*dt)"},
    {"name": "Step7",  "fix": "fix 1 all npt temp 300.0 300.0 $(100*dt) iso ${05P} ${05P} $(1000*dt)"},
    {"name": "Step8",  "fix": "fix 1 all npt temp 300.0 300.0 $(100*dt) iso ${01P} ${01P} $(1000*dt)"},
    {"name": "Step9",  "fix": "fix 1 all npt temp 300.0 300.0 $(100*dt) iso ${001P} ${001P} $(1000*dt)"},
    {"name": "Step10", "fix": "fix 1 all npt temp 300.0 300.0 $(100*dt) iso 1.0 1.0 $(1000*dt)"},
]
# ---------------------------------------------------------------------


# ------------------------ Utility Functions --------------------------
def compute_slope(data_pairs):
    """Return slope (density/ps) from linear regression of (time, density)."""
    if len(data_pairs) < 2:
        return 0.0
    x, y = zip(*data_pairs)
    slope, _ = np.polyfit(x, y, 1)
    return slope


def monitor_and_run(lmp, stage, check_interval, total_steps_run, time_ps):
    """Run simulation for one stage with density-based convergence check."""
    settings = DEFAULT_STAGE_SETTINGS
    name = stage["name"]
    fix_cmd = stage["fix"]
    dt_fs = settings["timestep"]
    slope_tol = settings["slope_tolerance"]
    history_len = settings["history_length"]

    lmp.command(f"timestep {dt_fs}")
    lmp.command(fix_cmd)
    if rank == 0:
        print(f"[{name}] start | fix={fix_cmd} | dt={dt_fs} fs | tol={slope_tol}")

    data_window = deque(maxlen=history_len)
    satisfied = 0
    stabilized = False

    while total_steps_run < max_total_steps:
        chunk = min(check_interval, max_total_steps - total_steps_run)
        lmp.command(f"run {chunk} post no")
        total_steps_run += chunk
        time_ps += (chunk * dt_fs) / 1000.0

        density = lmp.get_thermo("density")
        data_window.append((time_ps, density))

        if rank == 0:
            print(f"[{name}] step={total_steps_run}, time={time_ps:.3f} ps, density={density:.5f}")

        if len(data_window) == history_len:
            slope = compute_slope(list(data_window))
            if abs(slope) < slope_tol:
                satisfied += 1
                if rank == 0:
                    print(f"[{name}] slope OK {satisfied}/{consecutive_needed}")
            else:
                satisfied = 0
            if satisfied >= consecutive_needed:
                stabilized = True
                break

    lmp.command("unfix 1")
    lmp.command(f"write_data {name}-final.data")

    if rank == 0:
        msg = "stabilized" if stabilized else "reached step limit"
        print(f"[{name}] finished ({msg}).")

    return total_steps_run, time_ps
# ---------------------------------------------------------------------


# ----------------------------- Main Loop -----------------------------
def main():
    total_steps_run = 0
    time_ps = 0.0

    for stage in stages:
        if rank == 0:
            print(f"\n===== Starting {stage['name']} =====\n")
        if total_steps_run >= max_total_steps:
            if rank == 0:
                print("Global step limit reached. Stopping.")
            break
        total_steps_run, time_ps = monitor_and_run(lmp, stage, check_interval, total_steps_run, time_ps)

    if rank == 0:
        print("\nAll stages completed.\n")


if __name__ == "__main__":
    main()
