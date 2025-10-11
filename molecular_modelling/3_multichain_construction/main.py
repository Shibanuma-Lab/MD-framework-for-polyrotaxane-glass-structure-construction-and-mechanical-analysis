"""
Script: Construction of a periodic multi-chain simulation cell
--------------------------------------------------------------
This script constructs a periodic multi-chain PRg (polyrotaxane glass) system
by replicating and packing multiple single-chain configurations into a defined
simulation box.
"""

from __future__ import annotations
import os
from read_lammps_data import read_lammps_data
from generate_molecules import generate_molecules
from write_lammps_data import write_lammps_data

# ====================== User Parameters ======================
num_chains = 10                       # number of PRg chains to generate
box_half = 196.228301691482           # half box length (Å)
min_distance = 0.1                    # minimum allowed intermolecular distance
max_attempts = 100                    # maximum trial attempts per molecule
input_filename = "prg_singlechain.data"     # input single-chain data file
output_filename = "prg_systems.data" # output multi-chain data file
# =============================================================


def main() -> None:
    """Generate a multi-chain PRg system from a single-chain LAMMPS data file."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(current_dir, input_filename)
    output_path = os.path.join(current_dir, output_filename)
    box_size = [-box_half, box_half, -box_half, box_half, -box_half, box_half]

    template = read_lammps_data(input_path)
    generated = generate_molecules(
        template,
        num_chains,
        box_size,
        min_distance=min_distance,
        max_attempts=max_attempts,
    )
    write_lammps_data(output_path, template, box_size, generated)
    print(f"✅ Done. Wrote multi-chain data file: {output_path}")


if __name__ == "__main__":
    main()
