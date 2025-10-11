"""
Build a LAMMPS input file (combine_prg.in) to assemble a PEG chain with CD molecules.
- Extract two anchor x-coordinates from a PEG data file.
- Place `num_cd` CDs evenly between the two anchors along +x (y=z=0).
"""

from __future__ import annotations
import os
import random
from typing import List, Tuple, Optional

# -------------------- User parameters --------------------
front_monomer = 5                 # PEG monomer index from the front (1-based)
back_monomer = 6                  # PEG monomer index from the back  (1-based)
num_cd = 91                       # number of CD molecules to insert
input_filename = "peg_n728.data"  # LAMMPS data file containing the PEG
# ---------------------------------------------------------


def extract_x_coordinates(
    filename: str,
    front_mono: int,
    back_mono: int
) -> Tuple[float, float]:
    """
    Parse 'Atoms' section of a LAMMPS data file and return two x-coordinates:
      - x of the O atom in monomer #front_mono from the front
      - x of the O atom in monomer #back_mono from the back
    Assumptions:
      * Atom indexing matches the PEG layout used by upstream scripts:
        - Left end occupies atom IDs 1..35
        - PEG chain starts at atom ID 36
        - Each monomer has 7 atoms; the O atom has local ID = 1
        - Right end occupies the last 29 atoms
      * 'Atoms' lines use: id mol type charge x y z ...
    """
    with open(filename, "r") as fh:
        lines = fh.readlines()

    # Collect raw atom rows inside "Atoms" block.
    atoms_block: List[List[str]] = []
    in_atoms = False
    for line in lines:
        if line.strip().startswith("Atoms"):
            in_atoms = True
            continue
        if not in_atoms:
            continue
        if not line.strip():
            # skip blank lines
            continue
        # stop when hitting a non-numeric row (e.g., "Velocities", "Bonds", etc.)
        if any(c.isalpha() for c in line.split()[0]):
            break
        atoms_block.append(line.split())

    if not atoms_block:
        raise ValueError("No 'Atoms' section found or it is empty.")

    total_atoms = len(atoms_block)

    # Compute global atom IDs for the two O atoms (local O id = 1 in each monomer)
    left_end_atoms = 35
    monomer_atoms = 7
    right_end_atoms = 29

    first_atom_id = left_end_atoms + (front_mono - 1) * monomer_atoms + 1
    second_atom_id = total_atoms - right_end_atoms - (back_mono - 1) * monomer_atoms + 1  # local O=1

    first_x: Optional[float] = None
    second_x: Optional[float] = None

    for row in atoms_block:
        # columns: id mol type charge x y z ...
        try:
            atom_id = int(row[0])
            x_coord = float(row[4])
        except (ValueError, IndexError):
            continue
        if atom_id == first_atom_id:
            first_x = x_coord
        if atom_id == second_atom_id:
            second_x = x_coord
        if first_x is not None and second_x is not None:
            break

    if first_x is None or second_x is None:
        raise ValueError("Failed to locate one or both anchor atoms in 'Atoms' section.")
    return first_x, second_x


def linspace(a: float, b: float, n: int) -> List[float]:
    """Minimal replacement for numpy.linspace(a, b, n)."""
    if n <= 1:
        return [a]
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def generate_cd_anchor_points(x1: float, x2: float, n_cd: int) -> List[Tuple[float, float, float]]:
    """Evenly spaced anchors on x in [x1, x2], y=z=0."""
    return [(x, 0.0, 0.0) for x in linspace(x1, x2, n_cd)]


def generate_lammps_in_file(
    anchors: List[Tuple[float, float, float]],
    output_filename: str,
    data_filename: str
) -> None:
    """Write a LAMMPS input that reads PEG data and appends CDs at anchor positions."""
    header = f"""units           real
atom_style      full
bond_style      hybrid harmonic
angle_style     hybrid harmonic
dihedral_style  hybrid fourier
improper_style  hybrid cvff
pair_style      hybrid lj/charmm/coul/long 9.0 10.0 10.0
kspace_style    pppm 0.0001
pair_modify     mix arithmetic
special_bonds   amber

# ----------------- Atom Definition Section -----------------
read_data "{data_filename}"
include "system.in.settings"
"""

    with open(output_filename, "w") as fh:
        fh.write(header)
        for x, y, z in anchors:
            # Randomly pick a CD variant (1 or 2)
            data_file = f"cd_modified_groups_12_{random.choice([1, 2])}.data"
            fh.write(f'read_data "{data_file}" add append shift {x:.6f} {y:.6f} {z:.6f}\n')

        fh.write("\n# ----------------- Output -----------------\n")
        fh.write("write_data prg.data\n")


if __name__ == "__main__":
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(cur_dir, input_filename)
    out_in_path = os.path.join(cur_dir, "combine_prg.in")

    x_front, x_back = extract_x_coordinates(data_path, front_monomer, back_monomer)
    print(f"Front monomer #{front_monomer} O-atom x: {x_front}")
    print(f"Back  monomer #{back_monomer} O-atom x: {x_back}")

    cd_anchors = generate_cd_anchor_points(x_front, x_back, num_cd)
    print(f"Generated {len(cd_anchors)} CD anchors along x in [{x_front}, {x_back}]")

    generate_lammps_in_file(cd_anchors, out_in_path, input_filename)
    print(f"Wrote LAMMPS input: {out_in_path}")
