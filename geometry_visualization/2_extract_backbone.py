"""
Script: Extract PEG backbone atoms from PRg chain data
------------------------------------------------------
This script extracts backbone atoms (C and O types) from
the PEG portion of a PRg model.

Input:
    - PEG-only LAMMPS data file (e.g., 1_peg_xxx.data)

Output:
    - A simplified file containing backbone atom coordinates
      grouped by molecule ID.

Parameters:
    - target_atom_types = {5, 73}   # atom types corresponding to C and O atoms

"""

import os
import sys
from collections import defaultdict

# ====================== User Parameters ======================
target_atom_types = {5, 73}   # atom types to retain (C and O)
# =============================================================


def parse_atoms_section(lines):
    """
    Extract atom data from the 'Atoms' section of a LAMMPS data file.
    Returns a list of tuples: (atom_id, mol_id, atom_type, x, y, z)
    """
    atoms_data = []
    reading_atoms = False
    section_headers = {"Bonds", "Angles", "Dihedrals", "Impropers", "Velocities"}

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Atoms"):
            reading_atoms = True
            continue
        if reading_atoms:
            if stripped in section_headers or stripped.startswith(tuple(section_headers)):
                break
            parts = stripped.split()
            if len(parts) < 7:
                continue
            try:
                atom_id = int(parts[0])
                mol_id = int(parts[1])
                atom_type = int(parts[2])
                x, y, z = map(float, parts[4:7])
                if atom_type in target_atom_types:
                    atoms_data.append((atom_id, mol_id, atom_type, x, y, z))
            except ValueError:
                continue

    print(f"✅ Extracted {len(atoms_data)} backbone atoms matching {target_atom_types}.")
    return atoms_data


def group_atoms_by_molecule(atoms_data):
    """Group atom records by molecule ID and sort each group by atom ID."""
    grouped = defaultdict(list)
    for entry in atoms_data:
        atom_id, mol_id, atom_type, x, y, z = entry
        grouped[mol_id].append(entry)
    for mol_id in grouped:
        grouped[mol_id].sort(key=lambda x: x[0])
    return grouped


def write_backbone_file(grouped_atoms, output_path):
    """Write extracted backbone coordinates grouped by molecule."""
    new_atom_id = 1
    with open(output_path, "w") as out:
        for mol_id in sorted(grouped_atoms.keys()):
            out.write(f"molecule_{mol_id}\n")
            for atom in grouped_atoms[mol_id]:
                _, _, atom_type, x, y, z = atom
                out.write(f"{mol_id} {new_atom_id} {atom_type} {x:.6f} {y:.6f} {z:.6f}\n")
                new_atom_id += 1
            out.write("\n")
    print(f"✅ Backbone coordinates written to: {output_path}")


def main():
    if len(sys.argv) != 3:
        print("\nUsage: python3 2_extract_backbone.py <input_PEG_file> <output_file>\n")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(current_dir, input_file)
    output_path = os.path.join(current_dir, output_file)

    with open(input_path, "r") as f:
        lines = f.readlines()

    atoms_data = parse_atoms_section(lines)
    grouped_atoms = group_atoms_by_molecule(atoms_data)
    write_backbone_file(grouped_atoms, output_path)


if __name__ == "__main__":
    main()
