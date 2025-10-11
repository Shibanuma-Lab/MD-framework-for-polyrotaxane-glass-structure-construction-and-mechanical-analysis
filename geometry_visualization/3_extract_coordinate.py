"""
Script: Extract geometric centers of cyclic molecules (CDs) from LAMMPS data file
--------------------------------------------------------------------------------
This script computes the geometric centers of all cyclodextrin (CD) molecules
from a LAMMPS-style data file using the "full" atom style.

Features:
    - Automatically extracts box size and periodic dimensions.
    - Corrects coordinates for molecules crossing periodic boundaries.
    - Outputs one line per molecule: ID x_center y_center z_center.
    - Appends box dimension information at the top of the output file.

Input:
    - CD-only LAMMPS data file (e.g., 1_cd_xxx.data)

Output:
    - CD centers coordinate file (e.g., 3_cd_centers_xxx.data)

"""

import os
import sys
import numpy as np


# ====================== Helper Functions ======================
def extract_box_size(data_file):
    """Extract simulation box dimensions [Lx, Ly, Lz] from a LAMMPS data file."""
    with open(data_file, "r") as f:
        lines = f.readlines()

    xlo = xhi = ylo = yhi = zlo = zhi = None
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 4:
            if "xlo" in parts and "xhi" in parts:
                xlo, xhi = float(parts[0]), float(parts[1])
            elif "ylo" in parts and "yhi" in parts:
                ylo, yhi = float(parts[0]), float(parts[1])
            elif "zlo" in parts and "zhi" in parts:
                zlo, zhi = float(parts[0]), float(parts[1])
        if all(v is not None for v in (xlo, ylo, zlo)):
            break

    if None in (xlo, ylo, zlo):
        raise ValueError("Box boundaries not found (xlo xhi, ylo yhi, zlo zhi).")

    return [xhi - xlo, yhi - ylo, zhi - zlo]


def apply_periodic_boundary_condition(atom_positions, box_size):
    """
    Correct atoms that cross periodic boundaries.
    Keeps all atoms of a molecule in the same reference frame.
    """
    ref = atom_positions[0]
    groups = {k: [] for k in [
        "no", "x", "y", "z", "xy", "xz", "yz", "xyz"
    ]}

    for atom in atom_positions:
        cross = [abs(atom[i] - ref[i]) >= box_size[i] / 2 for i in range(3)]
        key = "".join(ax for ax, c in zip("xyz", cross) if c) or "no"
        groups[key].append(atom)

    # Determine main side (largest group)
    main_group = max(groups, key=lambda k: len(groups[k]))
    main_atoms = groups[main_group]

    # Shift minority atoms back to main side
    main_ref = main_atoms[0]
    corrected = []
    for group, atoms in groups.items():
        if group == main_group:
            corrected.extend(atoms)
            continue
        for atom in atoms:
            shifted = []
            for i in range(3):
                delta = atom[i] - main_ref[i]
                if delta >= box_size[i] / 2:
                    shifted.append(atom[i] - box_size[i])
                elif delta <= -box_size[i] / 2:
                    shifted.append(atom[i] + box_size[i])
                else:
                    shifted.append(atom[i])
            corrected.append(shifted)

    return np.array(corrected)


def extract_molecule_centers(data_file, box_size):
    """Compute geometric centers of all molecules in a data file."""
    with open(data_file, "r") as f:
        lines = f.readlines()

    atoms_section = False
    mol_atoms = {}

    for line in lines:
        if "Atoms" in line:
            atoms_section = True
            continue
        if "Velocities" in line:
            break
        if atoms_section:
            parts = line.strip().split()
            if len(parts) < 7:
                continue
            try:
                _, mol_id, *_ , x, y, z = map(float, parts[:7])
            except ValueError:
                continue
            mol_id = int(mol_id)
            mol_atoms.setdefault(mol_id, []).append([x, y, z])

    if not mol_atoms:
        raise ValueError("No molecular data found in Atoms section.")

    centers = {}
    for mol_id, atoms in mol_atoms.items():
        arr = np.array(atoms)
        arr = apply_periodic_boundary_condition(arr, box_size)
        centers[mol_id] = tuple(np.mean(arr, axis=0))

    return centers


def save_centers(centers, output_file):
    """Save molecule centers to file."""
    with open(output_file, "w") as f:
        for mol_id, (x, y, z) in sorted(centers.items()):
            f.write(f"{mol_id} {x:.6f} {y:.6f} {z:.6f}\n")
    print(f"✅ Molecule centers written to: {output_file}")


def extract_box_dimensions(data_file):
    """Return formatted box boundary info string."""
    with open(data_file, "r") as f:
        lines = f.readlines()

    xlo, xhi = lines[13].split()[:2]
    ylo, yhi = lines[14].split()[:2]
    zlo, zhi = lines[15].split()[:2]
    return f"# Box dimensions: xlo={xlo} xhi={xhi} ylo={ylo} yhi={yhi} zlo={zlo} zhi={zhi}"


def prepend_box_info(box_info, output_file):
    """Prepend box dimension info to the top of the output file."""
    if not os.path.exists(output_file):
        raise FileNotFoundError(f"Output file not found: {output_file}")

    with open(output_file, "r") as f:
        content = f.readlines()

    with open(output_file, "w") as f:
        f.write(box_info + "\n")
        f.writelines(content)

    print("✅ Box dimensions added to output file.")


# ====================== Main ======================
def main():
    if len(sys.argv) != 3:
        print("\nUsage: python3 3_extract_coordinate.py <input_CD_file> <output_file>\n")
        sys.exit(1)

    input_file, output_file = sys.argv[1:3]
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(current_dir, input_file)
    output_path = os.path.join(current_dir, output_file)

    try:
        box_size = extract_box_size(input_path)
        print(f"Box size: {box_size}")
        centers = extract_molecule_centers(input_path, box_size)
        save_centers(centers, output_path)
        box_info = extract_box_dimensions(input_path)
        prepend_box_info(box_info, output_path)
    except Exception as e:
        print(f"❌ Error: {e}")


if __name__ == "__main__":
    main()
