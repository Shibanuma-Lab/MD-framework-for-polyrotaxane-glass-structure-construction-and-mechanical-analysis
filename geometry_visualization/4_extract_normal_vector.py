"""
Script: Extract orientation (normal vectors) of cyclic molecules (CDs)
----------------------------------------------------------------------
This script computes the orientation (normal vectors) of each cyclodextrin (CD)
molecule from a LAMMPS-style data file ("full" atom style).

Algorithm:
    1. Extract atoms grouped by molecule ID.
    2. Select three representative atoms (indices 9, 45, 81) per molecule.
    3. Apply periodic boundary correction to align all points within the same box.
    4. Compute the normal vector using the cross product of two spanning vectors.

Input:
    - CD-only LAMMPS data file (e.g., 1_cd_xxx.data)

Output:
    - File containing normal vectors of all CDs (e.g., 4_cd_normals_xxx.data)
      Format: mol_id nx ny nz

"""

import os
import sys
import numpy as np


# ====================== User Parameters ======================
# Indices of three atoms used to define the ring plane (1-based index)
target_indices = [9, 45, 81]
# =============================================================


def extract_box_size(data_file):
    """Extract simulation box dimensions [Lx, Ly, Lz] from the data file."""
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

    if None in (xlo, xhi, ylo, yhi, zlo, zhi):
        raise ValueError("Incomplete box boundary information found.")

    return [xhi - xlo, yhi - ylo, zhi - zlo]


def extract_atoms_by_molecule(data_file):
    """
    Read atoms from 'Atoms' section and group by molecule ID.
    Automatically handles unsorted atom IDs.
    """
    molecule_data = {}
    in_atoms = False

    with open(data_file, "r") as f:
        for line in f:
            if "Atoms" in line:
                in_atoms = True
                continue
            if "Velocities" in line:
                break
            if in_atoms:
                parts = line.strip().split()
                if len(parts) >= 7:
                    atom_id = int(parts[0])
                    mol_id = int(parts[1])
                    x, y, z = map(float, parts[4:7])
                    molecule_data.setdefault(mol_id, []).append((atom_id, x, y, z))

    # Sort atoms by atom_id within each molecule
    for mol_id in molecule_data:
        molecule_data[mol_id].sort(key=lambda a: a[0])

    return molecule_data


def select_reference_atoms(molecule_data, indices):
    """
    Select three representative atoms per molecule.
    indices: list of target atom indices (1-based)
    """
    selected = {}
    for mol_id, atoms in molecule_data.items():
        if len(atoms) >= max(indices):
            chosen = [atoms[i - 1] for i in indices]  # convert to 0-based index
            selected[mol_id] = chosen
    return selected


def apply_pbc(coords, box_size):
    """Apply periodic boundary correction to align points within one box."""
    ref = np.array(coords[0])
    corrected = [ref]
    for c in coords[1:]:
        atom = np.array(c)
        delta = atom - ref
        for i in range(3):
            if delta[i] > box_size[i] / 2:
                atom[i] -= box_size[i]
            elif delta[i] < -box_size[i] / 2:
                atom[i] += box_size[i]
        corrected.append(atom)
    return corrected


def compute_normal(atom1, atom2, atom3):
    """Compute normalized vector perpendicular to the plane defined by three atoms."""
    AB = atom2 - atom1
    AC = atom3 - atom1
    n = np.cross(AB, AC)
    norm = np.linalg.norm(n)
    return n / norm if norm > 0 else np.array([0.0, 0.0, 0.0])


def save_normals(normals, output_file):
    """Write normal vectors to output file."""
    with open(output_file, "w") as f:
        for mol_id, vec in sorted(normals.items()):
            f.write(f"{mol_id} {vec[0]:.6f} {vec[1]:.6f} {vec[2]:.6f}\n")
    print(f"✅ Normal vectors written to: {output_file}")


def main():
    if len(sys.argv) != 3:
        print("\nUsage: python3 4_extract_normal_vector.py <input_CD_file> <output_file>\n")
        sys.exit(1)

    input_file, output_file = sys.argv[1:3]
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(current_dir, input_file)
    output_path = os.path.join(current_dir, output_file)

    # 1. Extract box size
    box_size = extract_box_size(input_path)

    # 2. Group atoms by molecule ID
    molecules = extract_atoms_by_molecule(input_path)

    # 3. Select representative atoms
    selected_atoms = select_reference_atoms(molecules, target_indices)

    # 4. Compute normal vectors
    normals = {}
    for mol_id, atoms in selected_atoms.items():
        coords = [np.array(a[1:]) for a in atoms]
        corrected = apply_pbc(coords, box_size)
        normals[mol_id] = compute_normal(*corrected)

    # 5. Save results
    save_normals(normals, output_path)


if __name__ == "__main__":
    main()
