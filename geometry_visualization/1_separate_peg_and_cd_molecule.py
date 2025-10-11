"""
Script: Separate chain and ring molecules from PRg data file
------------------------------------------------------------
This script separates the polymeric (PEG) chains and cyclic (CD) rings
from a LAMMPS data file of polyrotaxane (PRg).

Output:
    1. PEG atoms file  (end groups removed)
    2. CD atoms file   (cyclic molecules)

Input:
    - {specific_model}.data  : original LAMMPS data file

Parameters:
    - molecules_per_chain = 47   # each PRg chain: 1 PEG + multiple CDs
    - remove 35 atoms from head and 29 from tail for each PEG chain

"""

import os
import sys
from collections import defaultdict

# ====================== User Parameters ======================
molecules_per_chain = 47   # 1 PEG + multiple CD per PRg chain
remove_front = 35          # atoms to remove from PEG head
remove_tail = 29           # atoms to remove from PEG tail
# =============================================================


def split_data_file(input_file, output_peg, output_cd):
    """Separate PEG and CD molecules from a LAMMPS data file."""
    with open(input_file, "r") as f:
        lines = f.readlines()

    header, masses_section, atoms_section = [], [], []
    in_atoms = in_masses = False
    atoms_started = False

    for line in lines:
        if line.startswith("Atoms"):
            in_atoms, in_masses = True, False
            atoms_section.append(line)
            atoms_started = False
            continue
        elif line.startswith("Masses"):
            in_masses, in_atoms = True, False
            masses_section.append(line)
            continue
        elif line.strip() == "":
            if in_atoms:
                atoms_section.append(line)
            elif in_masses:
                masses_section.append(line)
            else:
                header.append(line)
            continue
        elif in_atoms:
            if not atoms_started:
                atoms_section.append(line)
                atoms_started = True
            elif line.startswith(("Bonds", "Angles", "Velocities", "Impropers", "Dihedrals")):
                break
            else:
                atoms_section.append(line)
        elif in_masses:
            masses_section.append(line)
        else:
            header.append(line)

    # Extract atom data
    atom_data = [line for line in atoms_section if line.strip() and line.strip().split()[0].isdigit()]
    mol_to_atoms = defaultdict(list)
    for line in atom_data:
        parts = line.strip().split()
        mol_id = int(parts[1])
        mol_to_atoms[mol_id].append(parts)

    # Classify PEG and CD molecules
    peg_mols, cd_mols = [], []
    for mol_id in sorted(mol_to_atoms.keys()):
        if (mol_id - 1) % molecules_per_chain == 0:
            peg_mols.append(mol_id)
        else:
            cd_mols.append(mol_id)

    # --- Write PEG data ---
    mol_remap_peg = {old: new for new, old in enumerate(sorted(peg_mols), start=1)}
    with open(output_peg, "w") as f:
        f.writelines(header + masses_section)
        f.write("Atoms\n\n")
        for old_mol in peg_mols:
            atoms = sorted(mol_to_atoms[old_mol], key=lambda x: int(x[0]))
            trimmed = atoms[remove_front:-remove_tail] if len(atoms) > (remove_front + remove_tail) else []
            for atom in trimmed:
                atom[1] = str(mol_remap_peg[old_mol])
                f.write(" ".join(atom) + "\n")

    # --- Write CD data ---
    mol_remap_cd = {old: new for new, old in enumerate(sorted(cd_mols), start=1)}
    with open(output_cd, "w") as f:
        f.writelines(header + masses_section)
        f.write("Atoms\n\n")
        for old_mol in cd_mols:
            for atom in mol_to_atoms[old_mol]:
                atom[1] = str(mol_remap_cd[old_mol])
                f.write(" ".join(atom) + "\n")

    print(f"✅ PEG atoms written to: {output_peg}")
    print(f"✅ CD atoms written to:  {output_cd}")


def main():
    if len(sys.argv) != 4:
        print("\nUsage: python3 1_separate_chain_and_ring.py <input> <peg_output> <cd_output>\n")
        sys.exit(1)

    input_name, peg_name, cd_name = sys.argv[1:4]
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(current_dir, f"{input_name}.data")
    output_peg = os.path.join(current_dir, peg_name)
    output_cd = os.path.join(current_dir, cd_name)

    split_data_file(input_path, output_peg, output_cd)


if __name__ == "__main__":
    main()
