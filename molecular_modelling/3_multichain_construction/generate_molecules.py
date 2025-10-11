import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.spatial import cKDTree

def generate_molecules(template_data, num_molecules, box_size, min_distance=1.5, max_attempts=100):
    xlo, xhi, ylo, yhi, zlo, zhi = box_size
    box_lengths = [xhi - xlo, yhi - ylo, zhi - zlo]

    all_atoms = []
    all_bonds = []
    all_angles = []
    all_dihedrals = []
    all_impropers = []

    all_positions = []

    atom_offset = 0
    bond_offset = 0
    angle_offset = 0
    dihedral_offset = 0
    improper_offset = 0

    coords_template = np.array([[atom['x'], atom['y'], atom['z']] for atom in template_data['atoms']])
    center = coords_template.mean(axis=0)
    local_coords = coords_template - center

    for mol_id in range(1, num_molecules + 1):
        success = False
        for attempt in range(max_attempts):
            if attempt == 0:
                print(f"[{mol_id:3d}/{num_molecules}] Trying to insert molecule...", flush=True)

            rotation = R.random().as_matrix()
            dx = np.random.uniform(xlo, xhi)
            dy = np.random.uniform(ylo, yhi)
            dz = np.random.uniform(zlo, zhi)

            rotated = (rotation @ local_coords.T).T + np.array([dx, dy, dz])

            image_flags = np.floor((rotated - np.array([xlo, ylo, zlo])) / np.array(box_lengths)).astype(int)
            wrapped_coords = (rotated - np.array([xlo, ylo, zlo])) % np.array(box_lengths) + np.array([xlo, ylo, zlo])

            if all_positions:
                tree = cKDTree(np.vstack(all_positions))
                distances, _ = tree.query(wrapped_coords, k=1)
                if np.any(distances < min_distance):
                    print(f"    Overlap detected. Retrying... (attempt {attempt+1})", flush=True)
                    continue

            print(f"    Success.", flush=True)
            success = True
            break


        if not success:
            raise RuntimeError(f"Failed to insert molecule {mol_id} after {max_attempts} attempts.")

        atom_map = {}
        for i, atom in enumerate(template_data['atoms']):
            old_id = atom['id']
            new_id = atom_offset + old_id
            pos = wrapped_coords[i]
            img = image_flags[i]

            new_atom = {
                'id': new_id,
                'mol': mol_id,
                'type': atom['type'],
                'charge': atom['charge'],
                'x': pos[0],
                'y': pos[1],
                'z': pos[2],
                'ix': img[0],
                'iy': img[1],
                'iz': img[2]
            }
            atom_map[old_id] = new_id
            all_atoms.append(new_atom)

        all_positions.append(wrapped_coords)

        for bond in template_data['bonds']:
            all_bonds.append({
                'id': bond_offset + bond['id'],
                'type': bond['type'],
                'atom1': atom_map[bond['atom1']],
                'atom2': atom_map[bond['atom2']]
            })

        for angle in template_data['angles']:
            all_angles.append({
                'id': angle_offset + angle['id'],
                'type': angle['type'],
                'a1': atom_map[angle['a1']],
                'a2': atom_map[angle['a2']],
                'a3': atom_map[angle['a3']]
            })

        for dih in template_data['dihedrals']:
            all_dihedrals.append({
                'id': dihedral_offset + dih['id'],
                'type': dih['type'],
                'a1': atom_map[dih['a1']],
                'a2': atom_map[dih['a2']],
                'a3': atom_map[dih['a3']],
                'a4': atom_map[dih['a4']]
            })

        for imp in template_data['impropers']:
            all_impropers.append({
                'id': improper_offset + imp['id'],
                'type': imp['type'],
                'a1': atom_map[imp['a1']],
                'a2': atom_map[imp['a2']],
                'a3': atom_map[imp['a3']],
                'a4': atom_map[imp['a4']]
            })

        atom_offset += len(template_data['atoms'])
        bond_offset += len(template_data['bonds'])
        angle_offset += len(template_data['angles'])
        dihedral_offset += len(template_data['dihedrals'])
        improper_offset += len(template_data['impropers'])

    return {
        'atoms': all_atoms,
        'bonds': all_bonds,
        'angles': all_angles,
        'dihedrals': all_dihedrals,
        'impropers': all_impropers
    }
