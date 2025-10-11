import re

def read_lammps_data(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    data = {
        'masses': {},
        'atoms': [],
        'bonds': [],
        'angles': [],
        'dihedrals': [],
        'impropers': [],
        'box': [],
        'counts': {
        'atom_types': 0,
        'bond_types': 0,
        'angle_types': 0,
        'dihedral_types': 0,
        'improper_types': 0
        }
    }

    section = None
    box_bounds = []

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        # Parse type counts
        if 'atom types' in line:
            data['counts']['atom_types'] = int(line.split()[0])
        elif 'bond types' in line:
            data['counts']['bond_types'] = int(line.split()[0])
        elif 'angle types' in line:
            data['counts']['angle_types'] = int(line.split()[0])
        elif 'dihedral types' in line:
            data['counts']['dihedral_types'] = int(line.split()[0])
        elif 'improper types' in line:
            data['counts']['improper_types'] = int(line.split()[0])

        # Box size
        if 'xlo xhi' in line or 'ylo yhi' in line or 'zlo zhi' in line:
            box_bounds.append(tuple(map(float, line.split()[:2])))
            continue

        # Detect sections
        if line.startswith('Masses'):
            section = 'masses'
            continue
        elif line.startswith('Atoms'):
            section = 'atoms'
            continue
        elif line.startswith('Bonds'):
            section = 'bonds'
            continue
        elif line.startswith('Angles'):
            section = 'angles'
            continue
        elif line.startswith('Dihedrals'):
            section = 'dihedrals'
            continue
        elif line.startswith('Impropers'):
            section = 'impropers'
            continue
        elif re.match(r'^\w+', line) and not line[0].isdigit():
            section = None  # reset for any other header
            continue

        # Parse data lines
        if section == 'masses':
            parts = line.split()
            atom_type = int(parts[0])
            mass = float(parts[1])
            data['masses'][atom_type] = mass

        elif section == 'atoms':
            parts = line.split()
            atom = {
                'id': int(parts[0]),
                'mol': int(parts[1]),
                'type': int(parts[2]),
                'charge': float(parts[3]),
                'x': float(parts[4]),
                'y': float(parts[5]),
                'z': float(parts[6])
            }
            data['atoms'].append(atom)

        elif section == 'bonds':
            parts = line.split()
            bond = {
                'id': int(parts[0]),
                'type': int(parts[1]),
                'atom1': int(parts[2]),
                'atom2': int(parts[3])
            }
            data['bonds'].append(bond)

        elif section == 'angles':
            parts = line.split()
            angle = {
                'id': int(parts[0]),
                'type': int(parts[1]),
                'a1': int(parts[2]),
                'a2': int(parts[3]),
                'a3': int(parts[4])
            }
            data['angles'].append(angle)

        elif section == 'dihedrals':
            parts = line.split()
            dih = {
                'id': int(parts[0]),
                'type': int(parts[1]),
                'a1': int(parts[2]),
                'a2': int(parts[3]),
                'a3': int(parts[4]),
                'a4': int(parts[5])
            }
            data['dihedrals'].append(dih)

        elif section == 'impropers':
            parts = line.split()
            imp = {
                'id': int(parts[0]),
                'type': int(parts[1]),
                'a1': int(parts[2]),
                'a2': int(parts[3]),
                'a3': int(parts[4]),
                'a4': int(parts[5])
            }
            data['impropers'].append(imp)

    # Store box info
    if len(box_bounds) == 3:
        data['box'] = [box_bounds[0][0], box_bounds[0][1],
                       box_bounds[1][0], box_bounds[1][1],
                       box_bounds[2][0], box_bounds[2][1]]
    else:
        raise ValueError("Box bounds not properly parsed.")

    return data
