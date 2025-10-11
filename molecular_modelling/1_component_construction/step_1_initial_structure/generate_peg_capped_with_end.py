"""
Generate a PEG molecule with end-groups and write a PDB file.
Build: [Left end] -- [PEG chain with N monomers] -- [Right end]
"""

from __future__ import annotations
import os
from typing import Dict, List, Tuple

# ============================ CONFIG (edit here) ============================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILENAME = os.path.join(SCRIPT_DIR, "peg_capped_with_end.pdb")
NUM_MONOMERS = 728  # number of PEG monomers
# ===========================================================================

Coord = Tuple[int, str, str, float, float, float]
ConMap = Dict[int, List[int]]

# Left end group
LEFT_END_COORDS: List[Coord] = [
    (1, "H0", "H", -8.149, 1.516, 2.594), (2, "H0", "H", -6.413, -0.355, 2.581),
    (3, "H0", "H", -8.836, -0.882, 2.803), (4, "C0", "C", -8.668, 1.022, 1.769),
    (5, "H0", "H", -9.73, 1.257, 1.867),  (6, "C0", "C", -8.454, -0.51, 1.85),
    (7, "C0", "C", -6.942, -0.82, 1.745), (8, "H0", "H", -6.074, 1.759, 1.091),
    (9, "H0", "H", -2.457, -0.997, 1.385), (10, "H0", "H", -6.786, -1.899, 1.831),
    (11, "C0", "C", -8.13, 1.546, 0.415),  (12, "C0", "C", -6.616, 1.233, 0.3),
    (13, "H0", "H", -8.284, 2.626, 0.355), (14, "H0", "H", -1.533, 0.962, 0.085),
    (15, "C0", "C", -2.539, -0.789, 0.315),(16, "C0", "C", -1.272, -0.047, -0.154),
    (17, "C0", "C", -3.819, -0.012, 0.064),(18, "C0", "C", -6.356, -0.303, 0.4),
    (19, "N0", "N", -4.937, -0.724, 0.352),(20, "O0", "O", -3.821, 1.141, -0.318),
    (21, "H0", "H", -4.784, -1.675, 0.649),(22, "C0", "C", -9.199, -1.206, 0.684),
    (23, "H0", "H", -10.272, -1.011, 0.758),(24, "H0", "H", -9.062, -2.289, 0.746),
    (25, "H0", "H", -6.256, 1.62, -0.657), (26, "H0", "H", -2.598, -1.754, -0.195),
    (27, "C0", "C", -8.884, 0.851, -0.746),(28, "H0", "H", -1.175, -0.692, -1.003),
    (29, "H0", "H", -9.952, 1.073, -0.684),(30, "C0", "C", -8.653, -0.678, -0.666),
    (31, "C0", "C", -7.134, -0.981, -0.767),(32, "H0", "H", -6.978, -2.063, -0.743),
    (33, "H0", "H", -8.529, 1.234, -1.707),(34, "H0", "H", -6.75, -0.626, -1.727),
    (35, "H0", "H", -9.179, -1.169, -1.488),
]
LEFT_END_CONECT: ConMap = {
    1: [4], 2: [7], 3: [6], 4: [1, 5, 11, 6], 5: [4],
    6: [3, 4, 22, 7], 7: [2, 6, 10, 18], 8: [12], 9: [15], 10: [7],
    11: [4, 13, 27, 12], 12: [8, 11, 25, 18], 13: [11], 14: [16],
    15: [9, 26, 17, 16], 16: [14, 15, 28], 17: [15, 20, 19],
    18: [7, 12, 19, 31], 19: [17, 18, 21], 20: [17], 21: [19],
    22: [6, 23, 24, 30], 23: [22], 24: [22], 25: [12], 26: [15],
    27: [11, 29, 33, 30], 28: [16], 29: [27],
    30: [22, 27, 35, 31], 31: [18, 30, 32, 34],
    32: [31], 33: [27], 34: [31], 35: [30],
}

# PEG monomers
MON1_COORDS: List[Coord] = [
    (1, "O1", "O", 0.0, 0.0, 0.0), (2, "C2", "C", 1.181, 0.002, -0.799),
    (3, "C3", "C", 2.397, -0.001, 0.11), (4, "H4", "H", 1.196, -0.883, -1.447),
    (5, "H5", "H", 1.197, 0.89, -1.442), (6, "H6", "H", 2.381, -0.89, 0.752),
    (7, "H7", "H", 2.382, 0.882, 0.759),
]
MON1_CON: ConMap = {1: [2], 2: [1, 3, 4, 5], 3: [2, 6, 7], 4: [2], 5: [2], 6: [3], 7: [3]}
MON2_COORDS: List[Coord] = [
    (1, "O1", "O", 0.0, 0.002, -0.698), (2, "C2", "C", 1.181, -0.003, 0.107),
    (3, "C3", "C", 2.397, 0.003, -0.806), (4, "H4", "H", 1.196, -0.894, 0.747),
    (5, "H5", "H", 1.197, 0.894, 0.747), (6, "H6", "H", 2.381, -0.88, -1.456),
    (7, "H7", "H", 2.382, 0.88, -1.456),
]
MON2_CON: ConMap = {1: [2], 2: [1, 3, 4, 5], 3: [2, 6, 7], 4: [2], 5: [2], 6: [3], 7: [3]}

# Right end group
RIGHT_END_COORDS: List[Coord] = [
    (1, "H0", "H", 6.188, 2.711, -1.774), (2, "H0", "H", 3.853, 2.461, -0.752),
    (3, "H0", "H", 5.886, 2.473, 0.7),    (4, "C0", "C", 6.49, 1.788, -1.273),
    (5, "C0", "C", 4.181, 1.536, -0.27),  (6, "C0", "C", 5.696, 1.617, 0.046),
    (7, "H0", "H", 7.559, 1.872, -1.058), (8, "H0", "H", 3.629, 1.452, 0.669),
    (9, "H0", "H", 4.392, 1.384, -3.032), (10, "C0", "C", 4.71, 0.482, -2.503),
    (11, "C0", "C", 6.225, 0.573, -2.194),(12, "C0", "C", 3.883, 0.317, -1.195),
    (13, "C0", "C", 1.4, 0.129, -0.736),  (14, "N0", "N", 2.459, 0.23, -1.57),
    (15, "O0", "O", 1.475, 0.089, 0.476), (16, "H0", "H", 6.781, 0.694, -3.126),
    (17, "C0", "C", 6.145, 0.315, 0.756), (18, "H0", "H", 2.269, 0.239, -2.561),
    (19, "H0", "H", 7.208, 0.373, 1.004), (20, "H0", "H", 5.601, 0.194, 1.697),
    (21, "H0", "H", 4.525, -0.362, -3.173),(22, "C0", "C", 6.68, -0.726, -1.485),
    (23, "C0", "C", 4.367, -0.987, -0.491),(24, "C0", "C", 5.882, -0.9, -0.169),
    (25, "H0", "H", 7.751, -0.681, -1.271),(26, "H0", "H", 3.813, -1.146, 0.437),
    (27, "H0", "H", 6.517, -1.588, -2.138),(28, "H0", "H", 4.175, -1.852, -1.129),
    (29, "H0", "H", 6.201, -1.816, 0.332),
]
RIGHT_END_CONECT: ConMap = {
    1: [4], 2: [5], 3: [6], 4: [1, 7, 11, 6], 5: [2, 12, 6, 8],
    6: [3, 4, 5, 17], 7: [4], 8: [5], 9: [10],
    10: [9, 21, 11, 12], 11: [4, 10, 22, 16], 12: [5, 10, 14, 23],
    13: [15, 14], 14: [12, 13, 18], 15: [13], 16: [11],
    17: [6, 19, 24, 20], 18: [14], 19: [17], 20: [17],
    21: [10], 22: [11, 25, 27, 24], 23: [12, 28, 24, 26],
    24: [17, 22, 23, 29], 25: [22], 26: [23], 27: [22], 28: [23], 29: [24],
}

def _add_conect(store: List[Tuple[int, int]], a1: int, a2: int) -> None:
    if a1 > a2:
        a1, a2 = a2, a1
    store.append((a1, a2))

def _get_monomer(i: int):
    return (MON1_COORDS, MON1_CON) if i % 2 == 0 else (MON2_COORDS, MON2_CON)

def generate_peg_with_ends(
    num_monomers: int = NUM_MONOMERS,
    output_filename: str = OUTPUT_FILENAME,
) -> None:
    # Internal defaults intentionally kept inside the function
    shift_per_monomer = (3.5, 0.0, 0.0)
    left_connect_id = 16
    right_connect_local_id = 13

    all_atoms: List[Tuple[int, str, str, float, float, float]] = []
    all_conects: List[Tuple[int, int]] = []

    # Left end
    for lid, name, elem, x0, y0, z0 in LEFT_END_COORDS:
        all_atoms.append((lid, name, elem, x0, y0, z0))
    for l1, lst in LEFT_END_CONECT.items():
        for l2 in lst:
            _add_conect(all_conects, l1, l2)

    # PEG chain
    atoms_per_monomer = 7
    chain_offset = len(LEFT_END_COORDS)
    current_offset = chain_offset
    for i in range(num_monomers):
        coords, conect = _get_monomer(i)
        dx = shift_per_monomer[0] * i
        dy = shift_per_monomer[1] * i
        dz = shift_per_monomer[2] * i

        for lid, name, elem, x0, y0, z0 in coords:
            gid = current_offset + lid
            all_atoms.append((gid, name, elem, x0 + dx, y0 + dy, z0 + dz))
        for l1, lst in conect.items():
            for l2 in lst:
                _add_conect(all_conects, current_offset + l1, current_offset + l2)

        if i > 0:
            _add_conect(all_conects, current_offset - atoms_per_monomer + 3, current_offset + 1)

        current_offset += atoms_per_monomer

    _add_conect(all_conects, left_connect_id, chain_offset + 1)

    # Last PEG C3 x
    last_offset = chain_offset + (num_monomers - 1) * atoms_per_monomer
    peg_tail_gid = last_offset + 3
    tail_x = next(X for (gid, _, _, X, _, _) in all_atoms if gid == peg_tail_gid)

    # Right end
    right_offset = current_offset
    for lid, name, elem, x0, y0, z0 in RIGHT_END_COORDS:
        gid = right_offset + lid
        all_atoms.append((gid, name, elem, x0 + tail_x, y0, z0))
    for l1, lst in RIGHT_END_CONECT.items():
        for l2 in lst:
            _add_conect(all_conects, right_offset + l1, right_offset + l2)
    _add_conect(all_conects, peg_tail_gid, right_offset + right_connect_local_id)

    # Output PDB
    conect_set = {(min(a, b), max(a, b)) for (a, b) in all_conects}
    with open(output_filename, "w") as f:
        f.write("REMARK  PEG + 2 end-groups generated\n")
        f.write(f"REMARK  num_monomers = {num_monomers}\n")
        for gid, aname, elem, X, Y, Z in sorted(all_atoms, key=lambda t: t[0]):
            f.write(
                f"HETATM{gid:5d} {aname:>4s} UNL A   1    "
                f"{X:8.3f}{Y:8.3f}{Z:8.3f}  1.00  0.00           {elem:>2s}\n"
            )
        for c1, c2 in sorted(conect_set):
            f.write(f"CONECT{c1:5d}{c2:5d}\n")
        f.write("END\n")

    print(f"Done! Wrote {output_filename} (PEG with {num_monomers} monomers).")

if __name__ == "__main__":
    generate_peg_with_ends()
