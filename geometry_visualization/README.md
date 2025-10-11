# Geometry-Based Molecular Visualization

## Overview  
This directory contains a complete set of **geometry-based visualization scripts** for polyrotaxane glass (PRg) systems.  
The workflow enables users to **extract geometric information** (molecular centers and orientations) from LAMMPS data files and visualize the PRg model using **PyVista**.  

The visualization is **simplified but physically representative**, using:
- **Tubes** for axis polymer  
- **Hollow frustums (truncated cones)** for cyclic molecule  
- **Spheres** for end-capping group  

---

## 🧱 Directory Contents  

| No. | Script | Function |
|-----|---------|-----------|
| 1 | `1_separate_peg_and_cd_molecule.py` | Separates PEG and CD atoms from the full LAMMPS data file |
| 2 | `2_extract_main_chain.py` | Extracts the PEG backbone atoms and writes simplified main-chain coordinates |
| 3 | `3_extract_central_coordinate.py` | Calculates the geometric centers of CD molecules (with PBC correction) |
| 4 | `4_extract_normal_vector.py` | Computes the orientation (normal vector) of each CD molecule |
| 5 | `5_plot_pyvista.py` | Visualizes the PRg model using PyVista (chain, chain ends, CDs) |
| 6 | `6_main.py` | Integrates all modules for batch processing and automated visualization |

---
