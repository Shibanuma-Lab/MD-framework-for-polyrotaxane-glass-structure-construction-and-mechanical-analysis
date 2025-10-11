# 3_multichain_construction

## Overview  
This directory contains the scripts required to **construct a periodic multi-chain simulation cell** for the polyrotaxane glass (PRg) system.  
The process replicates and randomly packs multiple single-chain PRg molecules into a defined simulation box, generating a ready-to-use LAMMPS data file for subsequent annealing and equilibration.

---

##  Purpose  
After assembling a single PRg chain in the previous stage, this step expands the model to a **bulk-like multi-chain system**.  
Each single-chain configuration is duplicated, spatially randomized, and packed into a cubic periodic box with defined boundaries and minimum intermolecular distances.

---

##  Scripts Overview  

| Script | Description |
|--------|--------------|
| **`main.py`** | Main control script that orchestrates the construction of the periodic multi-chain simulation cell. |
| **`read_lammps_data.py`** | Reads atom, bond, angle, dihedral, and improper information from a LAMMPS data file. |
| **`generate_molecules.py`** | Randomly positions and orients multiple PRg chains inside the periodic box while avoiding overlaps. |
| **`write_lammps_data.py`** | Writes the final packed configuration into a new LAMMPS-compatible data file. |

---

## 🔧 Main Script: `main.py`

### Function  
This script generates a **multi-chain periodic PRg system** from a single-chain configuration.

**Workflow:**
1. Read a single-chain LAMMPS data file (`prg_n1456.data`).  
2. Replicate and randomly place `num_chains` copies inside a cubic periodic box.  
3. Ensure no two molecules are closer than the specified `min_distance`.  
4. Write the final packed structure to a new `.data` file.

---

### 🧩 User Parameters  
Located at the top of `main.py`:

```python
# ====================== User Parameters ======================
num_chains = 10                       # number of PRg chains to generate
box_half = 196.228301691482           # half box length (Å)
min_distance = 0.1                    # minimum allowed intermolecular distance
max_attempts = 100                    # maximum trial attempts per molecule
input_filename = "prg_singlechain.data"     # input single-chain data file
output_filename = "prg_systems.data" # output multi-chain data file
# =============================================================

