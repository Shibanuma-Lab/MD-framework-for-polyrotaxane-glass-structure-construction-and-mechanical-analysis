# 2_single_chain_assembly

## Overview  
This directory provides the **assembly script for building a single polyrotaxane (PRg) chain** from its molecular components.  
It combines the previously generated **PEG chain** and **α-cyclodextrin (CD) rings**, positioning CDs evenly along the polymer axis to form a threaded structure.

---

## 🔧 Script: `define_ring_position.py`

### Function  
This script constructs a **LAMMPS input file (`combine_prg.in`)** to assemble a PEG chain with multiple CD molecules.  
It performs the following operations:

1. **Extracts two anchor positions** (x-coordinates) from a PEG data file — one near the front and one near the back of the chain.  
2. **Places a specified number of CD molecules** (`num_cd`) evenly between these anchors along the x-axis (with y = z = 0).  
3. **Writes a complete LAMMPS input file**, which loads the PEG structure and appends CDs at those anchor points.

---

### 🧩 User Parameters  
Located at the top of the script:

```python
# -------------------- User parameters --------------------
front_monomer = 5                 # PEG monomer index from the front (1-based)
back_monomer = 6                  # PEG monomer index from the back  (1-based)
num_cd = 91                       # number of CD molecules to insert
input_filename = "peg_n728.data"  # LAMMPS data file containing the PEG
# ---------------------------------------------------------

