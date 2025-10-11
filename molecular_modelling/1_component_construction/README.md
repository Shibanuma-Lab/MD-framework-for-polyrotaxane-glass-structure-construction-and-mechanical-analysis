# 1_component_construction

## Overview  
This directory contains the **molecular components** used to build the polyrotaxane glass (PRg) system.  
It includes three fundamental molecular units represented as `.pdb` files:

1. **PEG monomer** – the repeat unit of the linear polymer axis.  
2. **Propionylated α-cyclodextrin (α-CD)** – cyclic molecule with a **67% modification degree**, acting as the sliding ring.  
3. **Adamantane** – terminal end-capping group that prevents CD detachment.

These components form the basis of all subsequent modeling stages in the PRg construction workflow.

---

## 🔧 Script: `generate_peg_capped_with_end.py`

### Function  
This script automatically **generates a polyethylene glycol (PEG) chain** capped with adamantane end groups and outputs the molecular structure as a `.pdb` file.

**Topology:**
[Left adamantane end] -- [PEG chain with N monomers] -- [Right adamantane end]

### Features  
- Flexible control of **chain length** via the number of PEG monomers (`num_monomers`).  
- Generates a **chemically continuous** chain with both termini capped.  
- Automatically writes **CONECT records** to describe bonding topology.

---

## 🧩 User Parameters  
Located at the top of the script:

```python
# -------------------- User parameters --------------------
num_monomers = 728                # number of PEG monomers
pdb_file = "peg_with_endcap.pdb"  # output PDB filename
# ---------------------------------------------------------
