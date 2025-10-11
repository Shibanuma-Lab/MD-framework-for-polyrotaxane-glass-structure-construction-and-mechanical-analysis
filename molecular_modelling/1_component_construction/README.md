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

