# Molecular Modelling Workflow of Polyrotaxane Glass (PRg)

## Overview  
This repository provides a **complete molecular modelling workflow** for constructing, equilibrating, and mechanically testing polyrotaxane glass (PRg) systems using **Python** and **LAMMPS**.  
Each subfolder corresponds to a major stage in the workflow, from molecular component generation to mechanical deformation simulation.

---

## 🧬 Workflow Summary  

| Stage | Directory | Description |
|--------|------------|-------------|
| **1** | [`1_component_construction`](../1_component_construction) | Build the fundamental molecular components — PEG monomer, propionylated α-cyclodextrin (CD), and adamantane end-caps. |
| **2** | [`2_single_chain_assembly`](../2_single_chain_assembly) | Assemble a single polyrotaxane chain by threading multiple CDs onto a PEG axis with end-capping. |
| **3** | [`3_multichain_construction`](../3_multichain_construction) | Replicate and pack multiple PRg chains into a periodic simulation box to form a bulk-like structure. |
| **4** | [`4_annealing_equilibration`](../4_annealing_equilibration) | Perform automated multi-stage annealing and equilibration using a density-based convergence criterion. |
| **5** | [`5_mechanical_response`](../5_mechanical_response) | Conduct uniaxial tensile deformation to evaluate the mechanical response of the equilibrated PRg system. |

---
