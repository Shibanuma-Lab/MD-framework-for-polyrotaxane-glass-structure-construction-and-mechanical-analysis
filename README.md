## Overview

This repository provides a **comprehensive molecular dynamics (MD) framework** for the modelling, equilibration, and mechanical analysis of **polyrotaxane glass (PRG)**.  

To address this limitation, this repository provides a **unified and extensible MD workflow** that integrates:

1. Molecular component generation  
2. Single-chain assembly  
3. Multi-chain system construction  
4. Annealing and equilibration  
5. Mechanical response evaluation  

In addition, a **geometry-based molecular visualization module** is provided to intuitively and quantitatively represent structural configurations and their evolution during deformation.

---
## Repository Structure

The repository consists of two main modules, corresponding directly to the methodological components of the paper.

MD-framework-for-polyrotaxane-glass-structure-construction-and-mechanical-analysis/
├── molecular_modelling/
│ ├── 1_component_construction/
│ ├── 2_single_chain_assembly/
│ ├── 3_multichain_construction/
│ ├── 4_annealing_equilibration/
│ └── 5_mechanical_response/
│
├── geometry_visualization/
│
└── README.md

---

## Molecular Modelling Workflow

The `molecular_modelling/` directory implements the **core MD workflow**, organized explicitly according to the modelling steps described in the paper:

1. **Component construction**  
   Generation of molecular building blocks (axis polymer, cyclic molecules, end-capping groups).

2. **Single-chain assembly**  
   Construction of individual polyrotaxane chains with prescribed molecular architecture.

3. **Multi-chain system construction**  
   Replication and packing of multiple chains to generate periodic simulation cells.

4. **Annealing and equilibration**  
   Automated multi-stage annealing and equilibration, including a **density-based convergence protocol** to ensure thermodynamically stable and reproducible configurations.

5. **Mechanical response evaluation**  
   Deformation simulations and extraction of stress–strain behaviour and related mechanical quantities.

Each subdirectory contains its own scripts and documentation describing the corresponding modelling step.

---

## Geometry-Based Molecular Visualization

The `geometry_visualization/` directory provides a **geometry-based visualization framework** for PRG systems.  
Instead of atomistic rendering, molecular structures are represented using simplified but physically representative geometric objects:

- Tubes for axis polymers  
- Hollow frustums (truncated cones) for cyclic molecules  
- Spheres for end-capping groups  

This approach enables efficient visualization and quantitative analysis of molecular configurations, orientations, and structural evolution during deformation.

Detailed usage instructions and script descriptions are provided in the `geometry_visualization/README.md`.

---

## Requirements

- Python 
- LAMMPS 
- Standard scientific Python libraries (NumPy, etc.)  
- PyVista (for visualization module)

Specific dependencies and usage details are described in the corresponding subdirectories.

---

## Usage

This repository is intended to be used **step by step following the workflow order**:

1. Run scripts in `molecular_modelling/1_component_construction/`
2. Proceed sequentially through `molecular_modelling/2` to `5`
3. Use `geometry_visualization/` for post-processing and visualization of simulation results

Each directory contains a dedicated `README.md` explaining the expected inputs, outputs, and execution procedures.

---

## Citation

If you use this framework in your research, please cite:

> *Molecular dynamics modelling of polyrotaxane glass:  
> A workflow for structure construction and mechanical analysis*  

---

