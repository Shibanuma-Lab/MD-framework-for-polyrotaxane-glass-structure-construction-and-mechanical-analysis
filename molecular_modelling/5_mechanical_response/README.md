# 5_mechanical_response

## Overview  
This directory contains the **LAMMPS input script** used to perform **uniaxial tensile deformation** of an equilibrated polyrotaxane glass (PRg) system.  
The simulation applies a controlled strain rate along a specified direction while maintaining a constant temperature (NVT ensemble).  
This stage evaluates the mechanical response of PRg, producing stress–strain data and structural snapshots for post-analysis.

---

## 📜 Script: `uniaxial_tensile_deformation.in`

### Function  
Performs **tensile deformation** under periodic boundary conditions using LAMMPS.  
The deformation is applied incrementally via `fix deform`, and system configurations are output periodically for visualization and stress–strain analysis.

---

### ⚙️ Workflow Description  

**Main steps:**
1. **Initialization**  
   - Define units, force field styles, and GPU acceleration.  
   - Read equilibrated structure (`Step10-final.data`) and force field settings.  

2. **Define Variables**  
   - Specify the number of PRg chains, temperature, strain rate, and output frequency.  
   - Compute instantaneous engineering strain from the box deformation.  

3. **Apply Deformation**  
   - Use `fix deform` to stretch the simulation box along a chosen direction (x, y, or z).  
   - Maintain temperature control with `fix nvt`.  

4. **Loop Execution**  
   - Perform deformation in small increments (`output_int` steps).  
   - After each segment, output a data file (`tensile-segment-n.data`) for post-processing.  

5. **Thermodynamic Output**  
   - Record stress, strain, temperature, and volume at each step to analyze the tensile response.

---

## 🧩 User Parameters  

Located in the **User Parameters** section:

```LAMMPS
variable num_chains     equal 50          # number of PRg chains
variable temp_load      equal 300         # tensile test temperature (K)
variable strain_rate    equal 1.0e-5      # engineering strain rate (1/ps)
variable num_steps      equal 200000      # total steps
variable output_int     equal 500         # data output interval

