# 4_annealing_equilibration

## Overview  
This directory contains the Python-based automation script for performing **multi-stage annealing and equilibration** of the polyrotaxane glass (PRg) system.  
The procedure combines **GPU-accelerated LAMMPS simulations** with a **density-based convergence criterion**, ensuring that the system reaches a thermodynamically stable and physically realistic configuration before moving to the next stage.

---

## 🔧 Script: `annealing_equilibration.py`

### Function  
The script automates a **10-stage NPT annealing protocol** for PRg or similar polymeric systems.  
Each stage runs until the system’s density evolution satisfies a **linear-slope convergence condition**, allowing for self-termination without manual monitoring.

---

### ⚙️ Workflow Description  

**Main steps:**
1. **Initialization**  
   - Loads the input LAMMPS data file (`input_data`).  
   - Configures force fields, interaction styles, and GPU acceleration.  

2. **Annealing Stages**  
   - Sequentially executes multiple **temperature–pressure (NPT)** stages.  
   - Temperature and pressure vary according to a predefined protocol (heating, pressurization, cooling, depressurization).  

3. **Density-Based Convergence Check**  
   - At regular intervals (`check_interval`), the instantaneous system density is recorded.  
   - The slope of a linear fit to the last `history_length` points is calculated.  
   - If the slope remains below a specified tolerance (`slope_tolerance`) for `consecutive_needed` consecutive checks, the system is considered converged and the next stage begins.  

4. **Automatic Output**  
   - Each stage writes an output file `StepX-final.data` after convergence.  
   - The full log is recorded in `anneal_prg_system.lammps`.

---

## 🧩 User Parameters  

Located at the top of the script:

```python
# ========================== User Parameters ==========================
input_data = "prg_system.data"          # input LAMMPS data file
log_file = "anneal_prg_system.lammps"
max_total_steps = 1000000               # global maximum steps
check_interval = 50                     # steps per run chunk
consecutive_needed = 200                # consecutive satisfied checks required

DEFAULT_STAGE_SETTINGS = {
    "timestep": 0.5,                    # fs
    "slope_tolerance": 1.0e-3,          # density/ps
    "history_length": 401,              # sliding window length
}
# =====================================================================

