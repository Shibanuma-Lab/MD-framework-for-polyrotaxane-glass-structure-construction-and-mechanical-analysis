"""
Script: Automated visualization pipeline for PRg model
-------------------------------------------------------
This script sequentially executes all visualization modules to
process a single PRg data file and generate a 3D representation
of the system.

Modules executed:
    1. Separate chain and ring molecules
    2. Extract backbone atoms
    3. Extract ring centroid coordinates
    4. Extract ring orientation (normal vectors)
    5. Visualize PRg model using PyVista

Input format:
    - specific_model = "prg_nXXX"  (without the ".data" extension)

"""

import os
import subprocess
import sys
import time

# ====================== User Parameters ======================
specific_model = "prg_nXXX"  # target model name (without ".data")
# =============================================================


# ------------------- Working Directory -------------------
current_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_dir)

# ------------------- Required Files -------------------
required_files = [
    "1_separate_chain_and_ring.py",
    "2_extract_backbone.py",
    "3_extract_coordinate.py",
    "4_extract_normal_vector.py",
    "5_plot.py",
    f"{specific_model}.data",
]


# ------------------- Utility Functions -------------------
def check_files(files):
    """Verify that all required files exist before execution."""
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        print("\nMissing required files:")
        for f in missing:
            print(f"  - {f}")
        sys.exit(1)
    print("✅ All required files are present.\n")


def run_command(command):
    """Execute a shell command with error handling."""
    print(f"▶ Running: {command}")
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print(f"❌ Command failed: {command}")
        sys.exit(1)
    print(f"✅ Completed: {command}\n")


# ------------------- Execution Pipeline -------------------
def main():
    check_files(required_files)
    start_time = time.time()

    model_number = specific_model.split("_")[-1]

    # 1. Separate PEG chains and CD rings
    run_command(
        f"python3 1_separate_chain_and_ring.py {specific_model} "
        f"1_peg_{model_number}.data 1_cd_{model_number}.data"
    )

    # 2. Extract main backbone atoms
    run_command(
        f"python3 2_extract_backbone.py "
        f"1_peg_{model_number}.data 2_backbone_{model_number}.data"
    )

    # 3. Extract CD centroids
    run_command(
        f"python3 3_extract_coordinate.py "
        f"1_cd_{model_number}.data 3_cd_centers_{model_number}.data"
    )

    # 4. Extract CD orientation vectors
    run_command(
        f"python3 4_extract_normal_vector.py "
        f"1_cd_{model_number}.data 4_cd_normals_{model_number}.data"
    )

    # 5. Visualize with PyVista
    run_command(
        f"python3 5_plot.py "
        f"2_backbone_{model_number}.data 3_cd_centers_{model_number}.data "
        f"4_cd_normals_{model_number}.data PRg_visualization_{model_number}.png"
    )

    total_time = time.time() - start_time
    print(f"\n🎉 Visualization completed successfully! Total time: {total_time:.2f} s")


if __name__ == "__main__":
    main()
