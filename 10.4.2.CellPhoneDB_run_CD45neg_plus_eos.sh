#!/bin/bash
#SBATCH --job-name=cellphonedb
#SBATCH --time=100:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4

# Exit immediately if a command fails or an undefined variable is used
set -euo pipefail

module load miniforge3
source activate cpdb

# Project paths - edit these if you run the script from a different project location
export PROJECT_DIR="/home/khandl/data/Neonatal_eosinophils"
export CPDB_DIR="/home/khandl/data/cpdb/v5.0.0"
export CPDB_FILE="${CPDB_DIR}/cellphonedb.zip"
export TABLES_DIR="${PROJECT_DIR}/results/10.signaling_pathways/tables/CellPhoneDB"
export OUTPUT_DIR="${TABLES_DIR}/neo_eos_CD45neg_output"
export COUNTS_FILE="${TABLES_DIR}/neo_count.txt"
export META_FILE="${TABLES_DIR}/neo_meta.txt"

python3 << 'PYEOF' # runs everything in python until it says PYEOF 
import os
import sys

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# Read paths from the environment variables set above
cpdb_file_path = os.environ["CPDB_FILE"]
counts_file_path = os.environ["COUNTS_FILE"]
meta_file_path = os.environ["META_FILE"]
output_path = os.environ["OUTPUT_DIR"]

# Check that all required input files exist
for path, label in [(cpdb_file_path, "CellPhoneDB database"),
                    (counts_file_path, "counts file"),
                    (meta_file_path, "meta file")]:
    if not os.path.exists(path):
        sys.exit(f"ERROR: {label} not found: {path}")

# Create the output folder if it does not exist
os.makedirs(output_path, exist_ok=True)

# Print the settings so you can see what is being used
print("Running CellPhoneDB analysis")
print(f"  cpdb_file_path   = {cpdb_file_path}")
print(f"  counts_file_path = {counts_file_path}")
print(f"  meta_file_path   = {meta_file_path}")
print(f"  output_path      = {output_path}")

# Run the CellPhoneDB statistical analysis
cpdb_statistical_analysis_method.call(
    cpdb_file_path=cpdb_file_path,
    meta_file_path=meta_file_path,
    counts_file_path=counts_file_path,
    counts_data="ensembl",
    score_interactions=True,
    threshold=0.1,
    output_path=output_path,
)

print(f"\nFinished. Results written to: {output_path}")
PYEOF
