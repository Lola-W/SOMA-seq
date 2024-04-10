"""
Script to Submit Uncompleted SLURM Jobs for snRNA-seq Data Processing

Description:
This script identifies snRNA-seq samples that have not yet been processed based on the absence of corresponding SLURM log files.
It then submits SLURM jobs for these uncompleted samples using the generated SLURM scripts.

Usage:
python submit_uncompleted_jobs.py --res_folder <path_to_slurm_log_dir> --input_folder <path_to_input_data_dir> --slurm_scripts_dir <path_to_slurm_scripts_dir>

Arguments:
--res_folder: Directory where SLURM log files are stored. Used to identify completed jobs.
--input_folder: Directory containing snRNA-seq sample directories to be processed.
--slurm_scripts_dir: Directory where SLURM scripts for processing each sample are stored.

Example:
python submit_uncompleted_jobs.py --res_folder /path/to/slurm_logs --input_folder /path/to/input_data --slurm_scripts_dir /path/to/slurm_scripts
"""

import argparse
from pathlib import Path
import subprocess

def submit_jobs(res_folder, input_folder, slurm_scripts_dir):
    completed_list = [path.stem for path in Path(res_folder).glob("*L003.out")] # Please modify it if your file has a different structure
    uncompleted_list = [directory for directory in Path(input_folder).iterdir() if directory.name not in completed_list and directory.is_dir()]
    
    print(f"Number of uncompleted jobs: {len(uncompleted_list)}")
    
    for directory in uncompleted_list:
        script_path = Path(slurm_scripts_dir, f"{directory.name}.sh")
        print(f"Submitting job for {directory.name}: sbatch {script_path}")
        # Comment the following line to not actually submit the jobs
        subprocess.run(["sbatch", str(script_path)])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submit uncompleted SLURM jobs for snRNA-seq data processing.")
    parser.add_argument("--res_folder", required=True, help="Directory where SLURM log files are stored.")
    parser.add_argument("--input_folder", required=True, help="Directory containing snRNA-seq sample directories.")
    parser.add_argument("--slurm_scripts_dir", required=True, help="Directory where SLURM scripts are stored.")
    
    args = parser.parse_args()
    
    submit_jobs(args.res_folder, args.input_folder, args.slurm_scripts_dir)