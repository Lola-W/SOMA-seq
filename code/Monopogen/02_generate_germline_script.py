"""
Script: 02_generate_germline_script.py

Description:
This script generates SLURM batch scripts for processing RNA-seq data to identify germline variants
using the Monopogen pipeline. It creates a .lst file for each sample containing paths to BAM files and
a corresponding SLURM script to run the Monopogen analysis.

Usage:
python 02_generate_germline_script.py --input_dir <input_directory> --output_dir <output_directory> 
                                      --slurm_scripts_dir <slurm_scripts_directory> --threads <num_threads>

Arguments:
--input_dir           Directory containing STAR alignment results (BAM files).
--output_dir          Directory where Monopogen results should be saved.
--slurm_scripts_dir   Directory where SLURM scripts will be stored.
--threads             Number of threads (CPUs) to use for processing.

Example:
python 02_generate_germline_script.py --input_dir /path/to/STAR/results --output_dir /path/to/Monopogen/output 
                                      --slurm_scripts_dir /path/to/Monopogen --threads 80
"""

import argparse
from pathlib import Path

def generate_germline_scripts(input_dir, output_dir, slurm_scripts_dir, threads):
    bamlsts_dir = Path(output_dir, "bamlst")
    germline_script_dir = Path(slurm_scripts_dir, "slurm_scripts", "germline_scripts")
    germline_script_dir.mkdir(parents=True, exist_ok=True)

    for sample_dir in Path(input_dir).iterdir():
        if sample_dir.is_dir():
            bam_list_file = bamlsts_dir / f"{sample_dir.name}.bam.lst"
            with open(bam_list_file, "w") as f:
                f.write(f"{sample_dir.name}, {sample_dir}/Aligned.sortedByCoord.out.bam")

            slurm_script_file = germline_script_dir / f"{sample_dir.name}.sh"
            with open(slurm_script_file, "w") as f:
                f.write(f"""#!/bin/bash
#SBATCH --job-name=Germ{sample_dir.name}
#SBATCH --output={slurm_scripts_dir}/slurm_logs/germline/{sample_dir.name}.out
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={threads}
#SBATCH --mail-type=FAIL,END

# Set environment variables
export PATH="{slurm_scripts_dir}/apps:$PATH"
export LD_LIBRARY_PATH="{slurm_scripts_dir}/apps:$LD_LIBRARY_PATH"

mkdir -p {output_dir}/{sample_dir.name}

# Activate the Monopogen conda environment
conda activate somaseqenv

cd {slurm_scripts_dir}
./src/Monopogen.py preProcess \\
    -b {bam_list_file} \\
    -a apps \\
    -t {threads} \\
    -o {output_dir}/{sample_dir.name}

./src/Monopogen.py germline \\
    -a apps \\
    -t {threads} \\
    -r {slurm_scripts_dir}/resource/GRCh38.region.somatic.lst \\
    -p /path/to/1KG3_imputation_panel/ \\ # MODIFY IT: Point to path to 1KG3 imputation panel
    -g {slurm_scripts_dir}/resource/Homo_sapiens.GRCh38.dna.primary_assembly.modified.fa \\
    -s all \\
    -o {output_dir}/{sample_dir.name}
                """)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SLURM batch scripts for Monopogen germline analysis.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input BAM files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save Monopogen output.")
    parser.add_argument("--slurm_scripts_dir", required=True, help="Directory to Monopogen, where you will store SLURM scripts.")
    parser.add_argument("--threads", type=int, default=80, help="Number of threads (CPUs) to use.")
    
    args = parser.parse_args()
    
    generate_germline_scripts(args.input_dir, args.output_dir, args.slurm_scripts_dir, args.threads)
    print(f"Germline SLURM scripts generated in: {args.slurm_scripts_dir}")
