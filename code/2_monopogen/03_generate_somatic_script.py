"""
Script: 03_generate_somatic_script.py

Description:
Generates SLURM batch scripts to identify somatic variants using the Monopogen pipeline for single-cell data. 
It processes input directories containing sequencing data, using specified genomic resources, and outputs somatic variant analysis results.

Usage:
python 03_generate_somatic_script.py --input_dir <Monopogen_directory> 
                                     --slurm_scripts_dir <slurm_scripts_directory>
                                     --region_list_path <region_list_file>
                                     --genome_fa_path <genome_fasta_file>
                                     --threads <num_threads>

Arguments:
--input_dir          Directory containing input data for somatic variant analysis.
--slurm_scripts_dir  Monopogen Directory to save generated SLURM scripts for somatic variant analysis.
--genome_fa_path     Path to the Homo_sapiens.GRCh38.dna.primary_assembly.modified.fa generated in STARSolo
--threads            Number of threads (CPUs) to allocate for processing.

Example:
python 03_generate_somatic_script.py --input_dir /path/to/input_data 
                                     --slurm_scripts_dir /path/to/slurm_scripts
                                     --genome_fa_path /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.modified.fa
                                     --threads 64
"""

import argparse
from pathlib import Path

def generate_somatic_scripts(input_dir, output_dir, slurm_scripts_dir, region_list_path, genome_fa_path, threads):
    region_list_path = Path(input_dir, "resource","GRCh38.region.somatic.lst")
    somatic_script_dir = Path(slurm_scripts_dir, "slurm_scripts", "somatic_scripts")
    somatic_script_dir.mkdir(parents=True, exist_ok=True)

    for sample_dir in Path(input_dir).iterdir():
        if sample_dir.is_dir():
            barcode_counts_csv = f"{output_dir}/{sample_dir.name}/soloOut/GeneFull/filtered/barcode_counts.csv"
            slurm_script_path = somatic_script_dir / f"{sample_dir.name}.sh"
            with open(slurm_script_path, "w") as f:
                f.write(f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={threads}
#SBATCH --time=168:00:00
#SBATCH --mem=186G
#SBATCH --job-name=sm{sample_dir.name}
#SBATCH --output={slurm_scripts_dir}/slurm_logs/somatic/{sample_dir.name}.out
#SBATCH --mail-type=FAIL,END

# Set environment variables
export PATH="{slurm_scripts_dir}/apps:$PATH"
export LD_LIBRARY_PATH=path/to/.anaconda3/envs/somaseqenv/lib:$LD_LIBRARY_PATH  # MODIFY IT: Point to conda environment lib path
export LD_LIBRARY_PATH="{slurm_scripts_dir}/apps:$LD_LIBRARY_PATH"

# Load necessary modules
module load StdEnv/2023 java/17.0.6

# index vcf files
for file in {input_dir}/{sample_dir.name}/germline/*.vcf.gz; do
    {slurm_scripts_dir}/apps/tabix -f -p vcf "$file"
done

cd {slurm_scripts_dir}

# Activate the Monopogen conda environment
conda activate somaseqenv

# Somatic variant calling steps
# featureInfo
python {slurm_scripts_dir}/src/Monopogen.py somatic \
    -a {slurm_scripts_dir}/apps \
    -r {region_list_path} \
    -t {threads} \
    -i {output_dir}/{sample_dir.name} \
    -l {barcode_counts_csv} \
    -s featureInfo \
    -g {genome_fa_path}

# cellScan
python {slurm_scripts_dir}/src/Monopogen.py somatic \
    -a {slurm_scripts_dir}/apps \
    -r {region_list_path} \
    -t {threads} \
    -w 10MB \
    -i {output_dir}/{sample_dir.name} \
    -l {barcode_counts_csv} \
    -s cellScan \
    -g {genome_fa_path}

# LDrefinement
python {slurm_scripts_dir}/src/Monopogen_opt.py somatic \
    -a {slurm_scripts_dir}/apps \
    -r {region_list_path} \
    -t {threads} \
    -i {output_dir}/{sample_dir.name} \
    -l {barcode_counts_csv} \
    -s LDrefinement \
    -g {genome_fa_path}

""")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SLURM batch scripts for Monopogen somatic analysis.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input sample directories.")
    parser.add_argument("--output_dir", required=True, help="Directory to save Monopogen output.")
    parser.add_argument("--slurm_scripts_dir", required=True, help="Directory to store SLURM scripts.")
    parser.add_argument("--genome_fa", required=True, help="Path to the genome fasta file.")
    parser.add_argument("--threads", type=int, default=64, help="Number of threads (CPUs) to use.")
    
    args = parser.parse_args()
    
    generate_somatic_scripts(args.input_dir, args.output_dir, args.slurm_scripts_dir, args.region_list, args.genome_fa, args.threads)
    print(f"Somatic SLURM scripts generated in: {args.slurm_scripts_dir}")