"""
Script: 01_generate_sbatch_script.py

Description:
Generates SLURM batch scripts for processing RNA-seq data with STAR aligner. 
This script automates the creation of batch scripts for each sample directory found in the input directory.
It allows users to specify custom input and output directories, as well as the directory to store the generated SLURM scripts.

Users should:
1. Change the --genomeDir to point to the output of `00_generate_genomeDir.sh`.
2. Download and unzip the 3M-february-2018.txt whitelist file. Change the --soloCBwhitelist to its path

Usage:
python 01_generate_sbatch_script.py --input_dir <path_to_input_dir> --output_dir <path_to_output_dir> --slurm_scripts_dir <path_to_slurm_scripts_dir>

Arguments:
--input_dir: Directory containing sample directories to process, must contain {sample_dir.name}/{sample_dir.name}_R2_001.fastq.gz and {sample_dir.name}/{sample_dir.name}_R1_001.fastq.gz
--output_dir: Directory where processed outputs will be stored.
--slurm_scripts_dir: Directory where generated SLURM scripts will be saved.

Example:
python 01_generate_sbatch_script.py --input_dir /path/to/input --output_dir /path/to/output --slurm_scripts_dir /path/to/slurm_scripts
"""

import argparse
from pathlib import Path
            
def generate_sbatch_scripts(input_dir, output_dir, slurm_scripts_dir):
    """
    Generates SLURM batch scripts for each sample directory in the input directory.
    """
    slurm_scripts_dir_path = Path(slurm_scripts_dir, "slurm_scripts")
    slurm_scripts_dir_path.mkdir(parents=True, exist_ok=True)

    for sample_dir in Path(input_dir).iterdir():
        if sample_dir.is_dir():
            script_path = slurm_scripts_dir_path / f"{sample_dir.name}.sh"
            with open(script_path, "w") as f:
                f.write(f"""#!/bin/bash
#SBATCH --job-name={sample_dir.name}
#SBATCH --output={slurm_scripts_dir}/slurm_logs/{sample_dir.name}.out
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --mail-type=FAIL
source ~/.bashrc
cd {output_dir}/{sample_dir.name}
mkdir -p {output_dir}/{sample_dir.name}
STAR --runThreadN 80 \\
--genomeDir /path/to/Genomic_references/hg38/USE_THIS_GenomeDir \\ # MODIFY IT: Point to your GenomeDir
--readFilesCommand zcat \\
--soloType CB_UMI_Simple \\
--soloUMIdedup 1MM_CR \\
--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \\
--soloCBwhitelist ./resource/3M-february-2018.txt \\ # MODIFY IT: Download and specify the path to 3M-february-2018.txt
--soloUMIfiltering MultiGeneUMI_CR \\
--clipAdapterType CellRanger4 \\
--outFilterScoreMin 30 \\
--readFilesIn {input_dir}/{sample_dir.name}/{sample_dir.name}_R2_001.fastq.gz {input_dir}/{sample_dir.name}/{sample_dir.name}_R1_001.fastq.gz \\
--soloOutFileNames soloOut/ features.tsv barcodes.tsv matrix.mtx \\
--soloStrand Unstranded \\
--soloFeatures Gene SJ GeneFull \\
--outSAMtype BAM SortedByCoordinate \\
--outSAMattributes NH HI nM AS CR UR CB CY UY UB GX GN sS sQ sM 
module load samtools/1.18
samtools index -@ 80 {output_dir}/{sample_dir.name}/Aligned.sortedByCoord.out.bam
                """)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates SLURM batch scripts for RNA-seq data processing with STAR.")
    parser.add_argument("--input_dir", required=True, help="Directory containing sample directories to process")
    parser.add_argument("--output_dir", required=True, help="Directory where processed outputs will be stored")
    parser.add_argument("--slurm_scripts_dir", required=True, help="Directory where generated SLURM scripts will be saved")

    args = parser.parse_args()

    generate_sbatch_scripts(args.input_dir, args.output_dir, args.slurm_scripts_dir)