#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=2:00:00
#SBATCH --job-name generateGenomeDir
#SBATCH --output=/scratch/s/shreejoy/wengjiam/StarSolo/Output/generateGenomeDir_all.out # PLEASE Custom output file name,
#SBATCH --mail-type=ALL

# Sourced from https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38_2020A
# Description: This script prepares the genomic data for alignment by modifying FASTA headers
# to match GENCODE format and generating a STAR genome directory. It is the first step in the
# SOMA-seq alignment process.
#
# Prerequisites:
# - STAR aligner installed and available in PATH, or module available
# - Correctly formatted FASTA and GTF files in the specified build directory
#
# Usage:
# ./00_generate_genomeDir.sh [path_to_build_directory] [path_to_log_directory]
#   - path_to_build_directory: Absolute path where FASTA and GTF files are located.
#   - path_to_log_directory: Absolute path where output logs should be saved.
#
# Example:
# ./code/alignment/00_generate_genomeDir.sh /path/to/Genomic_references/hg38 /path/to/logs



# Module Environment Setup
# Load necessary modules and software environments. This may vary based on your HPC environment.
source ~/.bashrc
module load CCEnv
module load StdEnv/2023
module load star/2.7.11a

# User-defined variables for input and output directories
# Allows users to specify the base path for FASTA and GTF files, and the output log path via command line arguments
BUILD_PATH="/scratch/s/shreejoy/nxu/Genomic_references/hg38" # Default path to where FASTA and GTF are located. Replace with your path or pass as an argument.
LOG_PATH="/path/to/log" # Default path for logs. Replace with your path or pass as an argument.

# Check for command line arguments and override default paths if provided
if [ ! -z "$1" ]; then
  BUILD_PATH="$1"
fi

if [ ! -z "$2" ]; then
  LOG_PATH="$2"
fi

echo "Using build path: $BUILD_PATH"
echo "Log files will be stored in: $LOG_PATH"

# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
# Processing Ensembl FASTA file to match GENCODE format
# This step modifies the sequence headers in the FASTA file to ensure compatibility with downstream tools.
fasta_in="$BUILD_PATH/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fasta_modified="$BUILD_PATH/Homo_sapiens.GRCh38.dna.primary_assembly.modified.fa"
echo "Modifying FASTA headers to match GENCODE format..."

# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"

# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
gtf_in="$BUILD_PATH/gencode.v44.primary_assembly.annotation.gtf" # Note that I'm using v44 here, could modify as needed
gtf_modified="$BUILD_PATH/gencode.v44.primary_assembly.annotation.modified.gtf"
echo "Removing version suffix from GTF IDs..."

# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.
BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""

genome="GRCh38"
version="release-110"
# Genome Generation with STAR
echo "Generating genome with STAR..."
STAR \
    --runThreadN 80 \
    --runMode genomeGenerate \
GPT
To refine the provided script for publication on GitHub with enhanced documentation and the ability to specify file paths for FASTA and GTF locations, as well as the output log path through the command line, consider the following revised version of the script. This revision includes detailed comments to explain each step, making it more accessible for users with varying levels of expertise. It also incorporates parameterization for key paths to increase flexibility and ease of use.

bash
Copy code
#!/bin/bash

# SLURM Directives
# Provides SLURM with instructions on how to handle the job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=2:00:00
#SBATCH --job-name=generateGenomeDir
#SBATCH --output=%x_%j.out # Custom output file name, '%x' is job name, '%j' is job id
#SBATCH --mail-type=ALL # Notify user by email when the job starts, ends, or fails

# Module Environment Setup
# Load necessary modules and software environments. This may vary based on your HPC environment.
source ~/.bashrc
module load CCEnv
module load StdEnv/2023
module load star/2.7.11a

# User-defined variables for input and output directories
# Allows users to specify the base path for FASTA and GTF files, and the output log path via command line arguments
BUILD_PATH="/path/to/Genomic_references/hg38" # Default path to where FASTA and GTF are located. Replace with your path or pass as an argument.
LOG_PATH="/path/to/log" # Default path for logs. Replace with your path or pass as an argument.

# Check for command line arguments and override default paths if provided
if [ ! -z "$1" ]; then
  BUILD_PATH="$1"
fi

if [ ! -z "$2" ]; then
  LOG_PATH="$2"
fi

echo "Using build path: $BUILD_PATH"
echo "Log files will be stored in: $LOG_PATH"

# Processing Ensembl FASTA file to match GENCODE format
# This step modifies the sequence headers in the FASTA file to ensure compatibility with downstream tools.
fasta_in="$BUILD_PATH/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fasta_modified="$BUILD_PATH/Homo_sapiens.GRCh38.dna.primary_assembly.modified.fa"

echo "Modifying FASTA headers to match GENCODE format..."
sed -E 's/^>(\S+).*/>\1 \1/' $fasta_in | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' | sed -E 's/^>MT />chrM /' > $fasta_modified

# Adjusting GTF File to Remove Version Suffix from Identifiers
gtf_in="$BUILD_PATH/gencode.v44.primary_assembly.annotation.gtf"
gtf_modified="$BUILD_PATH/gencode.v44.primary_assembly.annotation.modified.gtf"

echo "Removing version suffix from GTF IDs..."
sed -E 's/gene_id "([^"]+)\.([0-9]+)";/gene_id "\1"; gene_version "\2";/' $gtf_in > $gtf_modified

# Genome Generation with STAR
echo "Generating genome with STAR..."
STAR --runThreadN 80 \
     --runMode genomeGenerate \
     --genomeDir $BUILD_PATH/USE_THIS_GenomeDir \
     --genomeFastaFiles $fasta_modified \
     --sjdbGTFfile $gtf_modified \
     --outFileNamePrefix $LOG_PATH/ \
    --clipAdapterType CellRanger4 --outFilterScoreMin 30
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR  

echo "Genome generation complete. Outputs are in $BUILD_PATH/USE_THIS_GenomeDir, log located in $LOG_PATH"