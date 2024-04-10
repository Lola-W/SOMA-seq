# SOMA-seq code repository
This repo is a place holder for running the SOMA-seq pipeline and reproducing the presented analyses.

## Table of Contents

[//]: # 

* [Introduction](#introduction)
* [Installation](#installation)
* [Somatic mutation calling](#Somatic-mutation-calling)
<!-- * [Quick Start](#quick-start)
  * [Data preprocess](#data-preprocess)
  * [Germline SNV calling](#germline-snv-calling)
* [Germline SNV calling from snRNA-seq](#germline-snv-calling-from-snRNA-seq)
  * [variant calling](#variant-calling)
  * [genotyping accuracy evaluation](#genotyping-accuracy-evaluation)
  * [ancestry identification](#ancestry-identification)
* [Somatic SNV calling from scRNA-seq](#somatic-snv-calling-from-scrna-seq)
  * [preprocess](#preprocess)
  * [germline calling](#germline-calling)
  * [ld refinement on putative somatic SNVs](#ld-refinement-on-putative-somatic-SNVs)
* [FAQs](#faqs) -->
* [Citation](#citation)

[//]: # 

## Introduction
Soma-seq is a pipeline designed for the detection and analysis of single-nucleus RNA sequencing (snRNA-seq) data to investigate somatic mutations across diverse cell types and diseases. 
<image src="./resources/OverviewPipeline.png" width="600"> 

The SOMA-seq pipeline is developed using `Python version 3.10.8`
`R version 4.3.1 (2023-06-16)`, `Platform: x86_64-pc-linux-gnu`and `aarch64-apple-darwin22.4.0, and `Running under: `x86_64, linux-gnu`, `aarch64`, and `darwin22.4.0`. Please note that the code are provided and optimized for SLURM users with HPC.

The pipeline consists of three main stages:

1. Somatic mutation calling: Utilizing Monopogen, a state-of-the-art SNV calling package developed by [Ken chen's lab](https://www.mdanderson.org/research/departments-labs-institutes/labs/ken-chen-laboratory.html), Soma-seq accurately identifies somatic mutations from snRNA-seq data. Monopogen is specifically designed to handle the unique challenges of single-cell sequencing datasets, such as sparsity and allelic dropout, ensuring high-confidence mutation calls.

2. Differential mutation analysis: Soma-seq performs a comprehensive analysis of the consequences of differential mutations across various cell types and disease groups. By employing advanced statistical methods and mixed-effect modeling, the pipeline accounts for potential confounding factors, such as age, read depth, and donor-specific effects, to uncover cell type-specific mutation patterns and rates.

3. Pathogenicity analysis: Leveraging the power of [AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492) predictions, Soma-seq assesses the pathogenic potential of identified somatic mutations. This stage provides valuable insights into the functional impact of mutations, aiding in the identification of potentially deleterious variants that may contribute to disease pathogenesis.

The output of Soma-seq includes raw csv summary of putative mutations, summaries of single base substitution patterns, mutation abundance across cell types, ranked enriched gene lists, pathway enrichment results, and a list of deleterious somatic mutations annotated by AlphaMissense. These comprehensive results enable researchers to gain a deeper understanding of the role of somatic mutations in various diseases and cell types, facilitating the discovery of potential disease biomarkers and the inference of underlying causal factors.

## Installation
SOMA-seq is avaiable on github
```
git clone https://github.com/Lola-W/SOMA-seq.git
cd SOMA-seq
```
Setup the conda environment after cloning this repo:

**Note:** Please modify the `prefix` field to match your local conda installation path in `environment.yml` with your desired conda environments directory path. For example:
prefix: `/{yourUsername}/.anaconda3/envs/somaseqenv`
```
conda env create -f environment.yml
conda activate somaseqenv
```
Then install Monopogen:
```
git clone https://github.com/KChen-lab/Monopogen.git
cd Monopogen 
pip install -e .
```

## Somatic mutation calling ##
Before starting with the analysis, ensure all dependencies are installed as listed in `environment.yml` and all necessary data is available in the data directory.

### Aligment
* **For 10x based snRNA seq data**, follow the steps below

First download the required files in [resources_to_download.txt](resources/resources_to_download.txt)

It is optional but you could set global variables like `export INPUTFOLDER="/path/to/input_data"`

#### 0. **Generate Genome Directory**: 

  The first step is preparing the GenomeDir for alignment. Run the [`00_generate_genomeDir.sh`](code/alignment/00_generate_genomeDir.sh) script in `code/alignment` to modify FASTA headers to match the GENCODE format and to generate a STAR genome directory.

  ```bash
  ./code/alignment/00_generate_genomeDir.sh /path/to/Genomic_references/hg38 /path/to/logs
  ```

  Replace `/path/to/Genomic_references/hg38` with the absolute path to where your FASTA and GTF files are located, and `/path/to/logs` with the path where you want the logs to be saved.

  This script requires the STAR aligner to be installed or you have a module for STAR 2.7.11a (see the source code).

#### 1. **sbatch Script Preparation**:

  Then we prepare SLURM batch scripts for processing each RNA-seq data sample with the STAR aligner. 

  Before running [`01_generate_sbatch_script.py`](code/alignment/01_generate_sbatch_script.py), ensure you have:
  - Modified the `--genomeDir` in the generated scripts to point to your $BUILD_PATH/USE_THIS_GenomeDir created in step 0.
  - Downloaded `3M-february-2018.txt` whitelist file as instructed in `resources_to_download.txt` and updated the `--soloCBwhitelist` path in the script accordingly.

  To generate sbatch scripts:

  ```bash
  python code/alignment/01_generate_sbatch_script.py --input_dir /path/to/input_samples --output_dir /path/to/alignment_output --slurm_scripts_dir /path/to/generated_slurm_scripts
  ```

  Modify `/path/to/input_samples`, `/path/to/alignment_output`, and `/path/to/generated_slurm_scripts` with your specific directories. The `input_dir` should contain directories for each of your samples that you wish to align.

#### 2. **Batch Submission of sbatch Scripts**:

  To Submit these scripts as batch jobs to the SLURM scheduler. The script `02_submit_slurm_jobs.py` automates the submission of uncompleted jobs by checking which samples haven't been processed based on the absence of their SLURM log files.

  To submit the jobs, run:

  ```bash
  python code/alignment/02_submit_slurm_jobs.py --res_folder /path/to/slurm_logs --input_folder /path/to/input_data --slurm_scripts_dir /path/to/slurm_scripts
  ```
