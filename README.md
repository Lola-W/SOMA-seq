# SOMA-seq code repository
This repo is a place holder for running the SOMA-seq pipeline and reproducing the presented analyses.

## Table of Contents

[//]: # 

* [Introduction](#introduction)
* [Installation](#installation)
* [Somatic mutation calling](#somatic-mutation-calling)
  * [Alignment](#aligment)
  * [Monopogen](#monopogen)
  * [Preprocess for Analysis](#preprocess-for-analysis)
* [Differential mutation analysis](#differential-mutation-analysis)
  * [Mixed-effects modeling](#mixed-effects-modeling-of-somatic-snv-burden)
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
* [Supplementary Scripts](#supplementary-scripts)
* [Citation](#citation)

[//]: # 

## Introduction
Soma-seq is a pipeline designed for the detection and analysis of single-nucleus RNA sequencing (snRNA-seq) data to investigate somatic mutations across diverse cell types and diseases. 
<image src="./resources/OverviewPipeline.png" width="600"> 

The SOMA-seq pipeline is developed using `Python version 3.10.8`
`R version 4.3.1 (2023-06-16)`, `Platform: x86_64-pc-linux-gnu`and `aarch64-apple-darwin22.4.0`, and Running under: `x86_64, linux-gnu`, `aarch64`, and `darwin22.4.0`. Please note that the code are provided and optimized for SLURM users with HPC.

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
cd ..
```

## Somatic mutation calling ##
Before starting with the analysis, ensure all dependencies are installed as listed in `environment.yml` and all necessary data is available in the data directory.

### Aligment ###

* **For 10x based snRNA seq data**, follow the steps below for using **STARSolo** to algn. You could also opt for CellRanger4 or other alignment tools.

First download the required files in [resources_to_download.txt](resources/resources_to_download.txt)

It is optional but you could set global variables like `export INPUTFOLDER="/path/to/input_data"`

#### 0. **Generate Genome Directory**: 

  The first step is preparing the GenomeDir for alignment. Run the [`00_generate_genomeDir.sh`](code/1_alignment/00_generate_genomeDir.sh) script in `code/alignment` to modify FASTA headers to match the GENCODE format and to generate a STAR genome directory.

  ```bash
  ./code/1_alignment/00_generate_genomeDir.sh /path/to/Genomic_references/hg38 /path/to/logs
  ```

  Replace `/path/to/Genomic_references/hg38` with the absolute path to where your FASTA and GTF files are located, and `/path/to/logs` with the path where you want the logs to be saved.

  This script requires the STAR aligner to be installed or you have a module for STAR 2.7.11a (see the source code).

#### 1. **sbatch Script Preparation**:

  Then we prepare SLURM batch scripts for processing each RNA-seq data sample with the STAR aligner. 

  Before running [`01_generate_STARSolo_script.py`](code/1_alignment/01_generate_STARSolo_script.py), ensure you have:

  - Modified the `--genomeDir` in the [`01_generate_STARSolo_script.py`](code/1_alignment/01_generate_STARSolo_script.py) point to your `$BUILD_PATH/USE_THIS_GenomeDir` created in step 0.

  - Downloaded `3M-february-2018.txt` whitelist file as instructed in `resources_to_download.txt` and updated the `--soloCBwhitelist` path in the script accordingly.

  To generate sbatch scripts:

  ```bash
  python code/1_alignment/01_generate_STARSolo_script.py --input_dir /path/to/input_samples --output_dir /path/to/alignment_output --slurm_scripts_dir /path/to/generated_slurm_scripts
  ```

  Modify `/path/to/input_samples`, `/path/to/alignment_output`, and `/path/to/generated_slurm_scripts` with your specific directories. The `input_dir` should contain directories for each of your samples that you wish to align.

#### **Batch Submission of sbatch Scripts**:

  To Submit these scripts as batch jobs to the SLURM scheduler. The script [`submit_slurm_jobs.py`](code/submit_slurm_jobs.py) automates the submission of uncompleted jobs by checking which samples haven't been processed based on the absence of their SLURM log files.

  To submit the jobs, run:

  ```bash
  python code/submit_slurm_jobs.py --res_folder /path/to/slurm_logs --input_folder /path/to/input_data --slurm_scripts_dir /path/to/slurm_scripts
  ```

* It is optional, but if you do NOT have the barcode available after running STARsolo, please use [`generateBarcode.py`](utils/generateBarcode.py)

### Monopogen ###

Please reference the github repo of [Monopogen](https://github.com/KChen-lab/Monopogen) for trouble-shooting.

Before running Monopogen, ensure you have:

  - The conda envirment `somaseqenv` ready. Monopogen and its dependencies installed.

  - Downloaded `1KG3 reference panel` and save into a folder as instructed in `resources_to_download.txt`.

#### 2. **Preprocess and Germline Mutation Calling**: 

First we prepare `.bam.lst` files that each contain one BAM file to be processed using the [`bam_list.sh`](code/2_monopogen/bam_list.sh).

```
path="XXX/Monopogen"  # where Monopogen is downloaded
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps
./code/2_monopogen/bam_list.sh /path/to/STARresult /path/to/Monopogen/bamlst
```

Then use the [`02_generate_germline_script.py`](code/2_monopogen/02_generate_germline_script.py) to create SLURM scripts for each sample. These scripts will preprocess the data and call germline mutations using Monopogen.

Before running, ensure you have:
- change the **`-p /path/to/1KG3_imputation_panel/`** pointing to your folder to 1KG3 reference panel. 

```bash
python code/2_monopogen/02_generate_germline_script.py --input_dir /path/to/Monopogen/bamlst --output_dir /path/to/Monopogen/germline_output --slurm_scripts_dir /path/to/slurm_scripts --threads 80
```

Modify `/path/to/Monopogen/bamlst`, `/path/to/Monopogen/germline_output`, and `/path/to/slurm_scripts` with your specific directories. Note that you could  adjust `--threads 80` with the number of threads you wish to allocate per job (adjust based on number of CPU you have).

Batch submit the SLURM scripts as [before](#batch-submission-of-sbatch-scripts)

#### 3. **Somatic Mutation Calling**: 

To call somatic mutation based on the germiline mutation results, we use the [`03_generate_somatic_script.py`](code/2_monopogen/03_generate_somatic_script.py) to generate SLURM scripts that will index VCF files from previous analyses and call somatic mutation using three different steps: `featureInfo`, `cellScan`, and `LDrefinement`.

Before running it, ensure you have:
- change the path in `export LD_LIBRARY_PATH=path/to/.anaconda3/envs/somaseqenv/lib:$LD_LIBRARY_PATH` to your conda lib. to fullfill the dependency requirements of Monopogen.
- This script requires the java to be installed or you have a module for java/17.0.6 (see the source code).

To generate the script:
```bash
python code/2_monopogen/03_generate_somatic_script.py --input_dir /path/to/input_data 
                                                    --output_dir /path/to/output_data
                                                    --slurm_scripts_dir /path/to/Monopogen
                                                    --genome_fa_path /path/to/STARSoloModifiedgenome.fa
                                                    --threads 64
```
Modify `/path/to/input_data`, `/path/to/output_data`, `/path/to/slurm_scripts` ,`/path/to/genome.fa` (the path to your modified genome FASTA file `Homo_sapiens.GRCh38.dna.primary_assembly.modified.fa`) and threads accordingly,

Then batch submit the SLURM scripts as [before](#batch-submission-of-sbatch-scripts).

### Preprocess for Analysis ###
#### 4. **vcf processing**: 

[`04_process_vcf_files.R`](code/3_preprocess_analysis/04_process_vcf_files.R) provide R functions to read the muliple results of Monopogen somatic mutation calling, combine with cell barcode information, use REDIPortal for RNA editing cite quality control, then incorporate Ensembl VEP (Variant Effect Predictor) results for genotyping annotation.

First, process groups of RDS files for downstream analysis using the `readAndCombineRDSWithDonorId` function.

```r
source("code/3_preprocess_analysis/04_process_vcf_files.R")
# Example of processing metadata and combining with Monopogen outputs
sample_info <- fread("path/to/metadata.csv")
barcode_info <- processBarcodeInfo("SOMA_seq/Monopogen/results", 
                                   "SOMA_seq/STARsolo/STAR_results/", 
                                   "/soloOut/GeneFull/filtered/barcode_counts.csv")

combinedDataList_whole_0.3 <- readAndCombineRDSWithDonorId("SOMA_seq/Monopogen/results",
sample_info)
processed_whole <- preprocessMonopogen(dataList = combinedDataList_whole, 
                                       info_meta = combinedMetadata, 
                                       barcode_info, 
                                       groupByFactors = "Subclass")

# Optional, save processed data for ploting
saveRDS(processed_whole, file = "path/to/processed_whole.rds")
```

For genomic studies, masking RNA editing sites is crucial as they can confound results of somatic mutation analyses. Integration of [REDIportal](http://srv00.recas.ba.infn.it/atlas/index.html) is recommended for this purpose. See instructions in [resources_to_download.txt](resources/resources_to_download.txt)

```r
processed_whole <- process_rna_editing("path/to/rna_editing_data.txt", processed_whole, sample_info)
```

Gene Information and VEP Integration: Combine genotyping information with gene annotations using data from Ensembl [VEP](https://useast.ensembl.org/Homo_sapiens/Tools/VEP?db=core;tl=VRi2Xtm6GBrxZ2Hb-9851365) (Variant Effect Predictor) to enhance the mutation analysis. The result including mutation type (e.g., missense) and gene-related impacts (e.g., intronic)

```r
# Save as input for VEP
prepare_input_for_vep(processed_whole, "/path/to/output/transformed.csv")
# Combine output
process_vep_results(processed_whole, "/path/to/ensembl_vep_output.txt", sample_info, "/path/to/output/processed_data.rds")
```

## Differential mutation analysis
### Somatic SNV burden ###
#### 5. **Mixed effect Model**

We aggregate mutation counts by cell type and donor to form pseudo-bulk data points. This pseudo-bulking is necessary to increase the statistical power of subsequent analyses. In our models, disease status, age, and other relevant covariates are treated as fixed effects, while variations attributable to individual donors are handled as random effects. This approach acknowledges the potential correlations among neurons derived from the same donor due to shared biological backgrounds.

**Model Formula**:
Note that the formula we used here is:
$$log10(total\; count) \sim (sub)Class + log10(read\;count) + log10(cell\;count) + log10(Age) + Sex + Disease + ROI+ (1|donor_{id})$$
Before running [`05_mixed_effect_model.R`](code/4_differential/05_mixed_effect_model.R), please ensure to customize the model's formula model <- lmer() to fit your dataset specifics and research questions.

Linear mixed-effects models were fitted, and P values from a t-test were calculated for each fixed effect as implemented in the lmerTest (v.3.1-3) R package (Kuznetsova et al., 2017). We also performed ANOVA to test the significance of individual fixed effects.
Pairwise comparison is performed using Tukey's HSD test to compare the significance of difference in all possible pairs of somatic mutation count means for cell types while controlling for the family-wise error rate.


```r
source("code/4_differential/04_process_vcf_files.R")
result <- perform_mixed_effect_model(meta_cell_info, processed_whole, sample_info, my_color)
model <- result$model # Extract the model from the results
# Output the summary, variance-covariance matrix, and ANOVA results of the model
summary(model)
vcov(model)   
anova(model)
# Perform pairwise comparisons using Tukey's Honest Significant Difference test
summary(glht(model, linfct = mcp(Subclass = "Tukey")))
```
**Statistical Tests Details**
- **Summary and Variance-Covariance Matrix**: Provides detailed output about the fixed effects and the random effects variance components.
- **ANOVA**: Tests the significance of the individual fixed effects included in the model.
- **Pairwise Comparisons**: Conducted using Tukey's HSD test, this is used to discern the significant differences in mutation counts between cell types, adjusted for multiple testing and controlling the family-wise error rate.

## Supplementary Scripts ##

- `generateBarcode.py`: Barcode Count Generator. This script produces Cellranger4-like `barcode_counts.csv` for STARSolo results. It also helps in assessing the quality of barcode tagging in sequencing experiments. Located in `utils/`, usage instructions are provided within the script.
- `colorControl.Rmd`: Plotting, Cell Type Mapping and Color Control. This helper R script is designed for snRNA-seq analysis to ensure consistent visualization of cell type classifications across publications. Located in `utils/`.
- `Figure`: This folder contain source code for plotting, note that this is for reproduction purpose, should be used along with the analysis code.