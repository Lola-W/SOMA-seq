# SOMA-seq code repository
This repo is a place holder for running the SOMA-seq pipeline and reproducing the presented analyses.

## Table of Contents

[//]: # 

* [Introduction](#introduction)
* [Installation](#installation)
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
`R version 4.3.1 (2023-06-16)`, `Platform: x86_64-pc-linux-gnu`and `aarch64-apple-darwin22.4.0, and `Running under: `x86_64, linux-gnu`, `aarch64`, and `darwin22.4.0`.

The pipeline consists of three main stages:

1. Somatic mutation calling: Utilizing Monopogen, a state-of-the-art SNV calling package developed by [Ken chen's lab](https://www.mdanderson.org/research/departments-labs-institutes/labs/ken-chen-laboratory.html), Soma-seq accurately identifies somatic mutations from snRNA-seq data. Monopogen is specifically designed to handle the unique challenges of single-cell sequencing datasets, such as sparsity and allelic dropout, ensuring high-confidence mutation calls.

2. Differential mutation analysis: Soma-seq performs a comprehensive analysis of the consequences of differential mutations across various cell types and disease groups. By employing advanced statistical methods and mixed-effect modeling, the pipeline accounts for potential confounding factors, such as age, read depth, and donor-specific effects, to uncover cell type-specific mutation patterns and rates.

3. Pathogenicity analysis: Leveraging the power of AlphaMissense predictions, Soma-seq assesses the pathogenic potential of identified somatic mutations. This stage provides valuable insights into the functional impact of mutations, aiding in the identification of potentially deleterious variants that may contribute to disease pathogenesis.

The output of Soma-seq includes raw csv summary of putative mutations, summaries of single base substitution patterns, mutation abundance across cell types, ranked enriched gene lists, pathway enrichment results, and a list of deleterious somatic mutations annotated by AlphaMissense. These comprehensive results enable researchers to gain a deeper understanding of the role of somatic mutations in various diseases and cell types, facilitating the discovery of potential disease biomarkers and the inference of underlying causal factors.

## Installation
SOMA-seq is avaiable on github
```
git clone https://github.com/Lola-W/SOMA-seq.git
cd SOMA-seq
```
Setup the conda environment after cloning this repo:

**Note:** If you encounter any issues with the environment path in the `environment.yml` file, please modify the `prefix` field to match your local conda installation path. This can be done by replacing the `prefix` path in `environment.yml` with your desired conda environments directory path. For example:
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

