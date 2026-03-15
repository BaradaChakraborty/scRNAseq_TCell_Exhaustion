# Computational Dissection of T-cell Exhaustion in Solid Tumors: A Single-Cell RNA-seq Analysis

**Author:** Barada Chakraborty  
**Affiliation:** Bioinformatics Research Project (Guided), BioGrademy (BioNovura Labs)  
**Certification:** Coding for Biologists Program (Certificate No. BL25CFB15048)  

## Project Overview
This repository contains a comprehensive computational pipeline engineered to analyze the tumor microenvironment (TME) of solid tumors (Melanoma). The primary objective is to map the continuous transition of CD8+ T-cells into terminal exhaustion and validate these transcriptomic signatures against clinical therapy outcomes. 

**Repository Versioning Note:** *An earlier iteration of this sandbox included extended independent modules (such as Random Forest machine learning feature selection and whole-transcriptome GSEA). To ensure strict adherence, transparency, and verifiable alignment with the specific parameters of the BioGrademy BL25CFB15048 guided syllabus, the repository history and file structure were formally reset to focus exclusively on the core scRNA-seq and ICB validation pipeline described below.*

## Core Pipeline Architecture

The framework is built across three distinct modules, demonstrating cross-language proficiency:

### 1. Data Retrieval (Linux/Bash)
* **Script:** `01_Linux_Data_Retrieval.sh`
* **Function:** Automated command-line retrieval of raw single-cell transcriptomic matrices (GSE115978) from server databases using `curl`, establishing secure local data directories.

### 2. Single-Cell Analytics (R)
* **Script:** `02_scRNAseq_QC_and_UMAP.R`
* **Function:** Core processing utilizing the `Seurat` ecosystem. 
    * Executed strict quality control (QC) filtering and LogNormalization.
    * Performed dimensionality reduction (PCA) and high-dimensional clustering.
    * Mapped canonical exhaustion signatures (e.g., *PDCD1*, *HAVCR2*, *LAG3*) onto the cellular manifold utilizing both UMAP and t-SNE algorithms.

### 3. Clinical Biomarker Validation (Python)
* **Script:** `03_Python_ICB_Validation.py`
* **Function:** Integrated molecular transcriptomic findings with clinical metadata.
    * Utilized `pandas`, `numpy`, and `scipy` to evaluate T-cell exhaustion as a predictive biomarker for Immune Checkpoint Blockade (Anti-PD1) response.
    * Generated statistical validation visualizations demonstrating that high baseline exhaustion heavily correlates with therapy non-responders (p < 0.001).

## System Requirements
* **Linux/Unix CLI** (Git Bash for Windows)
* **R (v4.0+)**: `Seurat`, `dplyr`, `ggplot2`
* **Python (3.9+)**: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`