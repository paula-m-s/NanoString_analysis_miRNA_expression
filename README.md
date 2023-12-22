# Nanostring Data Analysis - miRNA
## Overview
This repository contains R scripts for the analysis of miRNA Nanostring data, focusing on differential expression analysis. The analysis was conducted by Paula Morales-Sanchez, and the scripts are designed to provide insights into the dysregulation of miRNAs in various conditions, specifically tumors with DICER1 and DGCR8 mutations.
[![Canonical miRNA biogenesis pathway](https://github.com/paula-m-s/NanoString_analysis_d1_d8_thyroid/blob/main/Canonical_pathway_miRNA.png "Canonical miRNA biogenesis pathway")](http:/https://github.com/paula-m-s/NanoString_analysis_d1_d8_thyroid/blob/main/Canonical_pathway_miRNA.png/ "Canonical miRNA biogenesis pathway")
## Pipeline Information
### Authors
Paula Morales-Sanchez
### Pipeline Support
1. [NanoTube Bioconductor Package](http://https://www.bioconductor.org/packages/release/bioc/vignettes/NanoTube/inst/doc/NanoTube.html "NanoTube Bioconductor Package")
2. [CBCS Normalization GitHub Repository](http://https://github.com/bhattacharya-a-bt/CBCS_normalization/ "CBCS Normalization GitHub Repository")
3. [Published Article: An approach for normalization and quality control for NanoString RNA expression data.  Bhattacharya et al., 2021](http://https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8138885/ "Published Article: An approach for normalization and quality control for NanoString RNA expression data.  Bhattacharya et al., 2021")
### Date Information

Date Started: April 2023

Last Date Modified: December 2023

### R Version
The scripts are compatible with the R programming language. R version 4.1.3 (2022-03-10)

### Scripts
The scripts performs a comprehensive analysis of miRNA Nanostring data, focusing on differential expression analysis. It explores the dysregulation of miRNAs in tumors with DICER1 and DGCR8 mutations. Here is a brief explanation for each of the files:

1. 00_pipeline_helpers.R: This file contains functions and helper utilities used across other scripts in the pipeline. Typically, common functions, global configurations, or anything that can be shared among multiple scripts are placed here.
2. 01_Raw_data_processing_and_QC.R: This file is dedicated to raw data processing and quality control (QC). Here, you might find code related to loading raw data, initial cleaning, data exploration, and any steps necessary to prepare raw data for analysis.
3. 02_Normalization_and_BatchCorrection.R: In this script, normalization and batch correction of the data are performed. Normalization is an important step to compare samples and ensure that results are not biased by technical differences. Batch correction addresses potential batch effects that could influence results.
4. 03_DE_analysis.R: This file focuses on the differential expression (DE) analysis. Here, comparisons between groups are made, statistical tests are applied, and genes or features showing significant changes in expression under different conditions are identified. Also different types of graphs displayed across the manuscript are referenced here.

### Disclaimer
The analysis is based on an ongoing project leaded by Dr. Barbara Rivera. Relevant sources and references are provided in the script and the manuscript.

### Note
Feel free to reach out for any questions or clarifications. Contributions and feedback are welcome!

