BulkRNAseqShiny
================
Yassin Mreyoud
03-09-2023

- [1. Introduction to BulkRNAseqShiny](#1-introduction-to-bulkrnaseqshiny)
- [2. Main functionalities of BulkRNAseqShiny](#2-main-functionalities-of-bulkrnaseqshiny)
- [Installation](#installation)
  - [1. Install the package from github](#1-install-the-package-from-github)
  - [2. Launch package](#2-launch-package)
- [Walkthrough: Bulkseq analysis of three condition experiment](#walkthrough-bulkseq-analysis-of-three-condition-design)
  - [1. Introduction](#1-introduction)
  - [2. Download and load data](#2-download-and-load-data)
  - [3. Prepare data for analysis](#3-prepare-data-for-analysis)
  - [4. Fit model to data](#4-fit-model-to-data)
  - [5. Make desired contrasts](#5-make-desired-contrasts)
  - [6. Identify differentially expressed genes](#6-identify-differentially-expressed-genes)
  - [7. Generate figures](#7-generate-figures)
- [Contact](#contact)
##  1.  Introduction to BulkRNAseqShiny

This tool uses bioconductor packages [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to perform model fitting and differential gene expression analysis on bulk RNA-seq raw count data.

This tool also provides users the ability to generate several types of plots using [ggplot2](https://ggplot2.tidyverse.org/) and [pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap) packages. 

The sample data used to demonstrate this tool's use was generated from an experiment treating _Mycobacterium tuberculosis_ with DMSO or compound ZD at 40 and 80 micromolar concentrations.

##  2.  Main functionalities of BulkRNAseqShiny

The basic pipeline of BulkRNAseqShiny analysis involves:
- Loading raw count data from bulk RNA sequencing experiment.
- pre-processing data to define columns containing samples, genes, and conditions.
- Fitting a model to the data using limma or voom.
- Declaring which contrasts to make.
- Performing differential gene expression analysis.
- Generating figures.

# Installation

## 1. Install the package from github

This package has been tested using R version 4.2.2. This tool also requires instllation of devtools and BiocManager if not already installed.
``` r
  # Install devtools and BiocManager if not already installed
  require(devtools,BiocManager)
  
  # Install BulkRNAseqShiny
  devtools::install_github("YMreyoud/BulkRNAseqShiny")
```

##  2.  Launch package
To launch the interactive package, use the following command:
``` r
  BulkRNAseqShiny::shinyBulk()
```
This will open an browser-based window for interactive use of the tool.
Note: you may need to click 'open in browser' at the top of the window if the window is gray.

# Walkthrough: Bulkseq analysis of three condition experiment

##  1. Introduction
In this walkthrough we will perform a standard bulk RNA-sequencing analysis on data generated from an experiment treating _Mycobacterium tuberculosis_ with one of three conditions:
- ZD 40um
- ZD 80um
- DMSO
The samples have been sequenced using the Illumina Noviseq system, aligned using STAR, and quality controlled using FastQC.
This data has not yet been published but can be downloaded by ollowing this [link](https://github.com/YMreyoud/BulkRNAseqShiny/data/all.gene_counts.tsv)

## 2. Download and load data

Once the all.gene_counts.tsv file has been downloaded from this [link](https://github.com/YMreyoud/BulkRNAseqShiny/data/all.gene_counts.tsv) the data can be imported into the tool by using the file upload field. Next, select which type of file is being uploaded (.csv or .tsv) and upload the data by clicking the 'Upload' button.

After the data is uploaded, it can be browsed in the right window.

![](https://github.com/YMreyoud/BulkRNAseqShiny/images/input.gif)

## 3. Prepare data for analysis
The user must now define the column that contains the gene information. This can be gene name, gene ID, or any other gene-specific information you would like to use. For this walkthrough we will use the ensemble_gene_id column.
- Use the drop down menu to select the column containing your gene-specific information
- List all conditions you would like to analyze separated by commas. For this experiment our conditions are: ZD40, ZD80, DMSO
- Select the columns containing count data for the samples you would like to analyze
- Use the drop-down menus to label each sample with it's corresponding condition. 

![](https://github.com/YMreyoud/BulkRNAseqShiny/images/prep.gif)


