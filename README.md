RNA-seq Analysis Tool Set (RATS)
================
Yassin Mreyoud
03-09-2023

- [1. Introduction to BulkRNAseqShiny](#1-introduction-to-bulkrnaseqshiny)
- [2. Main functionalities of BulkRNAseqShiny](#2-main-functionalities-of-bulkrnaseqshiny)
- [Installation](#installation)
  - [1. Install the package from github](#1-install-the-package-from-github)
  - [2. Launch package](#2-launch-package)
- [Walkthrough: Bulkseq analysis of three condition experiment](#walkthrough-bulkseq-analysis-of-three-condition-experiment)
  - [1. Introduction](#1-introduction)
  - [2. Download and load data](#2-download-and-load-data)
  - [3. Prepare data for analysis](#3-prepare-data-for-analysis)
  - [4. Fit model to data](#4-fit-model-to-data)
  - [5. Make desired contrasts](#5-make-desired-contrasts)
  - [6. Identify differentially expressed genes](#6-identify-differentially-expressed-genes)
  - [7. Generate figures](#7-generate-figures)
 - [Walkthrough: Single-cell RNA seq analysis of two sample experiment](#walkthrough-single-cell-rna-seq-analysis-of-two-sample-experiment)
  - [1. Introduction](#1-introduction)
- [Contact](#contact)


##  1. Introduction to BulkRNAseqShiny

This tool uses bioconductor packages [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to perform model fitting and differential gene expression analysis on bulk RNA-seq raw count data.

This tool also provides users the ability to generate several types of plots using [ggplot2](https://ggplot2.tidyverse.org/) and [pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap) packages. 

The sample data used to demonstrate this tool's use was generated from an experiment treating _Mycobacterium tuberculosis_ with DMSO or compound ZD at 40 and 80 micromolar concentrations.

##  2. Main functionalities of BulkRNAseqShiny

The basic pipeline of BulkRNAseqShiny analysis involves:
- Loading raw count data from bulk RNA sequencing experiment.
- pre-processing data to define columns containing samples, genes, and conditions.
- Fitting a model to the data using limma or voom.
- Declaring which contrasts to make.
- Performing differential gene expression analysis.
- Generating figures.

# Installation

## 1. Install the package from github

This package has been tested using R version 4.2.2. This tool also requires instllation of R.utils, devtools, and BiocManager if not already installed.
``` r
  # Install devtools and BiocManager if not already installed
  install.packages(c('devtools','BiocManager', 'R.utils')
  devtools::install_version("Matrix", version = "1.5.3", force = TRUE) #Make sure Matrix version is 1.5.3 or newer!
  
  # Install BulkRNAseqShiny
  devtools::install_github("YMreyoud/BulkRNAseqShiny")
```

##  2. Launch package
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

![](https://github.com/YMreyoud/BulkRNAseqShiny/blob/main/images/input.gif)

## 3. Prepare data for analysis
The user must now define the column that contains the gene information. This can be gene name, gene ID, or any other gene-specific information you would like to use. For this walkthrough we will use the ensemble_gene_id column.
- Use the drop down menu to select the column containing your gene-specific information
- List all conditions you would like to analyze separated by commas. For this experiment our conditions are: ZD40, ZD80, DMSO
- Select the columns containing count data for the samples you would like to analyze
- Use the drop-down menus to label each sample with it's corresponding condition. 



![](https://github.com/YMreyoud/BulkRNAseqShiny/blob/main/images/prep.gif)

## 4. Fit model to data
Simply use the drop down to select which method to use for data fitting. Voom is recommended for most application, limma should only be used in cases where library size is consistent.

![](https://github.com/YMreyoud/BulkRNAseqShiny/blob/main/images/model.gif)

## 5. Make desired contrasts
Use the selection box to pick contrasts one at a time. For example, to compare ZD40 vs DMSO, select ZD40 then DMSO, and click 'Add Contrast' repeat with all the comparisons you would like to make then click 'Done'. If a mistake is made, use the 'Clear' button to reset the list of contrasts. Current contrasts will be shown in the right window separated by commas.

![](https://github.com/YMreyoud/BulkRNAseqShiny/blob/main/images/contrasts.gif)

## 6. Identify differentially expressed genes
Select the desired comparison from the drop down menu, and set the remaining parameters as desired. The recommended adjust method is BH (Benjamini-Hochberg). Click 'Done' to generate the list of differentially expressed genes.

A table containing all differentially expressed genes matching the user-defined p-value cutoff for the selected contrast will appear in the right window. To download a csf file of the table, click the 'Download table' button.

![](https://github.com/YMreyoud/BulkRNAseqShiny/blob/main/images/dge.gif)

## 7. Generate figures
This tool allows users to generate volcano plots and heatmaps. Use the 'Graph Type' drop down to select which type of graph you would like to generate. Generated graphs will appear on the right side, and can be downloaded using the 'Download' button. Additionally, download options are available for the volcano plots in the 'Download Options' tab. 

### Volcano plot:
- P value cut-off: select the p-value cutoff to differentiate between significant and non-significant genes
- Genes (optional): select genes that would like to be highlighted with a blue color
- Geneset (optional): upload a list of genes to be highlighted with a blue color
- Contrast: select the comparison you would like to plot
- Label: If TRUE, any genes provided in the previous fields will be labeled.



![](https://github.com/YMreyoud/BulkRNAseqShiny/blob/main/images/volcano.gif)

### Heatmap
- Genes: select the genes you would like to be displayed on the heatmap
- Conditions: select the conditions you would like to be displayed on the heatmap



![](https://github.com/YMreyoud/BulkRNAseqShiny/blob/main/images/heatmap.gif)


# Walkthrough: Single-cell RNA seq analysis of two sample experiment

# 1. Introduction:
This singe-cell RNA seq analysis tool creats a user friendly interface to allow exploration of single-cell rna sequencing data. In this walkthrough, I will demonstrate the use of this tool to analyze scRNAseq data from an experiment consisting of counts from two samples of GM-CSF cultured murine bone marrow cells derived from two genotypes (WT and Bhlhe40-/-). This tool functions using various packages for analysis including Seurat, EnrichR, and Monocle3.

# 2. Data Input:
scShiny() accepts three different types of data as input. 
- Option 1: Load 10x
  - For this input, you must select a parent folder containing a subfolder for each sample you would like to analyze. In these subfolders, there must be three files: barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz. These files are generated from CellRanger.
- Option 2: Load Loom
  - This option allows you to select loom files generated using velocyto. This filetype is useful for retaining spliced and unspliced read data for downstream analysis such as scVelo velocity analysis. 
- Option 2: Load Model
  - This option allows the user to upload an R object of a previously generated seurat object. This is useful for loading in any previously analyzed data for continued exploration. 

# 3. Data Merging:
scShiny() offers two methods for merging data samples into one seurat object. These usefulness of these two methods depends on your specific experiment.
- Simple Merge:
  - This merge methods, as the title suggests, simple merges the count matrices from the experiment without accounting for any variation due to batch effects. We recommend only using this type of merge if all samples were run together for sequencing.
- Integrate:
  - This merge method accounts for variation due to batch effects by finding integration anchors between samples and using this anchors to normalize the data to eliminate batch effects. This method should always be used for data that was sequenced in different batches. 
