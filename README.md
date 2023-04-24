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
  - [1. Introduction to shinyBulk](#1-introduction-to-shinyBulk)
  - [2. Download and load data](#2-download-and-load-data)
  - [3. Prepare data for analysis](#3-prepare-data-for-analysis)
  - [4. Fit model to data](#4-fit-model-to-data)
  - [5. Make desired contrasts](#5-make-desired-contrasts)
  - [6. Identify differentially expressed genes](#6-identify-differentially-expressed-genes)
  - [7. Generate figures](#7-generate-figures)
 - [Walkthrough: Single-cell RNA seq analysis of two sample experiment](#walkthrough-single-cell-rna-seq-analysis-of-two-sample-experiment)
  - [1. Introduction to shinySC](#1-introduction-to-shinysc)
  - [2. Data Input](#2-data-input)
  - [3. Data Merging](#3-data-merging)
  - [4. Quality Control](#4-quality-control)
  - [5. Expression Exploration](#5-expression-exploration)
  - [6. Cluster Labels](#6-cluster-labels)
  - [7. Differential Expression](#7-differential-expression)
  - [8. Cell Cycle Analysis](#8-cell-cycle-analysis)
  - [9. Module Score](#9-module-score)
  - [10. Subset](#10-subset)
  - [11. Pseudotime](#11-pseudotime)
  - [12. Graph Generator(#12-graph-generator)
- [Contact](#contact)


##  1. Introduction to RATS

The RNA-seq analysis tool set (RATS) is an interactive shiny-based tool that leverages packages such as Seurat, Limma, EdgeR, EnrichR, Monocle3, and more to enable researchers with any level of bioinformatics expertise to conduct basic analysis of bulk and single-cell RNA sequencing data. 
This tool uses bioconductor packages [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to perform model fitting and differential gene expression analysis on bulk RNA-seq raw count data.

This tool also provides users the ability to generate several types of plots using [ggplot2](https://ggplot2.tidyverse.org/) and [pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap) packages. 

The sample data used to demonstrate this tool's use was generated from an experiment treating _Mycobacterium tuberculosis_ with DMSO or compound ZD at 40 and 80 micromolar concentrations.

##  2. Main functionalities of BulkRNAseqShiny

The basic pipeline of RATS single-cell analysis involves:
- Loading in sparse matrix data in form of Cell Ranger output, volocyto loom file, or Seurat R object.
- Merging samples for analysis (if not loading in Seurat object)
- performing quality control and filtering on merged data
- Identifying and labeling clusters
- Performing desired analyses, including differential expression, cell cycle, gene set score, subclustering, and pseudotime
- Generating figures to visualize analysis

The basic pipeline of RATS bulk analysis involves:
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
  devtools::install_github("YMreyoud/RATS")
```

##  2. Launch package
To launch the interactive package, use the following command:
``` r
  RATS::shinySC() # To launcht the single-cell analysis tool
  RATS::shinyBulk() #To launch the bulk analysis tool
```
This will open an browser-based window for interactive use of the tool.
Note: you may need to click 'open in browser' at the top of the window if the window is gray.

# Walkthrough: Bulkseq analysis of three condition experiment

##  1. Introduction to shinyBulk
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

## 1. Introduction to shnySC
This singe-cell RNA seq analysis tool creats a user friendly interface to allow exploration of single-cell rna sequencing data. In this walkthrough, I will demonstrate the use of this tool to analyze scRNAseq data from an experiment consisting of counts from two samples of GM-CSF cultured murine bone marrow cells derived from two genotypes (WT and Bhlhe40-/-). This tool functions using various packages for analysis including Seurat, EnrichR, and Monocle3.

## 2. Data Input
scShiny() accepts three different types of data as input. 
- Option 1: Load 10x
  - For this input, you must select a parent folder containing a subfolder for each sample you would like to analyze. In these subfolders, there must be three files: barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz. These files are generated from CellRanger.
- Option 2: Load Loom
  - This option allows you to select loom files generated using velocyto. This filetype is useful for retaining spliced and unspliced read data for downstream analysis such as scVelo velocity analysis. 
- Option 2: Load Model
  - This option allows the user to upload an R object of a previously generated seurat object. This is useful for loading in any previously analyzed data for continued exploration. 

## 3. Data Merging
scShiny() offers two methods for merging data samples into one seurat object. These usefulness of these two methods depends on your specific experiment.
- Simple Merge:
  - This merge methods, as the title suggests, simple merges the count matrices from the experiment without accounting for any variation due to batch effects. We recommend only using this type of merge if all samples were run together for sequencing.
- Integrate:
  - This merge method accounts for variation due to batch effects by finding integration anchors between samples and using this anchors to normalize the data to eliminate batch effects. This method should always be used for data that was sequenced in different batches. 

After successfully merging, scShiny() will output a series of graphs exploring the quality of the data. These graphs include information regarding mitrochondrial read percentace, feature counts, RNA counts, PCA elbow plot, UMAP plot, and more. These graphs can then be used in the following section to set appropriate quality control parameters. 

## 4. Quality Control
Quality control is an important part of analyzing any sequencing data. The primary goal of this QC step is to remove low quality reads, multiplets, and dying cells. This is done through a seq of filters that act act on count number, feature number, and percent mitochondrial reads.
The filters are as follows:
- Cutoffs for number of counts:
  - Cells with counts that are too low or too high indicate low quality cells or multiplets. Set your cutoffs according to the histogram of RNA counts in the 'Merging' section.
- Cutoffs for number of features:
  - Once again, low feature counts indicate poor quality reads, and high feature counts indicate multiplets. Use the graphs in the 'Merging' section to appropriate set your fulters for your specific data.
- Max mitochondrial read percentage:
  - Cells with a high percentage of mitochondrial reads are assumed to be dying and should be removed. Use the histogram in the 'Merging' section to set this cutoff. 
- Resolution for clustering:
  - This is th eresolution from 0 to 1 to use for determining clusters within the data. A higher resolution results in a higher number of more specific clusters; this is useful for fidning very specific cell types with little variation form one another. A low resolution is less sensitive, but more useful for finding general cell-types with grater variation. 

## 5. Expression Exploration
This section is useful for exploring the gene expression within your clusters to help classify the different clusters based on their expression. The Featureplot subsection allows you to view the expression of genes on a UMAP representation. The grah for this subsection is shown under the 'Interactive Feature Plot' tab. The second functionality of this section is annotation using SingleR. 

There are two ways to annotate your cells:
- Celldex
  - This option allows you to select a reference from the Celldex database. This is the recommended method for cell annotation. Simply select the desired dataset from the dropdown and click 'Run SingleR Labeling'
- Custom reference
  - Alternatively, users can upload their own .rds reference file in SingleCellExperiment format to use for annotation.
Once the annotation is done, two tables will be populated in the output. These tables will show the number of cells belonging to each cell type in each cluster. This table can then be used in corroboration with the feature plot to determine cluster identities. Clusters can then be re-labeled in the next section.

## 6. Cluster Labels
Once the identifiy of your clusters has been determined from the previous section, the user can use this section to relabel the clusters based on their annotations. Simply type in the label for each cluster, then click 'Re-Label'. A series of UMAP plots will appear showing the SingleR annotation, the Clusters, and the user-defined labels. 

## 7. Differential Expression
This section allows the user to perform differential expression analysis and enrichment analysis using EnrichR. To perform DGE, the user must first select the idents hat contain their desired comparison. In most cases this will euther be 'seurat_clusters' or 'celltype.condition'. Next, the user my select the two identities whose expression the user would like to compare. the format for this is x vs y. Once these selections are made, the table of differentially expressed genes will be populated. The yser can download this table for further analysis. 

Next, the user can perform enrichment analysis on this differential expression table using the EnrichR subtab. Simply select the list of geneset databases you would like to use, and clicl 'Run EnrichR'. This will generate postive enrichment and negative enrichment graphs for each database selected. These can be viewed in the 'Positive Enrichment' and 'Negative Enrichment' output subtabs.

## 8. Cell Cycle Analysis
This section simply allows the user to perfom a cell-cycle analysis, which scores each cell based on its expression of cell cycle markers, and classifies it as being in S phase, G1, or G2. Simply click 'Run Analysis' to perform the analysis. Results can be viewed in dimplot format on the right side, and can be split by condition in the 'By Condition' output subtab.

## 9. Module Score
Similar to the cell cycle analysis, this section allows the user to upload two gene sets defining opposing cell states, define a module name, then score  each of the cells in the gene set on their expression levels of those two gene sets. This can then be visualized as a blended feature plot using the Visualize Module tab.

To perform this analysis follow these steps:
- Upload first geneset (M1.csv, found in the sample data folder) which contains genes that define pro-inflammatory macrophage polarization
- Upload second geneset (M2.csv, found in the sample data folder) which contains genes that define anti-inflammatory macrophage polarization
- Add desired module name. 'Polarization' for this example.
- Click 'Add Module'
- To visualize the results, click the 'Visualize Module' subtab and select the newly added idents, which, in this case, will be 'Polarization1', and 'Polarization2'.
- Click the 'Module Graph' output subtab to view the blended featureplot.

## 10. Subset
This section allows the user to subset their data for further analysis. For example, if a user would like to further subcluster and analyze their macrophage population, they would select 'cell_type' as the Ident, and then select 'Macrophages' as their values. The user also has the option to recluster by setting 'Recluster' to TRUE, and then setting the appropriate clustering resolution. This subset can then be used in any of the other analysis sections. To revert back to the original object, simply click 'Restore Parent'. (Note: the 'Restore Parent' feature only retains the most recent object in memory, thus, if you subset an object twice in a row without restoring the parent in between, the parent data will be lost. I recommend downloading your object from the quality control tab before subsetting.)

## 11. Pseudotime
This section allows the user to leverage the monocle3 package to perform pseudotime analysis. Simply select the gene you desire to track, and click 'Run Pseudotime'. This wil generate a pseudotime graph which can be downloaded. 

## 12. Graph Generator
This section gives the user the ability to generate many different types of graphs and figures. The 'Graph Options' subtab is used to define the parameters of the graph, while the 'Download Options' subtab is used to set the download options including size and file type. 
