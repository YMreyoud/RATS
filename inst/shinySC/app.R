#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# library(shiny)
# library(Seurat)
# library(SeuratWrappers)
# library(ggpubr)
# library(celldex)
# library(SingleR)
# library(patchwork)
# library(tidyverse)
# library(cowplot)
# library(stringr)
# library(enrichR)
# library(SeuratData)
# library(SeuratObject)
# library(scales)
# library(ggrepel)
# library(msigdbr)
# library(shinythemes)
library(enrichR)
library(Matrix)


read.loom.matrices <- function(file, engine='hdf5r') {
  if (engine == 'h5'){
    cat('reading loom file via h5...\n')
    f <- h5::h5file(file,mode='r');
    cells <- f["col_attrs/CellID"][];
    genes <- f["row_attrs/Gene"][];
    dl <- c(spliced="/layers/spliced",unspliced="/layers/unspliced",ambiguous="/layers/ambiguous");
    if("/layers/spanning" %in% h5::list.datasets(f)) {
      dl <- c(dl,c(spanning="/layers/spanning"))
    }
    dlist <- lapply(dl,function(path) {
      m <- as(f[path][],'dgCMatrix'); rownames(m) <- genes; colnames(m) <- cells; return(m)
    })
    h5::h5close(f)
    return(dlist)
  } else if (engine == 'hdf5r') {
    cat('reading loom file via hdf5r...\n')
    f <- hdf5r::H5File$new(file, mode='r')
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Gene"]][]
    dl <- c(spliced="layers/spliced",
            unspliced="layers/unspliced",
            ambiguous="layers/ambiguous")
    if("layers/spanning" %in% hdf5r::list.datasets(f)) {
      dl <- c(dl, c(spanning="layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[[path]][,]),'dgCMatrix')
      rownames(m) <- genes; colnames(m) <- cells;
      return(m)
    })
    f$close_all()
    return(dlist)
  }
  else {
    warning('Unknown engine. Use hdf5r or h5 to import loom file.')
    return(list())
  }
}



ReadVelocity <- function(file, engine = 'hdf5r', verbose = TRUE) {
  #CheckPackage(package = 'velocyto-team/velocyto.R', repository = 'github')
  if (verbose) {
    sink(file = stderr(), type = 'output')
    on.exit(expr = sink())
    ldat <- read.loom.matrices(file = file, engine = engine)
  } else {
    invisible(x = capture.output(ldat <- vread.loom.matrices(
      file = file,
      engine = engine
    )))
  }
  return(ldat)
}


`%>%` <- magrittr::`%>%`
mouse_mito <- readr::read_csv(file = "Mouse.MitoCarta3.0.csv") %>% # Import list of mouse mitochondrial genes
  dplyr::select(c(Symbol, MitoCarta2.0_Score)) %>%
  dplyr::arrange(desc(MitoCarta2.0_Score)) %>%
  head(200)
human_mito <- readr::read_csv(file = "Human.MitoCarta3.0.csv") %>% # Import list of mouse mitochondrial genes
  dplyr::select(c(Symbol, MitoCarta2.0_Score)) %>%
  dplyr::arrange(desc(MitoCarta2.0_Score)) %>%
  head(200)
mito.genes <- c(as.character(mouse_mito$Symbol), as.character(human_mito$Symbol))


renew_merged <- function(names, paths, read10x) {
  progress <- shiny::Progress$new(max = (2*length(paths)+7), min = 0)
  on.exit(progress$close())
  progress$set(message= "Merging data...", value = 0)
  list <- list()
  cell.ids <<- list()
  if (!read10x) {
    for (x in 1:length(paths)) {
      progress$inc(1, detail = paste("Reading file: ", names[x]))
      temp.data <- ReadVelocity(file = paths[x])
      temp <- Seurat::as.Seurat(temp.data)
      temp[["RNA"]] <- temp[["spliced"]]
      temp.mito <- mito.genes[mito.genes %in% rownames(temp)]
      temp[["percent.mt"]] <- colSums(temp[temp.mito, ]) / colSums(temp) * 100
      #temp[["percent.mt"]] <- Seurat::PercentageFeatureSet(temp, pattern = "^MT-")
      temp@meta.data[, "condition"] <- stringr::str_replace(names[x], ".loom", "")
      condition <- stringr::str_replace(names[x], ".loom", "")
      # assign(condition, temp)
      list <- c(list, temp)
      cell.ids <<- c(cell.ids, condition)
    }
  }
  else {
    for (x in 1:length(paths)) {
      progress$inc(1, detail = paste("Reading file: ", names[x]))
      temp.data <- Seurat::Read10X(data.dir = paths[x])
      temp <- Seurat::CreateSeuratObject(temp.data)
      temp.mito <- mito.genes[mito.genes %in% rownames(temp)]
      #temp[["percent.mt"]] <- colSums(temp[temp.mito, ]) / colSums(temp) * 100
      temp[["percent.mt"]] <- Seurat::PercentageFeatureSet(temp, pattern = "^MT-")
      temp@meta.data[, "condition"] <- names[x]
      condition <- names[x]
      # assign(condition, temp)
      list <- c(list, temp)
      cell.ids <<- c(cell.ids, condition)
    }
  }
  for (i in seq_along(list)) {
    progress$inc(1, detail = paste("Normalizing data for datset ", i, ' of ', length(list)))
    list[[i]] <- Seurat::NormalizeData(list[[i]], verbose = FALSE)
    list[[i]] <- Seurat::FindVariableFeatures(
      list[[i]],
      selection.method = "vst",
      nfeatures = 5000,
      verbose = FALSE
    )
  } # vst = variance stabilizing transformation

  # merged <- merge(list[[1]], list[-1],
  # add.cell.ids = cell.ids,
  # project = 'merged',
  # merge.data = TRUE)
  progress$inc(1, detail = 'Finding anchors...')
  anchors <- Seurat::FindIntegrationAnchors(object.list = list, dims = 1:30)
  progress$inc(1, detail = 'Merging data...')
  integrated <- Seurat::IntegrateData(anchorset = anchors, dims = 1:30)
  merged <- integrated
  # Run the standard workflow for visualization and clustering
  Seurat::DefaultAssay(merged) <- "integrated"
  progress$inc(1, detail = 'Scaling data...')
  merged <- Seurat::ScaleData(merged)
  progress$inc(1, detail = 'Running PCA...')
  merged <- Seurat::RunPCA(merged, npcs = 30, verbose = FALSE)
  progress$inc(1, detail = 'Running UMAP...')
  merged <- Seurat::RunUMAP(merged, reduction = "pca", dims = 1:30)
  progress$inc(1, detail = 'Finding Clusters...')
  merged <- Seurat::FindNeighbors(merged, reduction = "pca", dims = 1:30)
  merged <- Seurat::FindClusters(merged, resolution = 0.5)
  Seurat::DefaultAssay(merged) <- "RNA"
  progress$inc(1, detail = 'Finishing up...')
  merged <- Seurat::SCTransform(merged, verbose = TRUE)
  merged <- Seurat::RunPCA(merged, npcs = 30, verbose = FALSE)
  merged <- Seurat::RunUMAP(merged, reduction = "pca", dims = 1:30)
  # saveRDS(merged, "files/merged.rds")

  return(merged)
}

renew_integrated <- function(names, paths, read10x) {
  progress <- shiny::Progress$new(max = (2*length(paths)+7), min = 0)
  on.exit(progress$close())
  progress$set(message= "Integrating data...", value = 0)
  list <- list()
  cell.ids <<- list()
  if (!read10x) {
    for (x in 1:length(paths)) {
      progress$inc(1, detail = paste("Reading file: ", names[x]))
      temp.data <- ReadVelocity(file = paths[x])
      temp <- Seurat::as.Seurat(temp.data)
      temp[["RNA"]] <- temp[["spliced"]]
      temp[["percent.mt"]] <- Seurat::PercentageFeatureSet(temp, pattern = "^MT-")
      temp@meta.data[, "condition"] <- stringr::str_replace(names[x], ".loom", "")
      condition <- stringr::str_replace(names[x], ".loom", "")
      # assign(condition, temp)
      list <- c(list, temp)
      cell.ids <<- c(cell.ids, condition)
    }
  }
  else {
    for (x in 1:length(paths)){
      progress$inc(1, detail = paste("Reading file: ", names[x]))
      temp.data <- Seurat::Read10X(data.dir = paths[x])
      temp <- Seurat::CreateSeuratObject(temp.data)
      temp[["percent.mt"]] <- Seurat::PercentageFeatureSet(temp, pattern = "^MT-")
      temp@meta.data[,"condition"] <- names[x]
      condition <- names[x]
      list <- c(list, temp)
      cell.ids <<- c(cell.ids, condition)
    }
  }
  for (i in seq_along(list)) {
    progress$inc(1, detail = paste("Normalizing data for datset ", i, ' of ', length(list)))
    list[[i]] <- Seurat::SCTransform(list[[i]], verbose = TRUE)
  } # vst = variance stabilizing transformation
  progress$inc(1, detail = 'Selecting integration features...')
  features <- Seurat::SelectIntegrationFeatures(list, nfeatures = 5000)
  progress$inc(1, detail = 'Preparing integration...')
  list <- Seurat::PrepSCTIntegration(list, anchor.features = features)
  # integrated <- merge(list[[1]], list[-1],
  # add.cell.ids = cell.ids,
  # project = 'integrated',
  # merge.data = TRUE)
  progress$inc(1, detail = 'Finding integration anchors...')
  anchors <- Seurat::FindIntegrationAnchors(object.list = list, dims = 1:30, normalization.method = "SCT", anchor.features = features)
  progress$inc(1, detail = 'Integrating data...')
  integrated <- Seurat::IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "SCT")
  # Run the standard workflow for visualization and clustering
  Seurat::DefaultAssay(integrated) <- "integrated"
  progress$inc(1, detail = 'Running PCA...')
  integrated <- Seurat::RunPCA(integrated, npcs = 30, verbose = TRUE)
  progress$inc(1, detail = 'Running UMAP...')
  integrated <- Seurat::RunUMAP(integrated, reduction = "pca", dims = 1:30)
  progress$inc(1, detail = 'Clustering...')
  integrated <- Seurat::FindNeighbors(integrated, reduction = "pca", dims = 1:30)
  integrated <- Seurat::FindClusters(integrated, resolution = 0.5)
  # saveRDS(integrated, "files/integrated.rds")

  return(integrated)
}

qc <- function(seuratobject, featuremin = 500, featuremax = 6500, countmin = 1000, countmax = 7.5e4, percentmt = 5, res = 0.5) {

  progress <- shiny::Progress$new(max = 6, min = 0)
  on.exit(progress$close())
  progress$set(message= "Performing QC...", value = 0)

  # This removes cells that don't meet our filtering criteria. Upper bound filters out doublets, lower bound filters out empty bubbles
  print(c(featuremin, featuremax, countmin, countmax, percentmt))

  progress$inc(1, detail = 'Filtering cells...')

  seuratobject <- seuratobject %>%
    subset(nFeature_RNA > featuremin & nCount_RNA > countmin &
      percent.mt < percentmt & nFeature_RNA < featuremax & nCount_RNA < countmax)


  # Add cell cycle scores, can be useful for analysis.
  # s.genes <- cc.genes$s.genes
  # g2m.genes <- cc.genes$g2m.genes
  # seuratobject <- CellCycleScoring(seuratobject, search=TRUE, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  # Next three lines identify apoptotic cells.
  # death_gene_sets<-msigdbr(category = c("H"))
  # death_gene_sets<- death_gene_sets %>% filter(gs_name=="HALLMARK_APOPTOSIS") %>% pull(gene_symbol)
  # seuratobject <- AddModuleScore(seuratobject, name = "apop_score", assay = "SCT", features = list(death_gene_sets), search = TRUE)

  # Rerun PCA and UMAP as well as clustering
  progress$inc(1, detail = 'Re-scaling data...')
  seuratobject <- Seurat::SCTransform(seuratobject, verbose = TRUE)
  progress$inc(1, detail = 'Running PCA...')
  seuratobject <- Seurat::RunPCA(seuratobject, features = Seurat::VariableFeatures(seuratobject))
  progress$inc(1, detail = 'Running UMAP...')
  seuratobject <- Seurat::RunUMAP(seuratobject, reduction = "pca", dims = 1:30)
  progress$inc(1, detail = 'Reclustering...')
  seuratobject <- Seurat::FindNeighbors(seuratobject, reduction = "pca", dims = 1:30)
  seuratobject <- Seurat::FindClusters(seuratobject, resolution = res)
  return(seuratobject)
}


renew_labels <- function(seuratobject, reference) {
  # ref_se <- readRDS('ref_se_020222.rds') # imports mouse immune cell genesets for SingleR


  progress <- shiny::Progress$new(max = 4, min = 0)
  on.exit(progress$close())
  progress$set(message= "Labeling cells...", value = 0)


  counts <- Seurat::GetAssayData(seuratobject)
  progress$inc(1, detail = 'Annotating main labels...')
  cell_annotation <- SingleR::SingleR(counts, ref = reference, assay.type.test = 1, labels = reference$label.main) # performs annotation
  progress$inc(1, detail = 'Finished main annotation...')
  seuratobject <- Seurat::AddMetaData(object = seuratobject, metadata = data.frame(cell_annotation)) # adds metadata of new annotations

  progress$inc(1, detail = 'Annotating fine labels...')
  cell_annotation <- SingleR::SingleR(counts, ref = reference, assay.type.test = 1, labels = reference$label.fine)
  names(cell_annotation) <- BiocGenerics::paste(names(cell_annotation), "_fine", sep = "")
  progress$inc(1, detail = 'Finished fine annotation...')
  seuratobject <- Seurat::AddMetaData(object = seuratobject, metadata = data.frame(cell_annotation)) # adds metadata of new annotations
  seuratobject <- Seurat::SetIdent(seuratobject, value = seuratobject@meta.data$labels) # sets default method to labels

  return(seuratobject)
}

relabel <- function(seuratobject, new_ids) {
  progress <- shiny::Progress$new(max = 1, min = 0)
  on.exit(progress$close())
  progress$set(message= "Relabeling...", value = 0)
  Seurat::Idents(seuratobject) <- seuratobject$seurat_clusters
  seuratobject <- Seurat::RenameIdents(seuratobject, new_ids)
  seuratobject$cell_type <- Seurat::Idents(seuratobject)
  seuratobject$celltype.condition <- BiocGenerics::paste(Seurat::Idents(seuratobject), seuratobject$condition, sep = "_")
  progress$set(message= "Done!", value = 1)
  return(seuratobject)
}

getexp <- function(seuratobject, idents, id1, id2) {
  progress <- shiny::Progress$new(max = 1, min = 0)
  on.exit(progress$close())
  progress$set(message= "Calculating diff. exp. table...", value = 0)
  DET <- Seurat::FindMarkers(seuratobject, group.by = idents, ident.1 = id1, ident.2 = id2, min.pct = 0.25, logfc.threshold = 0.5)
  progress$set(message= "Done!", value = 1)
  return(DET)
}


refine_metadata_levels <- function(seurat_data){
  for (i in base::colnames(seurat_data@meta.data)){
    if (base::is.factor(seurat_data@meta.data[[i]])){
      base::print(base::paste("Re-evaluating levels for a factor column", i))
      base::print(
        base::paste(
          "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
      seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  # need to drop levels of the removed values
      base::print(
        base::paste(
          "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
    }
  }
  return (seurat_data)
}

run_pseudo <- function(seuratobject, root) {
  progress <- shiny::Progress$new(max = 6, min = 0)
  on.exit(progress$close())
  progress$set(message = "Running Pseudotime Analysis...", value = 0, detail = "Converting to cds...")
  cds <- SeuratWrappers::as.cell_data_set(seuratobject)
  progress$inc(1, detail = "Clustering...")
  cds <- monocle3::cluster_cells(cds)
  progress$inc(1, detail = "Learning graph...")
  cds <- monocle3::learn_graph(cds)
  progress$inc(1, detail = "Finding root...")
  max <- BiocGenerics::which.max(BiocGenerics::unlist(SeuratObject::FetchData(seuratobject, root)))
  max <- BiocGenerics::colnames(seuratobject)[max]
  progress$inc(1, detail = "Finishing up...")
  cds <- monocle3::order_cells(cds, root_cells = max)
  return(cds)
}

graph_params <- function(graph_type, object) {
  featurelist <- rownames(object)
  idents <- colnames(object@meta.data)
  params <- shiny::tagList()
  if (graph_type == "ridgeplot") {
    params[[1]] <- shiny::selectInput("gfeature", "Gene name", choices = featurelist, multiple = TRUE)
    params[[2]] <- shiny::selectInput("gcolors", "Colors", choices = colors(), multiple = TRUE)
    params[[3]] <- shiny::selectInput("gidents", "Idents", choices = idents, multiple = FALSE)
    params[[4]] <- shiny::selectInput("gsort", "Sort", choices = c("FALSE", "TRUE"), multiple = FALSE)
    params[[5]] <- shiny::selectInput("gstack", "Stack", choices = c("FALSE", "TRUE"), multiple = FALSE)
  } else if (graph_type == "dotplot") {
    params[[1]] <- shiny::selectInput("gfeature", "Gene name", choices = featurelist, multiple = TRUE)
    params[[2]] <- shiny::selectizeInput("gcolors", "Colors", choices = colors(), selected = c("lightgrey", "blue"), multiple = TRUE, options = list(maxItems = 2))
    params[[3]] <- shiny::selectInput("gidents", "Idents", choices = idents, multiple = FALSE)
    params[[4]] <- shiny::selectInput("gdotscale", "Scale", choices = c("TRUE", "FALSE"), multiple = FALSE)
  } else if (graph_type == "dimplot split") {
    params[[1]] <- shiny::numericInput("gptsize", "Point Size", value = 1, min = 0)
    params[[2]] <- shiny::selectInput("gcolors", "Colors", choices = colors(), multiple = TRUE)
    params[[3]] <- shiny::selectInput("gidents", "Idents", choices = idents, multiple = FALSE)
    params[[4]] <- shiny::selectInput("glabel", "Label", choices = c("FALSE", "TRUE"), multiple = FALSE)
    params[[5]] <- shiny::selectInput("grepel", "Repel", choices = c("FALSE", "TRUE"), multiple = FALSE)
    params[[6]] <- shiny::selectInput("gsplitby", "Split By (leave blank for aggregate)", choices = idents, multiple = FALSE)
  } else if (graph_type == "dimplot") {
    params[[1]] <- shiny::numericInput("gptsize", "Point Size", value = 1, min = 0)
    params[[2]] <- shiny::selectInput("gcolors", "Colors", choices = colors(), multiple = TRUE)
    params[[3]] <- shiny::selectInput("gidents", "Idents", choices = idents, multiple = FALSE)
    params[[4]] <- shiny::selectInput("glabel", "Label", choices = c("FALSE", "TRUE"), multiple = FALSE)
    params[[5]] <- shiny::selectInput("grepel", "Repel", choices = c("FALSE", "TRUE"), multiple = FALSE)
  } else if (graph_type == "heatmap") {
    params[[1]] <- shiny::selectInput("gfeature", "Gene name", choices = featurelist, multiple = TRUE)
    params[[2]] <- shiny::selectInput("gcolors", "Colors", choices = colors(), multiple = TRUE)
    params[[3]] <- shiny::selectInput("gidents", "Idents", choices = idents, multiple = FALSE)
  } else if (graph_type == "featureplot") {
    params[[1]] <- shiny::selectInput("gfeature", "Gene name", choices = c(featurelist, idents), multiple = TRUE)
    params[[2]] <- shiny::selectizeInput("gcolors", "Colors", choices = colors(), selected = c("lightgrey", "blue"), multiple = TRUE, options = list(maxItems = 2))
    params[[3]] <- shiny::numericInput("gptsize", "Point Size", value = 1, min = 0)
    params[[4]] <- shiny::selectizeInput("gidents", "Idents", choices = c(idents), multiple = TRUE, selected = NULL, options = list(maxItems = 1))
    params[[5]] <- shiny::selectInput("glabel", "Label", choices = c("FALSE", "TRUE"), multiple = FALSE)
    params[[6]] <- shiny::selectInput("grepel", "Repel", choices = c("FALSE", "TRUE"), multiple = FALSE)
  } else if (graph_type == "featureplot blend") {
    params[[1]] <- selectizeInput("gfeature", "Features", choices = c(featurelist, idents), multiple = TRUE, selected = NULL, options = list(maxItems = 2))
    params[[2]] <- selectizeInput("gcolors", "Colors", choices = colors(), selected = c("red", "blue"), multiple = TRUE, options = list(maxItems = 2))
    params[[3]] <- numericInput("gptsize", "Point Size", value = 1, min = 0)
    # params[[4]] <- selectizeInput('gidents', 'Idents', choices = c(idents), multiple = TRUE, selected = NULL, options = list(maxItems = 1))
    params[[4]] <- shiny::selectInput("glabel", "Label", choices = c("FALSE", "TRUE"), multiple = FALSE)
    params[[5]] <- shiny::selectInput("grepel", "Repel", choices = c("FALSE", "TRUE"), multiple = FALSE)
  }
  return(params)
}



# output$dynamic_labels <- renderUI({
#  tags <- tagList()
#  for (i in seq_len(length(levels(sobj$obj$seurat_clusters)))){
#    tags[[i]] <- textInput(paste0('n',i),
#                           paste0('Cluster',i-1),
#                           NULL)
#  }
#  tags
# })
options(spinner.image = 'ratdance.gif', spinner.image.height = 100)
options(shiny.maxRequestSize = 100 * 1024^3)
# Define UI for application that draws a histogram
# tags$script(HTML("var header = $('.navbar > .container-fluid'); header.append('<div style=\"float:right\"><a href=\"URL\"><img src=\"logo.png\" alt=\"alt\" style=\"float:right;width:33px;height:41px;padding-top:10px;\"> </a></div>'); console.log(header)"))
ui <- fluidPage(titlePanel(windowTitle = "Stallings Lab Single-Cell RNA Seq Analysis", title = a(img(src = "logo.png", height = 100, width = 100, style = "margin:10px 10px"), href = "http://stallingslab.wustl.edu/", "Stallings Lab Single-Cell RNA Seq Analysis")), theme = shinythemes::shinytheme("united"), navbarPage(" ",
  id = "inTabsetm",
  tabPanel("Data Input",
    shiny::fluidRow(sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel(
            "Load 10x",
            shinyFiles::shinyDirButton('testdir', 'Input directory', "upload"),
            shiny::actionButton("load10x", "Load")
          ),
          tabPanel("Load Loom", shiny::fileInput(
            inputId = "rawcounts",
            label = "Count Matrix",
            multiple = TRUE
          ), shiny::actionButton("loadloom", 'Load')),
          tabPanel(
            "Load Model",
            shiny::fileInput(
              inputId = "model",
              label = "Pre-load model"
            ),
            shiny::actionButton("loadmodel", "Load"),
          )
        )
      ),
      mainPanel(type = "tab", tabsetPanel(
        tabPanel(
          "Data",
          uiOutput('mdatatbl')
          #tableOutput("mdataTbl")
        )
      ))
    ))
  ),
  tabPanel("Merging",
    fluidRow(sidebarLayout(
      sidebarPanel(
        radioButtons("mergetype",
          "Type of merging to use",
          choices = c(
            "Simple Merge",
            "Integrate"
          ),
          selected = character(0)
        )
      ),
      mainPanel(tabsetPanel(tabPanel(
        "Data",
        plotOutput("histcounts") %>% shinycssloaders::withSpinner(),
        plotOutput("histfeatures"),
        plotOutput("histmito"),
        plotOutput("violin"),
        plotOutput("featurescattermt"),
        plotOutput("featurescatterfeatures"),
        plotOutput("Elbow"),
        plotOutput("Dimplot")
      )))
    ))
  ),
  tabPanel("Quality Control",
    fluidRow(sidebarLayout(
      sidebarPanel(
        sliderInput("ncounts", "Cutoffs for number of counts", min = 0, max = 100000, value = c(1000, 7.5e4)),
        sliderInput("nfeatures", "Cutoffs for number of features (genes)", min = 0, max = 25000, value = c(500, 6500)),
        numericInput("percentmt", "Maximum percentage of mitochondrial DNA allowed", value = 5),
        numericInput("res", "Resolution for clustering from 0 to 1", value = 0.5, min = 0, max = 1),
        actionButton("runqc", "Run QC")
      ),
      mainPanel(tabsetPanel(tabPanel(
        "QC'd Clustering",
        plotOutput("Dimplotqc") %>% shinycssloaders::withSpinner(),
        downloadButton("downloadSeur", "Download Seurat Object")
      )))
    ))
  ),
  tabPanel("Expression Exploration",
    fluidRow(sidebarLayout(
      sidebarPanel(tabsetPanel(
        tabPanel(
          "Featureplot",
          uiOutput("featureselect")
        ),
        tabPanel(
          "SingleR",
          selectInput('celldexdb', 'Select reference from Celldex', choices = c('ImmGen','BlueprintENCODE', 'DICE', 'Human Primary Cell Atlas', 'Monaco', 'Mouse RNA-seq','Novershtern'),multiple = FALSE),
          fileInput("ref_se", "Reference for Labels", multiple = FALSE),
          actionButton("runsingler", "Run SingleR Labeling")
        )
      )),
      mainPanel(tabsetPanel(
        tabPanel(
          "Interactive Feature Plot",
          plotOutput("interactiveFeature") %>% shinycssloaders::withSpinner()
        ),
        tabPanel(
          "SingleR Labels Main",
          DT::dataTableOutput("countsbycluster")
        ),
        tabPanel(
          "SingleR Labels Fine",
          DT::dataTableOutput("countsbyclusterfine")
        )
      ))
    ))
  ),
  tabPanel("Cluster Labels",
    fluidRow(sidebarLayout(
      sidebarPanel(
        uiOutput("dynamic_labels"),
        actionButton("relabel", "Re-Label")
      ),
      mainPanel(tabsetPanel(
        tabPanel(
          "singler_dimplot",
          plotOutput("SingleRDimPlot") %>% shinycssloaders::withSpinner(),
          plotOutput("SClusters")
        ),
        tabPanel(
          "DimPlot",
          plotOutput("Dimplotlabelled") %>% shinycssloaders::withSpinner()
        )
      ))
    ))
  ),
  tabPanel("Differential Expression",
    fluidRow(sidebarLayout(
      sidebarPanel(tabsetPanel(
        tabPanel(
          "Comparison",
          uiOutput("ident"),
          uiOutput("comparison")
        ),
        tabPanel(
          "EnrichR",
          selectInput("dbs", "EnrichR Databases to use", choices = BiocGenerics::sort(enrichR::listEnrichrDbs()$libraryName), multiple = TRUE),
          actionButton("enrich", "Run EnrichR")
        )
      )),
      mainPanel(tabsetPanel(
        tabPanel(
          "DE table",
          DT::dataTableOutput("diffexp") %>% shinycssloaders::withSpinner(),
          downloadButton("downloaddiffexp", "Download Differential Expression Table")
        ),
        tabPanel(
          "Positive Enrichment",
          uiOutput("penrichplots")
        ),
        tabPanel(
          "Negative Enrichment",
          uiOutput("nenrichplots")
        )
      ))
    ))
  ),
  tabPanel("Cell Cycle Analysis",
    fluidRow(sidebarLayout(
      sidebarPanel(
        actionButton("ccscore", " Run Analysis")
      ),
      mainPanel(tabsetPanel(
        tabPanel(
          "Total",
          plotOutput("ccdimplot") %>% shinycssloaders::withSpinner()
        ),
        tabPanel(
          "By Condition",
          plotOutput("ccdimplotsplit") %>% shinycssloaders::withSpinner()
        )
      ))
    ))
  ),
  tabPanel("Module Score",
    fluidRow(sidebarLayout(
      sidebarPanel(tabsetPanel(
        tabPanel(
          "Custom Score",
          fileInput("geneset1", "First Geneset"),
          fileInput("geneset2", "Second Geneset"),
          textInput("modulename", "Module Name"),
          actionButton("addmod", "Add Module")
        ),
        tabPanel(
          "Visualize Module",
          uiOutput("modfeature")
        )
      )),
      mainPanel(tabsetPanel(
        tabPanel(
          "Module Stats",
          DT::dataTableOutput("modstats")
        ),
        tabPanel(
          "Module Graph",
          plotOutput("modfplot") %>% shinycssloaders::withSpinner() # ,
          # downloadButton('downloadFplot', 'Download Plot')
        )
      ))
    ))
  ),
  tabPanel("Subset",
    fluidRow(sidebarLayout(
      sidebarPanel(
        uiOutput('subsetident'),
        uiOutput('subsetfeatures'),
        selectInput('subreclust', 'Recluster', choices = c('TRUE', 'FALSE'), multiple = FALSE, selected = 'FALSE'),
        uiOutput('subclustres'),
        actionButton('runsubset', 'Subset'),
        actionButton('restore','Restore Parent')
      ),
      mainPanel(
        plotOutput("subsetdimplot")
      )
    ))
  ),
  tabPanel("Pseudotime",
    fluidRow(sidebarLayout(
      sidebarPanel(
        uiOutput('pseudoroot'),
        actionButton('runpseudo', 'Run Pseudotime')
      ),
      mainPanel(
        plotOutput('pseudoplot') %>% shinycssloaders::withSpinner(),
        downloadButton("downloadpseudo", "Download")
      )
    ))
  ),
  tabPanel("Graph Generator",
    fluidRow(sidebarLayout(
      sidebarPanel(tabsetPanel(
        tabPanel(
          "Graph options",
          selectInput("graphtomake", label = "Graph Type:", choices = c("ridgeplot", "dimplot", "dimplot split", "dotplot", "heatmap", "featureplot", "featureplot blend"), multiple = FALSE),
          uiOutput("graphparams"),
          actionButton("graphbutton", label = "Graph")
        ),
        tabPanel(
          "Download Options",
          selectInput("gunits", "Units", choices = c("in", "cm", "mm", "px"), selected = "in"),
          numericInput("gwidth", "Width", value = 7),
          numericInput("gheight", "Height", value = 7),
          numericInput("gscale", "Scale", value = 1),
          selectInput("gmode", "File type", choices = c(".png", ".svg", ".jpg"), selected = ".png")
        )
      )),
      mainPanel(
        plotOutput("graph") %>% shinycssloaders::withSpinner(),
        downloadButton("downloadgraph", "Download")
      )
    ))
  )
))




# Define server logic required to draw a histogram
server <- function(input, output) {
  raw_inputs <- reactive({
    if (is.null(input$rawcounts)&(input$testdir == 0)&is.null(input$model)) {
      return(NULL)
    }
    else if (!(is.null(input$rawcounts)) | !(is.null(input$model))) {
      return(input$rawcounts)
    }
    else if (!(input$testdir==0)) {
      return(parseFilePaths(volumes, input$testdir))
    }
  })
  output$mdataTbl <- renderTable({
    raw_inputs()
  })
  output$mdatatbl <- renderUI({
    if (input$load10x == 1) {
      output$testing <- renderText(subdirs())
      textOutput("testing")
    }
    else if (!(is.null(input$rawcounts)) | !(is.null(input$model))){
      output$loomin <- renderTable({raw_inputs()})
      tableOutput('loomin')
    }
  })
  sobj <- reactiveValues(obj = NULL, tag = "raw")

  ftype <- reactiveValues(read10x = NULL)

  observeEvent(input$loadmodel, {
    req(input$model)
    progress <- shiny::Progress$new(max = 1, min = 0)
    on.exit(progress$close())
    progress$set(message= "Loading object...", value = 0)
    sobj$obj <- readRDS(input$model$datapath)
    progress$inc(1, message= "Done!")
    Sys.sleep(0.5)
  })

  observeEvent(input$loadloom, {
    ftype$read10x <- FALSE
  })

  observeEvent(input$load10x, {
    ftype$read10x <- TRUE
  })

  observe({
    req(input$mergetype, !(is.null(ftype$read10x)))
    if (!is.null(sobj$obj)) {
    } else if (input$mergetype == "Simple Merge") {
      
        if (ftype$read10x) {
          req(input$testdir)
          sobj$obj <- renew_merged(names(), subdirs(), read10x = TRUE)}
        else {
          req(input$rawcounts)
          sobj$obj <- renew_merged(raw_inputs()$name, raw_inputs()$datapath, read10x = FALSE)}
    } else if (input$mergetype == "Integrate") {
      if (ftype$read10x) {
        req(input$testdir)
        sobj$obj <- renew_integrated(names(), subdirs(), read10x = TRUE)}
      else {
        req(input$rawcounts)
        sobj$obj <- renew_integrated(raw_inputs()$name, raw_inputs()$datapath, read10x = FALSE)}
    }

    sobj$tag <- "merged"
    output$violin <- renderPlot({
      req(input$mergetype)
      Seurat::VlnPlot(sobj$obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    })

    output$histfeatures <- renderPlot({
      req(input$mergetype)
      hist(sobj$obj$nCount_RNA)
    })

    output$histmito <- renderPlot({
      req(input$mergetype)
      hist(sobj$obj$nFeature_RNA)
    })

    output$histcounts <- renderPlot({
      req(input$mergetype)
      hist(sobj$obj$percent.mt)
    })

    output$featurescattermt <- renderPlot({
      req(input$mergetype)
      Seurat::FeatureScatter(sobj$obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    })

    output$featurescatterfeatures <- renderPlot({
      req(input$mergetype)
      Seurat::FeatureScatter(sobj$obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    })

    output$Elbow <- renderPlot({
      req(input$mergetype)
      Seurat::ElbowPlot(sobj$obj, ndims = 40, reduction = "pca")
    })

    output$Dimplot <- renderPlot({
      req(input$mergetype)
      Seurat::DimPlot(sobj$obj, group.by = "condition")
    })
  })

  # seuratObjQC <- reactive({

  #  if (is.null(input$model)){
  #    req(seuratObj())
  #    qc(seuratObj(), featuremin = input$nfeatures[1], featuremax = input$nfeatures[2], countmin = input$ncounts[1], countmax = input$ncounts[2], percentmt = input$percentmt, res = input$res)
  #  }
  #  else {
  #    req(input$model)
  #    readRDS(input$model$datapath)
  #  }
  # })

  observeEvent(input$runqc, {
    if (is.null(input$model)) {
      sobj$obj <- qc(sobj$obj, featuremin = input$nfeatures[1], featuremax = input$nfeatures[2], countmin = input$ncounts[1], countmax = input$ncounts[2], percentmt = input$percentmt, res = input$res)
    }
    output$Dimplotqc <- renderPlot({
      req(input$runqc)
      Seurat::DimPlot(sobj$obj, group.by = "seurat_clusters")
    })
    sobj$tag <- "qc"
  })



  output$downloadSeur <- downloadHandler(
    filename = function() {
      BiocGenerics::paste("SeuratObjectQC", ".rds", sep = "")
    },
    content = function(con) {
      saveRDS(sobj$obj, con)
    }
  )

  output$interactiveFeature <- renderPlot({
    Seurat::FeaturePlot(sobj$obj, features = input$feature, raster = FALSE)
  })

  featurelist <- reactive({
    BiocGenerics::sort(BiocGenerics::rownames(sobj$obj))
  })

  output$featureselect <- renderUI({
    selectInput("feature", "Select Feature", choices = featurelist())
  })

  ref_se <- reactive({
    if (!is.null(input$ref_se)) {
      return(readRDS(input$ref_se$datapath))
    }
    if (!is.null(input$celldexdb)) {
      ifelse(input$celldexdb == 'ImmGen', return(celldex::ImmGenData()),
             ifelse(input$celldexdb =='BlueprintENCODE', return(celldex::BlueprintEncodeData()),
                    ifelse(input$celldexdb == 'DICE', return(celldex::DatabaseImmuneCellExpressionData()),
                          ifelse(input$celldexdb == 'Human Primart Cell Atlas', return(celldex::HumanPrimaryCellAtlasData()),
                                 ifelse(input$celldexdb == 'Monaco', return(celldex::MonacoImmuneData()),
                                        ifelse(input$celldexdb == 'Mouse RNA-seq', return(celldex::MouseRNAseqData()),
                                            return(celldex::NovershternHematopoieticData())))))))
    }
  })
  # seuratObjQClabelled <- reactive({
  # req(seuratObjQC, ref_se())
  # out <- renew_labels(seuratObjQC(), ref_se())
  # if (!is.null(new_ids)){
  # Idents(out) <- out$seurat_clusters
  # out <- RenameIdents(out, new_ids)
  # out$cell_type <- Idents(out)
  ## }
  # out
  # })

  observeEvent(input$runsingler, {
    req(ref_se())
    sobj$obj <- renew_labels(sobj$obj, ref_se())
    output$SingleRDimPlot <- renderPlot({
      Seurat::DimPlot(sobj$obj, group.by = "labels", label = TRUE, repel = TRUE)
    })
    sobj$tag <- "labeled"
    output$SClusters <- renderPlot({
      Seurat::DimPlot(sobj$obj, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
    })
  })

  observeEvent(input$relabel, {
    ids <- c()
    for (i in seq_len(length(levels(sobj$obj$seurat_clusters)))) {
      ids <- c(ids, input[[paste0("n", i)]])
    }
    names(x = ids) <- levels(x = sobj$obj$seurat_clusters)
    sobj$obj <- relabel(sobj$obj, ids)
    output$Dimplotlabelled <- renderPlot({
      req(input$relabel)
      Seurat::DimPlot(sobj$obj, label = TRUE, repel = TRUE)
    })

    sobj$tag <- "relabeled"
  })

  counts_by_cluster <- reactive({
    table(
      sobj$obj@meta.data$labels,
      sobj$obj@meta.data$seurat_clusters
    )
  })
  counts_by_cluster_fine <- reactive({
    table(
      sobj$obj@meta.data$labels_fine,
      sobj$obj@meta.data$seurat_clusters
    )
  })
  output$countsbycluster <- DT::renderDataTable({
    DT::datatable(
      as.data.frame.matrix(counts_by_cluster()),
      options = list(scrollX = TRUE, pageLength = 100)
    )
  })
  output$countsbyclusterfine <- DT::renderDataTable({
    DT::datatable(
      as.data.frame.matrix(counts_by_cluster_fine()),
      options = list(scrollX = TRUE, pageLength = 100)
    )
  })

  output$dynamic_labels <- renderUI({
    tags <- tagList()
    for (i in seq_len(length(levels(sobj$obj$seurat_clusters)))) {
      tags[[i]] <- textInput(
        paste0("n", i),
        paste0("Cluster", i - 1),
        NULL
      )
    }
    tags
  })

  # observeEvent(input$relabel,{
  # req(new_ids())
  # Idents(seuratObjQClabelled()) <- seuratObjQClabelled()$seurat_clusters
  # seuratObjQClabelled() <- RenameIdents(seuratObjQClabelled(), new_ids())
  # seuratObjQClabelled()$cell_type <- Idents(seuratObjQClabelled())
  # print(seuratObjQClabelled()$cell_type)
  # })

  # new_ids <- eventReactive(input$relabel, {
  #  ids <- c()
  #  for (i in seq_len(length(levels(seuratObjQClabelled()$seurat_clusters)))){
  #    ids <- c(ids, input[[paste0('n',i)]])
  #  }
  #  names(x = ids) <- levels(x = seuratObjQClabelled()$seurat_clusters)
  #  ids
  # })

  # seuratObjQCrelabelled <- reactive({
  #  return(relabel(seuratObjQClabelled(), new_ids()))
  # })

  identlist <- reactive({
    BiocGenerics::colnames(sobj$obj@meta.data)
  })

  output$ident <- renderUI({
    selectizeInput("idents", "Select Idents", choices = identlist(), selected = NULL, options = list(maxItems = 1))
  })

  output$comparison <- renderUI({
    selectizeInput("comp", "Comparison to Make", choices = levels(as.factor(BiocGenerics::unique(sobj$obj@meta.data[input$idents])[[input$idents]])), selected = NULL, options = list(maxItems = 2))
  })

  DE <- reactive({
    req(input$comp[1], input$comp[2])
    Seurat::FindMarkers(sobj$obj,
      ident.1 = input$comp[1],
      ident.2 = input$comp[2], group.by = input$idents, verbose = FALSE,
      min.pct = 0.25,
      logfc.threshold = 0.5
    )
  })


  pos <- reactive({
    req(DE())
    DE()[sign(DE()$avg_log2FC) == 1, ]
  })

  neg <- reactive({
    req(DE())
    DE()[sign(DE()$avg_log2FC) == -1, ]
  })
  output$diffexp <- DT::renderDataTable({
    DT::datatable(
      as.data.frame.matrix(DE()),
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })

  output$downloaddiffexp <- downloadHandler(
    filename = function() {
      BiocGenerics::paste(input$comp[1], "_vs_", input$comp[2], ".csv", sep = "")
    },
    content = function(con) {
      utils::write.csv(DE(), con)
    }
  )


  posenrichment <- reactive({
    req(input$enrich)
    progress <- shiny::Progress$new(max = 1, min = 0)
    on.exit(progress$close())
    progress$set(message= "Positive Enrichment...", value = 1)
    enrichR::enrichr(BiocGenerics::rownames(pos()), databases = input$dbs)
  })

  output$penrichplots <- renderUI({
    plot_output_list <- BiocGenerics::lapply(1:length(input$dbs), function(i) {
      plotname <- BiocGenerics::paste("pplot", i, sep = "")
      plotOutput(plotname) %>% shinycssloaders::withSpinner()
    })
    BiocGenerics::do.call(tagList, plot_output_list)
  })

  observeEvent(input$enrich, {
    BiocGenerics::lapply(1:length(input$dbs), function(i) {
      local({
        my_i <- i
        plotname <- BiocGenerics::paste("pplot", i, sep = "")

        output[[plotname]] <- renderPlot({
          enrichR::plotEnrich(posenrichment()[[i]], title = BiocGenerics::paste("Positively Enriched Pathways in ", input$dbs[i], sep = ""))
        })
      })
    })
  })
  negenrichment <- reactive({
    req(input$enrich)
    progress <- shiny::Progress$new(max = 1, min = 0)
    on.exit(progress$close())
    progress$set(message= "Negative Enrichment...", value = 1)
    enrichR::enrichr(BiocGenerics::rownames(neg()), databases = input$dbs)
  })

  output$nenrichplots <- renderUI({
    plot_output_list <- BiocGenerics::lapply(1:length(input$dbs), function(i) {
      plotname <- BiocGenerics::paste("nplot", i, sep = "")
      plotOutput(plotname) %>% shinycssloaders::withSpinner()
    })
    BiocGenerics::do.call(tagList, plot_output_list)
  })

  observeEvent(input$enrich, {
    BiocGenerics::lapply(1:length(input$dbs), function(i) {
      local({
        my_i <- i
        plotname <- BiocGenerics::paste("nplot", i, sep = "")

        output[[plotname]] <- renderPlot({
          enrichR::plotEnrich(negenrichment()[[i]], title = BiocGenerics::paste("Negatively Enriched Pathways in ", input$dbs[i], sep = ""))
        })
      })
    })
  })

  observeEvent(input$ccscore, {
    req(sobj$obj)
    progress <- shiny::Progress$new(max = 1, min = 0)
    on.exit(progress$close())
    progress$set(message= "Running Cell Cycle Scoring...", value = 0)
    sobj$obj <- Seurat::CellCycleScoring(object = sobj$obj, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
    progress$inc(1, message = "Done!")
  })

  output$ccdimplot <- renderPlot({
    req(sobj$obj, input$ccscore)
    Seurat::DimPlot(sobj$obj, group.by = "Phase")
  })

  output$ccdimplotsplit <- renderPlot({
    req(sobj$obj, input$ccscore)
    Seurat::DimPlot(sobj$obj, group.by = "Phase", split.by = "condition")
  })

  observeEvent(input$addmod, {
    req(input$geneset1, input$geneset2, input$modulename)
    progress <- shiny::Progress$new(max = 1, min = 0)
    on.exit(progress$close())
    progress$set(message= "Calculating Scores...", value = 0)
    gset1 <- readr::read_csv(input$geneset1$datapath)$x
    gset2 <- readr::read_csv(input$geneset2$datapath)$x
    modname <- input$modulename
    gset1 <- sobj$obj@assays[["SCT"]]@var.features[sobj$obj@assays[["SCT"]]@var.features %in% gset1]
    gset2 <- sobj$obj@assays[["SCT"]]@var.features[sobj$obj@assays[["SCT"]]@var.features %in% gset2]
    Module_features <- list("1" = gset1, "2" = gset2)
    sobj$obj <- Seurat::AddModuleScore(sobj$obj, name = modname, assay = Seurat::DefaultAssay(sobj$obj), features = Module_features, search = TRUE)
    progress$inc(1, message = "Done!")
  })

  output$modfeature <- renderUI({
    selectizeInput("idents2", "Select TWO Idents", choices = identlist(), selected = NULL, options = list(maxItems = 2))
  })

  fplot <- reactive({

  })

  output$modfplot <- renderPlot({
    req(input$idents2[2])
    print("Making Plot")
    print(input$idents2[1])
    print(input$idents2[2])
    print(c(BiocGenerics::paste(input$idents2[1]), BiocGenerics::paste(input$idents2[2])))
    FeaturePlot(sobj$obj, features = c(BiocGenerics::paste(input$idents2[1]), BiocGenerics::paste(input$idents2[2])), blend = TRUE, label = FALSE, repel = TRUE, pt.size = 3)
  })

  output$downloadFplot <- downloadHandler(
    filename = function() {
      BiocGenerics::paste("Feature_plot", ".png", sep = "")
    },
    content = function(file) {
      ggplot2::ggsave(file, plot = fplot(), width = 45, height = 15, units = "in", scale = 0.5)
    }
  )


  output$subsetident <- renderUI({
    selectizeInput("subsetidents", "Select Ident to subset on", choices = identlist(), selected = NULL, options = list(maxItems = 1))
  })

  output$subsetfeatures <- renderUI({
    selectizeInput("subsetvals", "Select which values to subset", choices = levels(as.factor(BiocGenerics::unique(sobj$obj@meta.data[input$subsetidents])[[input$subsetidents]])), selected = NULL, multiple = TRUE)
  })

  observeEvent(input$runsubset, {
    print(input$subsetidents)
    print(input$subsetvals)
    progress <- shiny::Progress$new(max = 1, min = 0)
    on.exit(progress$close())
    progress$set(message= "Saving current object...", value = 0)
    sobj$prev <- sobj$obj
    progress$inc(1, message = "Subsetting object...")
    Seurat::Idents(sobj$obj) <- sobj$obj[[input$subsetidents]]
    sobj$obj <- subset(sobj$obj, idents = input$subsetvals)
    sobj$obj <- refine_metadata_levels(sobj$obj)
    progress$close()
    if (input$subreclust == 'TRUE') {
      progress <- shiny::Progress$new(max = 5, min = 0)
      on.exit(progress$close())
      progress$set(message= "Re-clustering...", value = 0)
      progress$inc(1, detail = "Rescaling...")
      sobj$obj <- Seurat::SCTransform(sobj$obj)
      progress$inc(1, detail = "Running PCA...")
      sobj$obj <- Seurat::RunPCA(sobj$obj)
      progress$inc(1, detail = "Running UMAP...")
      sobj$obj <- Seurat::RunUMAP(sobj$obj, reduction = 'pca', dims = 1:30)
      progress$inc(1, detail = "Finding neighbors...")
      sobj$obj <- Seurat::FindNeighbors(sobj$obj)
      progress$inc(1, detail = "Clustering...")
      sobj$obj <- Seurat::FindClusters(sobj$obj, resolution = input$subres)
    }
  })

  output$subclustres <- renderUI({
    if (input$subreclust == 'TRUE'){
      output_rename = tagList()
      output_rename[[1]] = numericInput("subres", "Resolution for clustering from 0 to 1", value = 0.5, min = 0, max = 1)
      output_rename
    }
  })

  observeEvent(input$restore, {
    progress <- shiny::Progress$new(max = 1, min = 0)
    on.exit(progress$close())
    progress$set(message= "Saving current object...", value = 0)
    temp <- sobj$obj
    progress$inc(1, message = "Restoring parent...")
    sobj$obj <- sobj$prev
    sobj$prev <- temp
    rm(temp)
  })

  output$subsetdimplot <- renderPlot({
    Seurat::DimPlot(sobj$obj)
  })

  output$pseudoroot <- renderUI({
    selectInput('pseudogene', label = 'Select gene to use as root.', choices = featurelist(), multiple = FALSE)
  })

  pplot <- eventReactive(input$runpseudo, {
    print('Running pseudo')
    cds <- run_pseudo(sobj$obj, input$pseudogene)
    monocle3::plot_cells(cds, color_cells_by = 'pseudotime', label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + ggplot2::ggtitle(paste("Pseudotime map with max ",input$pseudogene, "expression as root"))
  })

  output$pseudoplot <- renderPlot({pplot()})

  output$graphparams <- renderUI({
    req(input$graphtomake, sobj$obj)
    params <- graph_params(input$graphtomake, sobj$obj)
    params
  })

  output$downloadpseudo <- downloadHandler(
    filename = "Pseudotime.svg",
    content = function(file) {
      ggplot2::ggsave(file, plot = pplot())
    }
  )

  graph <- eventReactive(input$graphbutton, {
    req(sobj$obj, input$graphtomake)
    if (!(is.null(input$gidents))) {
      Seurat::Idents(sobj$obj) <- input$gidents
    }
    if (input$graphtomake == "ridgeplot") {
      Seurat::RidgePlot(sobj$obj, features = input$gfeature, cols = input$gcolors, sort = input$gsort, stack = input$gstack)
    } else if (input$graphtomake == "dotplot") {
      if (length(input$gfeature) > 1) {
        Seurat::DotPlot(sobj$obj, features = input$gfeature, cols = input$gcolors, scale = input$gdotscale)
      }
    } else if (input$graphtomake == "dimplot") {
      Seurat::DimPlot(sobj$obj, pt.size = input$gptsize, cols = input$gcolors, group.by = input$gidents, label = input$glabel, repel = input$grepel)
    } else if (input$graphtomake == "dimplot split") {
      Seurat::DimPlot(sobj$obj, pt.size = input$gptsize, cols = input$gcolors, group.by = input$gidents, label = input$glabel, repel = input$grepel, split.by = input$gsplitby)
    } else if (input$graphtomake == "heatmap") {
      Seurat::DoHeatmap(sobj$obj, features = input$gfeature, group.colors = input$gcolors, group.by = input$gidents)
    } else if (input$graphtomake == "featureplot") {
      Seurat::FeaturePlot(sobj$obj, features = input$gfeature, cols = input$gcolors, pt.size = input$gptsize, shape.by = input$gidents, label = input$glabel, repel = input$grepel)
    } else if (input$graphtomake == "featureplot blend") {
      Seurat::FeaturePlot(sobj$obj, features = input$gfeature, cols = input$gcolors, pt.size = input$gptsize, label = input$glabel, repel = input$grepel, blend = TRUE)
    }
  })
  output$graph <- renderPlot({
    req(input$graphbutton)
    graph()
  })

  output$downloadgraph <- downloadHandler(
    filename = function() {
      BiocGenerics::paste(input$graphtomake, input$gmode, sep = "")
    },
    content = function(file) {
      ggplot2::ggsave(file, plot = graph(), width = input$gwidth, height = input$gheight, units = input$gunits, scale = input$gscale)
    }
  )

  volumes <- shinyFiles::getVolumes()()
  shinyFiles::shinyDirChoose(
    input,
    'testdir',
    #roots = c(home = '~')
    roots = volumes
  )
  dirname <- reactive({shinyFiles::parseDirPath(volumes, input$testdir)})
  subdirs <- reactive({list.dirs(dirname(), recursive = FALSE)})
  names <- reactive ({list.dirs(dirname(), full.names = FALSE, recursive = FALSE)})
  observe({
    if(!is.null(dirname)){
      #print(dirname())
      #print(subdirs())
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
