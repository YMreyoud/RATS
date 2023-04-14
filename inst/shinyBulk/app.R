#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#library(shiny)
#library(DT)
#library(tidyverse)
#library(readxl)
#library(edgeR)
#library(limma)
#library(ggrepel)
#library(pheatmap)

`%>%` <- magrittr::`%>%`

contrastsvector <- c()

fitlimma <- function(dge, design.mat, condition) {
  logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  colnames(logCPM) <- condition
  fit <- limma::lmFit(logCPM, design.mat)
  fit <- limma::eBayes(fit)
  return(fit)
}

fitvoom <- function(dge, design.mat) {
  v <- limma::voomWithQualityWeights(dge, design.mat, plot = TRUE)
  vfit <- limma::lmFit(v, design.mat)
  vfit <- limma::eBayes(vfit)
  return(vfit)
}

save_pheatmap <- function(hm, filename, width=7, height=7){
  filename <- paste(filename, '.pdf', sep = '')
  pdf(filename, width=width, height=height)
  print(hm)
  dev.off()
}

graph_params <- function(graph_type, object, model = NULL) {
  genes <- object$genes$genes
  contrast <- colnames(model$contrasts)
  params <- tagList()
  if (graph_type == "volcano") {
    params[[1]] <- numericInput("gpvalue", "P value cut-off", value = 0.05, min = 0, max = 1)
    params[[2]] <- selectInput("ggenes", "Genes", choices = genes, multiple = TRUE)
    params[[3]] <- fileInput("ggeneset", "Optional Geneset", multiple = FALSE)
    params[[4]] <- selectInput("gcontrast", "Contrast", choices = contrast, multiple = FALSE)
    params[[5]] <- selectInput("glabel", "Label", choices = c("FALSE", "TRUE"), multiple = FALSE)
  } else if (graph_type == "heatmap") {
    params[[1]] <- selectInput("ggenes", "Genes", choices = genes, multiple = TRUE)
    params[[2]] <- selectizeInput("gcondition", "Conditions", choices = object$samples$group, multiple = TRUE)
    #params[[3]] <- selectInput("gidents", "Idents", choices = idents, multiple = FALSE)
    #params[[4]] <- selectInput("gdotscale", "Scale", choices = c("TRUE", "FALSE"), multiple = FALSE)
  }
  return(params)
}


options(shiny.maxRequestSize = 100 * 1024^3)

# Define UI for application that draws a histogram
ui <-  fluidPage(titlePanel(windowTitle = "Stallings Lab Bulk RNA Seq Analysis", title = a(img(src = "logo.png", height = 100, width = 100, style = "margin:10px 10px"), href = "http://stallingslab.wustl.edu/", "Stallings Lab Bulk RNA Seq Analysis")),  navbarPage(" ",
  id = "inTabsetm",
  tabPanel("Data Input",
    shiny::fluidRow(sidebarLayout(
      sidebarPanel(
          shiny::fileInput(
            inputId = "rawcounts",
            label = "Count Matrix",
            multiple = FALSE
          ),
          shiny::radioButtons('extension', 'File extension', choices = c('tsv','csv')),
          shiny::actionButton("loaddata", "Upload")
      ),
      mainPanel(type = "tab", tabsetPanel(
        tabPanel(
          "Data",
          #tableOutput("mdataTbl"),
          DT::dataTableOutput("counts")
        )
      ))
    ))
  ),
  tabPanel("Data Preparation",
    shiny::fluidRow(sidebarLayout(
      sidebarPanel(
        uiOutput('selecting'),
        #selectInput('genes', 'Select the column containing Gene Name information'),
        #selectInput('count', 'Select the count columns you would like to use in the analysis', multiple = TRUE),
        #actionButton('selection1', 'Done'),
        uiOutput('renaming'),
        #textInput('genes', 'List column number containing gene names (e.g 1)'),
        #textInput('count', 'List column number of samples you would like to analyze separated by commas (e.g 2,3)'),
        #textInput('group', 'List the label for each sample selected separated by commas (e.g KO_Naive,WT_Naive)'),
        #actionButton('selection', 'Done')
      ),
      mainPanel(tabsetPanel(
        tabPanel(
          "Raw Data",
          DT::dataTableOutput("counts2")
        ),
        tabPanel(
          "Data Stats",
          plotOutput('bplot'),
          plotOutput('mds')
        )
      )
      )
    ))
  ),
  tabPanel("Model Fitting",
    shiny::fluidRow(sidebarLayout(
      sidebarPanel(
            selectInput('model', 'Select model to use', c('Limma','Voom', selected = 'Voom')),
            actionButton('fit', 'Fit')
      ),
      mainPanel(
        fluidRow( splitLayout(cellWidths = c("50%", "50%"), plotOutput("SA"), plotOutput("MD")),
        plotOutput('vplot')
      ))
    ))),
  tabPanel("Making Contrasts",
           fluidRow(sidebarLayout(
             sidebarPanel(
               uiOutput('selectcon'),
               actionButton('add', 'Add Contrast'),
               actionButton('clear','Clear'),
               actionButton('done', 'Done')
             ),
             mainPanel(
               textOutput('currcons')
             )
           ))),
  tabPanel("DGE",
           fluidRow(sidebarLayout(
             sidebarPanel(
               uiOutput('comparison'),
               #numericInput('n', 'Top N DGEs to show', value = Inf),
               selectInput('adjust', 'Adjust Method', choices = c('none', 'BH', 'BY', 'holm'), selected = 'BH'),
               selectInput('sortby', 'Sort By', choices = c('logFC', 'AveExpr', 'T','P','B','none'), selected = 'P'),
               numericInput('pcutoff', 'P-value cutoff', value = 0.05, min = 0, max = 1),
               actionButton('rundge', 'Done')
             ),
             mainPanel(
               tabPanel("DE output",
               DT::dataTableOutput("DEGS"),
               uiOutput('downloaddge')
             ))
           ))),
  tabPanel("Graphs",
           fluidRow(sidebarLayout(
             sidebarPanel(tabsetPanel(
               tabPanel(
                 "Graph options",
                 selectInput("graphtomake", label = "Graph Type:", choices = c("volcano", "heatmap"), multiple = FALSE),
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
               plotOutput('graph'),
               downloadButton("downloadgraph", "Download")
             )
           )))
))

# Define server logic required to draw a histogram
server <- function(input, output) {

    raw_counts <- reactiveValues(obj = NULL, dge = NULL, design.mat = NULL, model = NULL, tag = "raw", contrasts = c(), contrastmat = NULL, DEGS = NULL)
    raw_inputs <- reactive({
        if (is.null(input$rawcounts)) {
            return(NULL)
        }
        return(input$rawcounts)
    })
    output$mdataTbl <- renderTable({
        raw_inputs()
    })
    sep <- reactive({
      if (input$extension == 'tsv') {
        '\t'
      }
      else {
        ','
      }
    })
    observeEvent(input$loaddata, {
        raw_counts$obj <- utils::read.csv(input$rawcounts$datapath, sep = sep())
        output$counts <- DT::renderDataTable({
          DT::datatable(
            as.data.frame.matrix(raw_counts$obj),
            options = list(scrollX = TRUE, pageLength = 10)
          )
        })
        output$counts2 <- DT::renderDataTable({
          DT::datatable(
            as.data.frame.matrix(raw_counts$obj),
            options = list(scrollX = TRUE, pageLength = 10)
          )
        })

    })
    output$selecting <- renderUI({
      columns = colnames(raw_counts$obj)
      output = tagList()
      output[[1]] = selectInput('genes', 'Select the column containing Gene Name information', choices = columns)
      output[[2]] = textInput('conditions', 'List all conditions used, separated by commas (e.g. KO_Naive,WT_Naive,...).')
      output[[3]] = selectInput('count', 'Select the count columns you would like to use in the analysis', multiple = TRUE, choices = columns)

      #output[[3]] = actionButton('selection1', 'Done')
      output
    })
    output$renaming <- renderUI({
      groups <- strsplit(stringr::str_replace_all(input$conditions, ' ', ''), ',')[[1]]
      if (!is.null(input$count)){
        output_rename = tagList()
        for (i in seq_len(length(input$count))) {
          output_rename[[i]] = selectInput(
            paste0('col',input$count[[i]]),
            input$count[[i]],
            choices = groups
            )
        }
        output_rename[[length(input$count)+1]] = actionButton('selection', 'Done')
        output_rename
      }
    })
    observeEvent(input$selection, {
      selected <- input$count
      genes <- input$genes
      #selected <- as.numeric(strsplit(input$count, ',')[[1]])
      #genes <- as.numeric(input$genes)
      #groups <- strsplit(str_replace_all(input$group, ' ', ''), ',')[[1]]
      groups <- c()
      for (i in seq_len(length(input$count))) {
        groups <- c(groups, input[[paste0("col", input$count[[i]])]])
      }
      raw_counts$dge <- edgeR::DGEList(counts = raw_counts$obj[,selected], genes = raw_counts$obj[,genes],
                                group = groups) %>% edgeR::calcNormFactors(method = "TMM")
      raw_counts$design.mat <- stats::model.matrix(~ 0 + raw_counts$dge$samples$group)
      colnames(raw_counts$design.mat) <- levels(raw_counts$dge$samples$group)
      keep <- edgeR::filterByExpr(raw_counts$dge, design = raw_counts$design.mat)
      raw_counts$dge <- raw_counts$dge[keep,]
      output$bplot <- renderPlot({
        barplot(raw_counts$dge$samples$lib.size, names=colnames(raw_counts$dge),las=2)
      })
      output$mds <- renderPlot({
        colConditions <- rainbow(length(levels(raw_counts$dge$samples$group)))
        colConditions <- colConditions[BiocGenerics::match((raw_counts$dge$samples$group), levels(raw_counts$dge$samples$group))]
        limma::plotMDS(raw_counts$dge, col = colConditions, labels = raw_counts$dge$samples$group)
      })

    })
    observeEvent(input$fit, {
      if(input$model == 'Limma'){
          raw_counts$model <- fitlimma(raw_counts$dge, raw_counts$design.mat, raw_counts$dge$samples$group)
      }
      else {
        output$vplot <- renderPlot({
          raw_counts$model <- fitvoom(raw_counts$dge, raw_counts$design.mat)
        })
      }
      output$SA <- renderPlot({
        limma::plotSA(raw_counts$model)
      })
      output$MD <- renderPlot({
        limma::plotMD(raw_counts$model)
      })
    })
    output$selectcon <- renderUI({
      selectizeInput("selected", "Select Contrast", choices = levels(raw_counts$dge$samples$group), selected = NULL, options = list(maxItems = 2))
    })
    observeEvent(input$add, {
      raw_counts$contrasts <- BiocGenerics::append(raw_counts$contrasts, stringr::str_c(input$selected[1],'-',input$selected[2]))
    })
    observeEvent(input$clear, {
      raw_counts$contrasts <- c()
      #output$currcons <- renderText({
      #  paste('Current selections:')
      #})
    })
    output$currcons <- renderText({
      BiocGenerics::paste(BiocGenerics::append(raw_counts$contrasts, 'Current selections:', 0), collapse = ", ")
    })
    observeEvent(input$done, {
      if (purrr::is_empty(raw_counts$contrasts)) return()
      raw_counts$contrastmat <- limma::makeContrasts(contrasts = raw_counts$contrasts, levels = raw_counts$design.mat)
      raw_counts$model <- limma::contrasts.fit(raw_counts$model, raw_counts$contrastmat)
      raw_counts$model <- limma::eBayes(raw_counts$model)
    })
    output$comparison <- renderUI({
      selectizeInput("comparisons", "Select Comparison", choices = colnames(raw_counts$model$contrasts), selected = NULL, options = list(maxItems = 1))
    })
    observeEvent(input$rundge, {
      raw_counts$DEGS <- limma::topTable(raw_counts$model, coef = input$comparisons, number = Inf, adjust.method = input$adjust, sort.by = input$sortby, p.value = input$pcutoff)
      print(raw_counts$DEGS)
      output$DEGS <- DT::renderDataTable({
        DT::datatable(
          as.data.frame.matrix(raw_counts$DEGS),
          options = list(scrollX = TRUE, pageLength = 10)
        )
      })
      output$downloaddge <- renderUI({
        downloadButton("downloaddiffexp", "Download Table")
      })
    })
    output$downloaddiffexp <- downloadHandler(
      filename = function() {
        paste(input$comparisons, ".csv", sep = "")
      },
      content = function(con) {
        write.csv(raw_counts$DEGS, con)
      }
    )
    output$graphparams <- renderUI({
      req(input$graphtomake, raw_counts$dge, raw_counts$model)
      params <- graph_params(input$graphtomake, raw_counts$dge, raw_counts$model)
      params
    })
    graph <- eventReactive(input$graphbutton, {
      req(raw_counts$dge, input$graphtomake, raw_counts$model)
      if (input$graphtomake == "volcano") {
        top.table <- limma::topTable(raw_counts$model, coef = input$gcontrast, number = Inf, adjust.method = 'BH', p.value = 1)
        top.table <- BiocGenerics::as.data.frame(top.table)
        differential<-top.table
        mutateddf <- ggpubr::mutate(differential, sig = ifelse(differential$P.Value<input$gpvalue, "Sig", "Not Sig"))
        volc_in <- BiocGenerics::cbind(gene=mutateddf[,1], mutateddf)
        print((input$ggeneset$datapath))
        print(is.null(input$ggenes))
        if (is.null(input$ggeneset) & is.null(input$ggenes)){
          volc = ggplot2::ggplot(volc_in[2:length(volc_in$gene),], ggplot2::aes(logFC, -log10(P.Value))) +
            ggplot2::geom_point(ggplot2::aes(col=sig)) +
            ggplot2::scale_color_manual(values = c('Sig' = "red",  'Not Sig' = 'black')) +
            ggplot2::ggtitle(input$gcontrast)
          volc
        }
        else {
          if (!is.null(input$ggeneset)) {
          geneset <- readr::read_tsv(input$ggeneset$datapath)
          geneset <- geneset %>% ggpubr::mutate(SYMBOL = stringr::str_trim(SYMBOL))
          geneset <- as.list(geneset['SYMBOL'])
          sset<- volc_in[volc_in$genes %in% geneset$SYMBOL,]
        }
        else if (!is.null(input$ggenes)) {
          sset <- volc_in[volc_in$genes %in% input$ggenes,]
        }
        volc = ggplot2::ggplot(volc_in[2:length(volc_in$gene),], ggplot2::aes(logFC, -log10(P.Value))) +
            ggplot2::geom_point(ggplot2::aes(col=sig)) +
            ggplot2::scale_color_manual(values = c('Sig' = "red",  'Not Sig' = 'black', "Gene set" = "blue")) +
            ggplot2::geom_point(data = sset, ggplot2::aes(logFC, -log10(P.Value)), color = 'blue')+
            #geom_point(data = sset, aes(logFC, -log10(PValue)), color = 'green')+
            ggplot2::ggtitle(input$gcontrast)
        if (input$glabel) {
          volc + ggrepel::geom_text_repel(data = sset, ggplot2::aes(label=gene, segment.curvature = -0.2, segment.inflect = TRUE, segment.angle = 20, segment.ncp = 0, segment.linetype = 6),
          min.segment.length = ggplot2::unit(0,'lines'), nudge_x = 5, nudge_y = 4, arrow = ggplot2::arrow(length = ggplot2::unit(0.015,'npc')))
        }
        else {volc}
        }}
      else if (input$graphtomake == 'heatmap') {
        hm <- edgeR::cpm(raw_counts$dge, log = TRUE, prior.count = 1)
        groups <- c()
        for (i in seq_len(length(input$count))) {
          groups <- c(groups, input[[paste0("col", input$count[[i]])]])
        }
        colnames(hm) <- groups
        hm <- hm[,colnames(hm) %in% input$gcondition]
        rownames(hm) <- raw_counts$dge$genes$genes
        hm <- hm[rownames(hm) %in% input$ggenes,]
        pheatmap::pheatmap(hm)
      }
    })
    output$graph <- renderPlot({
      graph()
    })
    output$downloadgraph <- downloadHandler(
      filename = function() {
        if (input$graphtomake == 'volcano'){
          paste(input$graphtomake, input$gmode, sep = "")
        }
        else {
          paste(input$graphtomake, '.pdf', sep="")
        }
      },
      content = function(file) {
      if (input$graphtomake == 'heatmap') {
        hm <- edgeR::cpm(raw_counts$dge, log = TRUE, prior.count = 1)
        groups <- c()
        for (i in seq_len(length(input$count))) {
          groups <- c(groups, input[[paste0("col", input$count[[i]])]])
        }
        colnames(hm) <- groups
        hm <- hm[,colnames(hm) %in% input$gcondition]
        rownames(hm) <- raw_counts$dge$genes$genes
        hm <- hm[rownames(hm) %in% input$ggenes,]
        pdf(file, width=input$gwidth, height=input$gheight)
        pheatmap::pheatmap(hm)
        dev.off()

          #save_pheatmap(graph(), file, input$gwidth, input$gheight)
          #ggsave(file, plot = graph(), width = input$gwidth, height = input$gheight, units = input$gunits, scale = input$gscale)
      }
      else {
        ggplot2::ggsave(file, plot = graph(), width = input$gwidth, height = input$gheight, units = input$gunits, scale = input$gscale)
      }
      }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
