#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(tidyverse)
library(readxl)
library(edgeR)
library(limma)

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
        textInput('genes', 'List column number containing gene names (e.g 1)'),
        textInput('count', 'List column number of samples you would like to analyze separated by commas (e.g 2,3)'),
        textInput('group', 'List the label for each sample selected separated by commas (e.g KO_Naive,WT_Naive)'),
        actionButton('selection', 'Done')
      ),
      mainPanel(
        plotOutput('bplot'),
        plotOutput('mds')
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
    observeEvent(input$loaddata, {
        raw_counts$obj <- utils::read.csv(input$rawcounts$datapath)
        output$counts <- DT::renderDataTable({
          DT::datatable(
            as.data.frame.matrix(raw_counts$obj),
            options = list(scrollX = TRUE, pageLength = 10)
          )
        })

    })
    observeEvent(input$selection, {
      selected <- as.numeric(strsplit(input$count, ',')[[1]])
      genes <- as.numeric(input$genes)
      groups <- strsplit(str_replace_all(input$group, ' ', ''), ',')[[1]]
      raw_counts$dge <- DGEList(counts = raw_counts$obj[,selected], genes = raw_counts$obj[,genes],
                                group = groups) %>% calcNormFactors(method = "TMM")
      raw_counts$design.mat <- model.matrix(~ 0 + raw_counts$dge$samples$group)
      colnames(raw_counts$design.mat) <- levels(raw_counts$dge$samples$group)
      keep <- filterByExpr(raw_counts$dge, design = raw_counts$design.mat)
      raw_counts$dge <- raw_counts$dge[keep,]
      output$bplot <- renderPlot({
        barplot(raw_counts$dge$samples$lib.size, names=colnames(raw_counts$dge),las=2)
      })
      output$mds <- renderPlot({
        colConditions <- rainbow(length(levels(raw_counts$dge$samples$group)))
        colConditions <- colConditions[match((raw_counts$dge$samples$group), levels(raw_counts$dge$samples$group))]
        plotMDS(raw_counts$dge, col = colConditions, labels = raw_counts$dge$samples$group)
      })

    })
    observeEvent(input$fit, {
      if(input$model == 'Limma'){
          raw_counts$model <- fitlimma(raw_counts$dge, raw_counts$design.mat, groups)
      }
      else {
        output$vplot <- renderPlot({
          raw_counts$model <- fitvoom(raw_counts$dge, raw_counts$design.mat)
        })
      }
      output$SA <- renderPlot({
        plotSA(raw_counts$model)
      })
      output$MD <- renderPlot({
        plotMD(raw_counts$model)
      })
    })
    output$selectcon <- renderUI({
      selectizeInput("selected", "Select Contrast", choices = levels(raw_counts$dge$samples$group), selected = NULL, options = list(maxItems = 2))
    })
    observeEvent(input$add, {
      raw_counts$contrasts <- append(raw_counts$contrasts, str_c(input$selected[1],'-',input$selected[2]))
    })
    observeEvent(input$clear, {
      raw_counts$contrasts <- c()
      #output$currcons <- renderText({
      #  paste('Current selections:')
      #})
    })
    output$currcons <- renderText({
      paste(append(raw_counts$contrasts, 'Current selections:', 0), collapse = ", ")
    })
    observeEvent(input$done, {
      if (is_empty(raw_counts$contrasts)) return()
      raw_counts$contrastmat <- makeContrasts(contrasts = raw_counts$contrasts, levels = raw_counts$design.mat)
      raw_counts$model <- contrasts.fit(raw_counts$model, raw_counts$contrastmat)
      raw_counts$model <- eBayes(raw_counts$model)
    })
    output$comparison <- renderUI({
      selectizeInput("comparisons", "Select Comparison", choices = levels(raw_counts$dge$samples$group), selected = NULL, options = list(maxItems = 2))
    })
    observeEvent(input$rundge, {
      raw_counts$DEGS <- topTable(raw_counts$model, coef = str_c(input$comparisons[1], '-', input$comparisons[2]), number = Inf, adjust.method = input$adjust, sort.by = input$sortby, p.value = input$pcutoff)
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
        paste(input$comparisons[1], "_vs_", input$comparisons[2], ".csv", sep = "")
      },
      content = function(con) {
        write.csv(raw_counts$DEGS, con)
      }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
