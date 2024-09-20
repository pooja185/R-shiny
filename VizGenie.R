# Load required libraries
library(shiny)
library(shinythemes)
library(ggplot2)
library(DESeq2)
library(plotly)  # For interactive plots

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  
  titlePanel("VizGenie"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("counts_file", "Upload Raw Counts CSV", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      fileInput("metadata_file", "Upload Metadata CSV", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      
      selectInput("design_column", "Select Design Column", choices = NULL),
      selectInput("ref_level", "Select Reference Level", choices = NULL),  # New dropdown for reference level
      
      actionButton("run_deseq", "Run DESeq2 Analysis")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Results", verbatimTextOutput("deseq_results")),
        tabPanel("Volcano Plot", plotOutput("volcano_plot", width = "100%", height = "400px"), downloadButton("download_volcano", "Download Volcano Plot")),
        tabPanel("PCA Plot", plotOutput("pca_plot", width = "100%", height = "400px"), downloadButton("download_pca", "Download PCA Plot"))
      )
    )
  )
)

# Define server
server <- function(input, output, session) {
  
  # Reactive values to store uploaded data
  data <- reactiveValues(counts = NULL, metadata = NULL)
  
  # Observe the counts upload
  observeEvent(input$counts_file, {
    req(input$counts_file)
    data$counts <- read.csv(input$counts_file$datapath, row.names = 1)
  })
  
  # Observe the metadata upload
  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    data$metadata <- read.csv(input$metadata_file$datapath, row.names = 1)
    
    # Update design column choices based on the metadata
    updateSelectInput(session, "design_column", choices = colnames(data$metadata))
  })
  
  # Update reference level choices based on selected design column
  observeEvent(input$design_column, {
    req(data$metadata)
    unique_levels <- unique(data$metadata[[input$design_column]])
    updateSelectInput(session, "ref_level", choices = unique_levels)
  })
  
  # Run DESeq2 analysis when the button is clicked
  results_deseq <- eventReactive(input$run_deseq, {
    req(data$counts, data$metadata, input$design_column, input$ref_level)
    
    # Check if column names match
    if (!all(colnames(data$counts) %in% rownames(data$metadata))) {
      stop("Column names in counts data do not match row names in metadata!")
    }
    
    if (!all(colnames(data$counts) == rownames(data$metadata))) {
      stop("Column names in counts data and row names in metadata are not in the same order!")
    }
    
    # Construct DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(countData = data$counts,
                                  colData = data$metadata,
                                  design = as.formula(paste("~", input$design_column)))
    
    # Prefiltering to remove low count genes
    filtered <- rowSums(counts(dds)) >= 10
    dds <- dds[filtered,]
    
    # Set factor levels based on user selection
    dds[[input$design_column]] <- relevel(dds[[input$design_column]], ref = input$ref_level)
    
    # Run DESeq
    dds <- DESeq(dds)
    
    # Return results
    return(list(dds = dds, results = results(dds)))
  })
  
  # Output DESeq2 results
  output$deseq_results <- renderPrint({
    req(results_deseq())
    summary(results_deseq()$results)
  })
  
  # Render Volcano Plot
  output$volcano_plot <- renderPlot({
    req(results_deseq())
    res <- as.data.frame(results_deseq()$results)
    
    plot(log2FoldChange ~ -log10(pvalue), data = res, pch = 20, main = "Volcano Plot", xlim = c(-3, 3))
    with(subset(res, padj < .01), points(log2FoldChange, -log10(pvalue), pch = 20, col = "blue"))
    with(subset(res, padj < .01 & abs(log2FoldChange) > 1), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))
  })
  
  # Download Volcano Plot
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste("volcano_plot-", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      png(file, width = 800, height = 600)  # Specify dimensions
      print(renderPlot({ output$volcano_plot() }))
      dev.off()
    }
  )
  
  # Render PCA Plot
  output$pca_plot <- renderPlot({
    req(results_deseq())  # Ensure DESeq analysis has been run
    
    dds <- results_deseq()$dds  # Access the dds object
    req(dds)  # Ensure dds is not NULL
    
    # Perform variance stabilizing transformation
    vsdata <- vst(dds, blind = FALSE)
    
    # Ensure input$design_column is available
    req(input$design_column)
    
    # Create PCA plot
    plotPCA(vsdata, intgroup = input$design_column)
  })
  
  # Download PCA Plot
  output$download_pca <- downloadHandler(
    filename = function() {
      paste("pca_plot-", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      png(file, width = 800, height = 600)  # Specify dimensions
      print(renderPlot({ output$pca_plot() }))
      dev.off()
    }
  )
  
}

# Run the app
shinyApp(ui = ui, server = server)
