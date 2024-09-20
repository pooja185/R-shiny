# Load required libraries
library(shiny)
library(shinythemes)
library(ggplot2)

# Load Dataset
metadata <- read.csv("Deidentified_example_RNAseq_metadata.csv", row.names = "sample")
norm.counts <- read.csv("Deidentified_example_RNAseq_normalized_log2_counts.csv", row.names = "gene")

# Check if the row names of metadata match the column names of gene.df
if (!all(rownames(metadata) == colnames(norm.counts))) {
  stop("The metadata rows and normalized count columns do not match!")
}

# Transpose gene data for easier handling
gene.df <- t(norm.counts)

# Define UI
ui <- fluidPage(
  theme = shinytheme("united"),
  
  navbarPage(
    title = "VizWise",
    header = tags$style(".navbar-brand { font-size: 25px; }"),
    
    sidebarLayout(
      sidebarPanel(
        br(),
        selectInput("x_dataset", "Select dataset for X-axis", choices = c("metadata", "gene.df")),
        selectizeInput("x_variable", "Select the desired X-variable", choices = NULL),
        
        selectInput("y_dataset", "Select dataset for Y-axis", choices = c("metadata", "gene.df")),
        selectizeInput("y_variable", "Select the desired Y-variable", choices = NULL),
        
        selectInput("group_variable", "Select Grouping Variable", choices = NULL),
        
        actionButton("generate_summary", "Generate Summary"),
        actionButton("generate_plot", "Generate Plot")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel("Boxplot", plotOutput("boxplot")),
          tabPanel("Violin Plot", plotOutput("violin_plot")),
          tabPanel("Bar Graph", plotOutput("bar_graph")),
          tabPanel("Scatter Plot", plotOutput("scatter_plot")),
          tabPanel("Data Summary", verbatimTextOutput("summary_output"))
        )
      )
    )
  )
)

# Define server
server <- function(input, output, session) {
  
  # Update x variable choices based on selected dataset with server-side selectize
  observeEvent(input$x_dataset, {
    choices <- if (input$x_dataset == "metadata") colnames(metadata) else colnames(gene.df)
    updateSelectizeInput(session, "x_variable", choices = choices, server = TRUE)
  })
  
  # Update y variable choices based on selected dataset with server-side selectize
  observeEvent(input$y_dataset, {
    choices <- if (input$y_dataset == "metadata") colnames(metadata) else colnames(gene.df)
    updateSelectizeInput(session, "y_variable", choices = choices, server = TRUE)
  })
  
  # Update group variable choices based on metadata only
  observe({
    updateSelectInput(session, "group_variable", choices = colnames(metadata))
  })
  
  # Generate summary based on the selected x and y variables
  summary_data <- reactive({
    req(input$x_variable, input$y_variable, input$group_variable)
    
    x_data <- if (input$x_dataset == "metadata") metadata[, input$x_variable, drop = TRUE] else gene.df[, input$x_variable, drop = TRUE]
    y_data <- if (input$y_dataset == "metadata") metadata[, input$y_variable, drop = TRUE] else gene.df[, input$y_variable, drop = TRUE]
    
    group_data <- if (input$x_dataset == "metadata") metadata[, input$group_variable, drop = TRUE] else rep(NA, nrow(gene.df))
    
    data <- data.frame(X = x_data, Y = y_data, Group = factor(group_data))
    
    return(data)
  })
  
  # Render summary output
  output$summary_output <- renderPrint({
    req(input$generate_summary)
    summary(summary_data())
  })
  
  # Render Boxplot
  output$boxplot <- renderPlot({
    req(input$generate_plot)
    data <- summary_data()
    
    if (nrow(data) == 0) {
      return(NULL)
    }
    
    ggplot(data, aes(x = Group, y = Y, fill = Group)) +
      geom_boxplot() +
      geom_point(position = position_jitter(width = 0.2), size = 1.5) +
      labs(title = "Boxplot", x = input$group_variable, y = input$y_variable) +
      scale_fill_brewer(palette = "Set1")
  })
  
  # Render Violin Plot
  output$violin_plot <- renderPlot({
    req(input$generate_plot)
    data <- summary_data()
    
    if (nrow(data) == 0) {
      return(NULL)
    }
    
    ggplot(data, aes(x = Group, y = Y, fill = Group)) +
      geom_violin(trim = FALSE) +
      geom_point(position = position_jitter(width = 0.2), size = 1.5) +
      labs(title = "Violin Plot", x = input$group_variable, y = input$y_variable) +
      scale_fill_brewer(palette = "Set1")
  })
  
  # Render Bar Graph
  output$bar_graph <- renderPlot({
    req(input$generate_plot)
    data <- summary_data()
    
    if (nrow(data) == 0) {
      return(NULL)
    }
    
    # Summarize Y by group
    data_summary <- aggregate(Y ~ Group, data = data, FUN = mean)
    
    ggplot(data_summary, aes(x = Group, y = Y, fill = Group)) +
      geom_bar(stat = "identity") +
      labs(title = "Bar Graph", x = input$group_variable, y = input$y_variable) +
      scale_fill_brewer(palette = "Set1")
  })
  
  # Render Scatter Plot
  output$scatter_plot <- renderPlot({
    req(input$generate_plot)
    data <- summary_data()
    
    if (nrow(data) == 0) {
      return(NULL)
    }
    
    ggplot(data, aes(x = X, y = Y, color = Group, fill = Group)) +
      geom_point(size = 2) +
      labs(title = "Scatter Plot", x = input$x_variable, y = input$y_variable) +
      scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1")
  })
}

# Run the app
shinyApp(ui = ui, server = server)
