library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(DT)
library(shinyFiles)
library(gridExtra)
library(cowplot)
library(scales)
library(RColorBrewer)

# Modern UI with shinydashboard
ui <- dashboardPage(
  # Header
  dashboardHeader(
    title = div(
      "Population Genomics Viewer"
    ),
    titleWidth = 250
  ),
  
  # Sidebar
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      # Directory Selection
      div(
        style = "padding: 10px;",
        shinyDirButton("folder", 
                       "Choose Data Directory", 
                       "Please select a directory",
                       class = "btn-primary btn-block")
      ),
      verbatimTextOutput("selected_dir"),
      
      # Metadata Upload
      div(
        style = "padding: 10px;",
        fileInput("metadata", 
                  "Upload Metadata File",
                  accept = c("text/csv", "text/comma-separated-values", ".csv"),
                  width = "100%")
      ),
      
      # MAG Selection
      div(
        style = "padding: 10px;",
        selectInput("mag", "Select MAG:", choices = NULL, width = "100%")
      ),
      
      # Plot Type Selection
      div(
        style = "padding: 10px;",
        selectInput("plot_type", "Select Plot Type:", 
                    choices = c(
                      "FST Heatmap" = "fst",
                      "Correlation Matrix" = "correlation",
                      "Intradiversity Analysis" = "intradiv"
                    ), 
                    width = "100%")
      ),
      
      # Conditional Panels
      conditionalPanel(
        condition = "input.plot_type == 'correlation'",
        div(
          style = "padding: 10px;",
          checkboxGroupInput("corr_vars", 
                             "Select Variables for Correlation:", 
                             choices = NULL,
                             width = "100%")
        )
      ),
      
      conditionalPanel(
        condition = "input.plot_type == 'intradiv'",
        div(
          style = "padding: 10px;",
          selectInput("color_var", "Color by:", choices = NULL, width = "100%"),
          checkboxGroupInput("plot_vars", "Select Variables to Plot:", choices = NULL, width = "100%"),
          sliderInput("point_size", "Point Size:", 
                      min = 1, max = 10, value = 3, step = 0.5,
                      width = "100%"),
          selectInput("color_palette", "Color Palette:", 
                      choices = c(
                        "Husl" = "husl",
                        "Set1" = "Set1",
                        "Set2" = "Set2",
                        "Set3" = "Set3"
                      ), 
                      width = "100%")
        )
      ),
      
      # Download Button
      div(
        style = "padding: 10px;",
        downloadButton("downloadPlot", "Download Plot", class = "btn-success btn-block")
      )
    )
  ),
  
  # Body
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .box { 
          box-shadow: 0 4px 6px rgba(0,0,0,0.1); 
          border-top: 3px solid #2c3e50;
        }
        .main-sidebar { 
          background-color: #34495e !important; 
        }
        .sidebar-menu > li.active > a { 
          background-color: #2c3e50 !important; 
        }
        .nav-tabs-custom > .nav-tabs {
          background-color: #ecf0f5;
        }
        body { 
          background-color: #ecf0f5; 
        }
        .dataTables_wrapper {
          overflow-x: auto;
          max-width: 100%;
        }
      "))
    ),
    fluidRow(
      box(
        title = "Plot", 
        width = 12, 
        status = "primary", 
        solidHeader = TRUE,
        plotOutput("plot", height = "600px")
      ),
      box(
        title = "Data Table", 
        width = 12, 
        status = "primary", 
        solidHeader = TRUE,
        div(
          style = "overflow-x: auto; width: 100%;",
          DTOutput("table")
        )
      )
    )
  )
)

# Server remains the same as in the previous version
server <- function(input, output, session) {
  
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  shinyDirChoose(input, "folder", roots = volumes, session = session)
  
  data_dir <- reactive({
    getwd()
  })
  
  metadata <- reactive({
    df <- read.csv("metadata.csv", stringsAsFactors = FALSE)
    names(df) <- make.names(names(df), unique = TRUE)
    return(df)
  })
  
  # Update variable choices when metadata is uploaded
  observe({
    req(metadata())
    
    # Get numeric columns (excluding Sample_ID and Pond_Name)
    numeric_cols <- names(metadata())[sapply(metadata(), is.numeric)]
    
    # Update correlation variable choices
    updateCheckboxGroupInput(session, "corr_vars",
                             choices = numeric_cols,
                             selected = numeric_cols[1:min(length(numeric_cols), 5)])
    
    # Update intradiversity plot choices
    updateCheckboxGroupInput(session, "plot_vars",
                             choices = numeric_cols,
                             selected = numeric_cols[1:min(length(numeric_cols), 3)])
    
    # Update color variable choices (all columns)
    updateSelectInput(session, "color_var",
                      choices = names(metadata()))
  })
  
  output$selected_dir <- renderText({
    req(data_dir())
    paste("Selected directory:", data_dir())
  })
  
  observe({
    req(data_dir())
    mag_files <- list.files(data_dir(), pattern = "*.fst.txt$")
    mag_names <- unique(gsub("\\.fst\\.txt$", "", mag_files))
    updateSelectInput(session, "mag", choices = mag_names)
  })
  
  read_fst_data <- function(mag_name) {
    file_path <- file.path(data_dir(), paste0(mag_name, ".fst.txt"))
    df <- read.table(file_path, header = TRUE, row.names = 1)
    df[] <- lapply(df, function(x) as.numeric(as.character(x)))
    mat <- as.matrix(df)
    mat[upper.tri(mat, diag = FALSE)] <- NA
    return(mat)
  }
  
  read_intradiv_data <- function(mag_name) {
    file_path <- file.path(data_dir(), paste0(mag_name, ".intradiv.txt"))
    df <- read.table(file_path, header = TRUE)
    return(df)
  }
  
  output$plot <- renderPlot({
    if(input$plot_type == "correlation") {
      req(metadata(), input$corr_vars)
      
      # Select only the chosen numeric columns
      numeric_data <- metadata()[, input$corr_vars, drop = FALSE]
      
      corr_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
      mask <- upper.tri(corr_matrix, diag = FALSE)
      corr_matrix[mask] <- NA
      
      colors <- colorRampPalette(c("#4575B4", "white", "#D73027"))(100)
      
      pheatmap(corr_matrix,
               color = colors,
               display_numbers = TRUE,
               number_format = "%.2f",
               fontsize_number = 10,
               main = "Correlation Matrix",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               na_col = "white",
               breaks = seq(-1, 1, length.out = 101))
      
    } else if(input$plot_type == "intradiv") {
      req(input$mag, metadata(), input$plot_vars, input$color_var)
      
      # Read intradiversity data
      intradiv_data <- read_intradiv_data(input$mag)
      
      # Merge using Pond_Name from metadata and Sample from intradiv
      combined_data <- merge(
        intradiv_data,
        metadata(),
        by.x = "Sample",
        by.y = "Pond_Name",
        all.x = TRUE
      )
      
      # Check if merge was successful
      if(nrow(combined_data) == 0) {
        stop("No matching samples found between metadata and intradiversity data")
      }
      
      # Create list to store plots
      plots <- list()
      
      # Create color palette
      n_colors <- length(unique(combined_data[[input$color_var]]))
      if(input$color_palette == "husl") {
        colors <- scales::hue_pal()(n_colors)
      } else {
        colors <- RColorBrewer::brewer.pal(min(n_colors, 9), input$color_palette)
      }
      
      # Create individual plots
      for(var in input$plot_vars) {
        p <- ggplot(combined_data, aes_string(x = var, y = "Norm_intra_pi", color = input$color_var)) +
          geom_point(size = input$point_size) +
          theme_minimal() +
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                legend.position = "right") +
          scale_color_manual(values = colors) +
          labs(x = var,
               y = if(var == input$plot_vars[1]) "Normalized intra pi" else "")
        
        plots[[var]] <- p
      }
      
      # Arrange plots in grid
      do.call(grid.arrange, c(plots, 
                              list(ncol = min(3, length(plots)), 
                                   top = paste("Intradiversity Analysis -", input$mag))))
      
    } else if(input$plot_type == "fst") {
      req(input$mag)
      df <- read_fst_data(input$mag)
      
      colors <- colorRampPalette(c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", 
                                   "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026"))(100)
      
      pheatmap(df,
               color = colors,
               display_numbers = TRUE,
               number_format = "%.3f",
               main = paste("FST Heatmap -", input$mag),
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               na_col = "white",
               breaks = seq(0, 0.5, length.out = 101),
               show_diag = TRUE,
               na_show = FALSE)
    }
  })
  
  # Enhanced table rendering to ensure responsiveness
  output$table <- renderDT({
    req(input$plot_type)
    
    # Existing data retrieval logic
    data <- if(input$plot_type == "correlation") {
      req(metadata())
      metadata()
    } else if(input$plot_type == "fst") {
      req(input$mag)
      df <- read_fst_data(input$mag)
      as.data.frame(df)
    } else if(input$plot_type == "intradiv") {
      req(input$mag)
      read_intradiv_data(input$mag)
    }
    
    # Enhanced datatable rendering
    datatable(
      data,
      extensions = 'Responsive',
      options = list(
        scrollX = TRUE,  # Horizontal scrolling
        scrollCollapse = TRUE,
        pageLength = 10,
        lengthMenu = c(10, 25, 50, 100)
      )
    )
  })

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0(input$mag, "_", input$plot_type, ".pdf")
    },
    content = function(file) {
      pdf(file, width = 12, height = 10)
      if(input$plot_type == "correlation") {
        req(metadata(), input$corr_vars)
        numeric_data <- metadata()[, input$corr_vars, drop = FALSE]
        corr_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
        mask <- upper.tri(corr_matrix, diag = FALSE)
        corr_matrix[mask] <- NA
        colors <- colorRampPalette(c("#4575B4", "white", "#D73027"))(100)
        
        pheatmap(corr_matrix,
                 color = colors,
                 display_numbers = TRUE,
                 number_format = "%.2f",
                 fontsize_number = 10,
                 main = "Correlation Matrix",
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 na_col = "white",
                 breaks = seq(-1, 1, length.out = 101))
        
      } else if(input$plot_type == "intradiv") {
        req(input$mag, metadata(), input$plot_vars, input$color_var)
        intradiv_data <- read_intradiv_data(input$mag)
        combined_data <- merge(
          intradiv_data,
          metadata(),
          by.x = "Sample",
          by.y = "Pond_Name",
          all.x = TRUE
        )
        
        n_colors <- length(unique(combined_data[[input$color_var]]))
        if(input$color_palette == "husl") {
          colors <- scales::hue_pal()(n_colors)
        } else {
          colors <- RColorBrewer::brewer.pal(min(n_colors, 9), input$color_palette)
        }
        
        plots <- list()
        for(var in input$plot_vars) {
          p <- ggplot(combined_data, aes_string(x = var, y = "Norm_intra_pi", color = input$color_var)) +
            geom_point(size = input$point_size) +
            theme_minimal() +
            theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14),
                  legend.position = "right") +
            scale_color_manual(values = colors) +
            labs(x = var,
                 y = if(var == input$plot_vars[1]) "Normalized intra pi" else "")
          
          plots[[var]] <- p
        }
        
        do.call(grid.arrange, c(plots, 
                                list(ncol = min(3, length(plots)), 
                                     top = paste("Intradiversity Analysis -", input$mag))))
        
      } else if(input$plot_type == "fst") {
        df <- read_fst_data(input$mag)
        colors <- colorRampPalette(c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", 
                                     "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026"))(100)
        
        pheatmap(df,
                 color = colors,
                 display_numbers = TRUE,
                 number_format = "%.3f",
                 main = paste("FST Heatmap -", input$mag),
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 na_col = "white",
                 breaks = seq(0, 0.5, length.out = 101),
                 show_diag = TRUE,
                 na_show = FALSE)
      }
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
