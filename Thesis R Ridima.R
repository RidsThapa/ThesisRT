# Libraries
library(shiny)
library(plotly)
library(DT)
library(dplyr)

# User Interface
ui <- fluidPage(
  titlePanel("Gene essentiality data visualisation at domain level"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("tnseq", "Upload Tn-seq Data (.csv)"),
      fileInput("domains", "Upload Pfam/InterPro Annotations (.csv)"),
      downloadButton("report", "Results"),
      downloadButton("bedExport", "Export BED for JBrowse2")
    ),
    
    mainPanel(
      h3("Automated Visualisation"),
      plotlyOutput("genomePlot", height = "400px"),
      DTOutput("essentialTable"),
      
      tags$hr(),
      
      h3("Genome Browser (JBrowse2)"),
      tags$iframe(
        src = "http://localhost:3000",  
        style = "width:100%;height:800px;border:none;"
      )
    )
  )
)

server <- function(input, output) {
  
  # Load tnseq data
  tnseq_data <- reactive({
    req(input$tnseq)
    read.csv(input$tnseq$datapath)
  })
  
  # Load domain data
  domain_data <- reactive({
    req(input$domains)
    read.csv(input$domains$datapath)
  })
  
  # Genome plot
  output$genomePlot <- renderPlotly({
    req(tnseq_data(), domain_data())
    
    df <- tnseq_data()
    doms <- domain_data()
    
    # Set colors for essentiality
    essentiality_colors <- c(
      "essential" = "red",
      "non-essential" = "blue",
      "unknown" = "grey"
    )
    
    p <- plot_ly(df, x = ~position, y = ~insertions, type = "bar",
                 color = ~essentiality,
                 colors = essentiality_colors,
                 text = ~paste("Gene:", gene_id, "<br>Essentiality:", essentiality)) %>%
      layout(title = "Insertion Profile with Essentiality Calls",
             xaxis = list(title = "Genomic Position"),
             yaxis = list(title = "Insertions"))
    
    # Add domains
    if (nrow(doms) > 0) {
      for (i in 1:nrow(doms)) {
        showlegend <- ifelse(i == 1, TRUE, FALSE) 
        p <- add_trace(p,
                       x = c(doms$domain_start[i], doms$domain_end[i]),
                       y = c(max(df$insertions), max(df$insertions)),
                       type = "scatter", mode = "lines",
                       line = list(color = 'red', width = 5),
                       name = "Domain",
                       showlegend = showlegend,
                       inherit = FALSE)
      }
    }
    
    p
  })
  
  # Essentiality + domain table
  output$essentialTable <- renderDT({
    req(tnseq_data(), domain_data())
    
    # Summary of insertions per gene
    tnseq_summary <- tnseq_data() %>%
      group_by(gene_id, essentiality) %>%
      summarise(total_insertions = sum(insertions, na.rm = TRUE), .groups = "drop")
    
    # Add domains
    inner_join(tnseq_summary, domain_data(), by = "gene_id") %>%
      dplyr::select(
        gene_id,
        total_insertions,
        essentiality,
        domain_name,
        domain_start,
        domain_end
      )
  },
  options = list(pageLength = 25, autoWidth = TRUE, dom = "Blfrtip"),
  filter = "top")
  
  # CSV summary export
  output$report <- downloadHandler(
    filename = function() { "essentiality_report.csv" },
    content = function(file) {
      tnseq_summary <- tnseq_data() %>%
        group_by(gene_id, essentiality) %>%
        summarise(total_insertions = sum(insertions, na.rm = TRUE), .groups = "drop")
      
      df <- inner_join(tnseq_summary, domain_data(), by = "gene_id") %>%
        dplyr::select(
          gene_id,
          total_insertions,
          essentiality,
          domain_name,
          domain_start,
          domain_end
        )
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # BED file export
  output$bedExport <- downloadHandler(
    filename = function() { "essential_domains.bed" },
    content = function(file) {
      tnseq_summary <- tnseq_data() %>%
        group_by(gene_id, essentiality) %>%
        summarise(total_insertions = sum(insertions, na.rm = TRUE), .groups = "drop")
      
      df <- inner_join(tnseq_summary, domain_data(), by = "gene_id")
      bed <- df %>%
        transmute(chr = "chromosome1", start = domain_start, end = domain_end, name = gene_id)
      write.table(bed, file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  )
}

shinyApp(ui, server)
