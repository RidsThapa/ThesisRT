# Libraries
library(shiny)
library(plotly)
library(DT)
library(dplyr)


# User Interface
ui <- fluidPage(
  
  titlePanel("Gene Essentiality Data Visualisation at Domain Level"),
  
  sidebarLayout(sidebarPanel(
      
      fileInput("tnseq", "Upload Tn-seq Data (.csv)"),
      
      fileInput("domains", "Upload Pfam/InterPro Annotations (.csv)"),
      
      selectInput("selected_gene", "Select Gene:", choices = NULL),
      
      selectInput("essential_filter", "Essentiality Filter:",
        choices = c("All","Potentially essential","Low insertion","No insertion")),
      
      downloadButton("report", "Download Results")),
    
    
    mainPanel(tabsetPanel(tabPanel("Domain Visualisation", 
          br(), h3("Selected Gene"),
          verbatimTextOutput("geneSummary"),
          plotlyOutput("genePlot", height = "600px"), tags$hr(),
          h3("Domain-level Results"),
          DTOutput("essentialTable")),
        
        tabPanel("Data Summary", br(), h3("Dataset Summary"),
          verbatimTextOutput("Summary")
        )
      )
    )
  )
)



# Server
server <- function(input, output, session) {
  
  # Load Tn-seq data
  tnseq_data <- reactive({req(input$tnseq)
    tnseq <- read.csv(input$tnseq$datapath, stringsAsFactors = FALSE)
    tnseq <- tnseq %>%
      mutate(gene_id = tolower(gene_id), position = as.numeric(position), insertions = as.numeric(insertions)
      ) %>%
      filter(!is.na(gene_id), !is.na(position))
    tnseq
  })
  
  # Load domain data
  domain_data <- reactive({req(input$domains)
    domains <- read.csv(input$domains$datapath, stringsAsFactors = FALSE)
    
  # Use Pfam or InterPro if no domain_name
    if (!"domain_name" %in% names(domains)) {
      
      if ("Pfam" %in% names(domains)) {domains$domain_name <- domains$Pfam} 
      if ("InterPro" %in% names(domains)) {domains$domain_name <- domains$InterPro} 
      else {domains$domain_name <- "Domain"}
    }
    
    domains <- domains %>%
      mutate(gene_id = tolower(gene_id),
      domain_start = as.numeric(domain_start),
      domain_end = as.numeric(domain_end)
      ) %>%
      filter(!is.na(gene_id), !is.na(domain_start), !is.na(domain_end))
    domains
  })
  
  # Find common genes in both datasets
  common_genes <- reactive({
    intersect(unique(tnseq_data()$gene_id), unique(domain_data()$gene_id))
  })
  
  # Update gene selection
  observe({genes <- sort(common_genes())
    
    updateSelectInput(session, "selected_gene", choices = genes)
  })
  
  # Match insertions with domains
  domain_insertions <- reactive({
    tnseq <- tnseq_data()
    domains <- domain_data()
    data <- inner_join(domains, tnseq, by = "gene_id")
    data <- data %>%
      filter(position >= domain_start, position <= domain_end)
    data
  })
  
  # Calculate domain-level essentiality
  domain_summary <- reactive({
    domains <- domain_data()
    insertions <- domain_insertions()
    
  # Count insertions in each domain
    counts <- insertions %>%
      group_by(gene_id, domain_name, domain_start, domain_end) %>%
      summarise(
        total_insertions = sum(insertions, na.rm = TRUE),
        insertion_sites = sum(insertions > 0, na.rm = TRUE), .groups = "drop")
    
  # Add insertion counts to domains
    results <- domains %>%
      distinct(gene_id, domain_name, domain_start, domain_end, .keep_all = TRUE) %>%
      left_join(counts, by = c("gene_id", "domain_name", "domain_start", "domain_end"))

  # Classify domain insertion results
    results <- results %>%
      mutate(domain_length = domain_end - domain_start + 1, 
        
        domain_status = case_when(
          is.na(total_insertions) ~ "No TnSeq data",
          total_insertions == 0 ~ "Potentially essential",
          total_insertions <= 2 ~ "Low insertion", TRUE ~ "No insertion"),
        total_insertions = ifelse(is.na(total_insertions), 0, total_insertions),
        insertion_sites = ifelse(is.na(insertion_sites), 0, insertion_sites)
      )  
    results
  })
      
  # Filter domain results
  filtered_domains <- reactive({results <- domain_summary()
  if (input$essential_filter != "All") 
    {results <- results %>%
    filter(domain_status == input$essential_filter)}
  results})
  
  # Selected gene Tn-seq data
  selected_tnseq <- reactive({ req(input$selected_gene)
    tnseq_data() %>%
      filter(gene_id == input$selected_gene)
  })
  
  
  
  # Selected gene domain data
  selected_domains <- reactive({ req(input$selected_gene)
    filtered_domains() %>%
      filter(gene_id == input$selected_gene)
  })
  
  # Gene summary
  output$geneSummary <- renderText({ req(input$selected_gene)
    
    domains <- selected_domains()
    paste0("Gene: ", input$selected_gene,
      "\nDomains: ", nrow(domains),
      "\nTotal domain insertions: ", sum(domains$total_insertions),
      "\nPotentially essential domains: ", sum(domains$domain_status == "Potentially essential"),
      "\nDomains with no insertions: ", sum(domains$total_insertions == 0)
    )
  })
  
  
  
  # Domain-level plot
  output$genePlot <- renderPlotly({ req(input$selected_gene)
    
    tnseq <- selected_tnseq()
    domains <- selected_domains()
    
  # Check if domains are available
    validate(need(nrow(domains) > 0, "No domains found.")
    )
    
    
  # Keep insertion positions within filtered domains
    tnseq <- tnseq %>%
      filter(position %in% unlist(lapply(1:nrow(domains),
            function(i) {seq(domains$domain_start[i],domains$domain_end[i])}
          )
        )
      )
    
    
  # Tn-seq insertion plot
    insertion_plot <- plot_ly(tnseq, x = ~position, y = ~insertions, type = "bar", name = "Insertions") %>%
      layout(yaxis = list(title = "Insertions"))
    
  # Protein domain plot
    domain_plot <- plot_ly()
    
    for (i in 1:nrow(domains)) {
      
      domain_plot <- domain_plot %>%
        add_segments(
          
          x = domains$domain_start[i], xend = domains$domain_end[i],
          
          y = domains$domain_name[i], yend = domains$domain_name[i],
          
          line = list(width = 15),
          
          name = domains$domain_name[i],
          
          text = paste0(
            "Gene: ", domains$gene_id[i],
            "<br>Domain: ", domains$domain_name[i],
            "<br>Start: ", domains$domain_start[i],
            "<br>End: ", domains$domain_end[i],
            "<br>Insertions: ", domains$total_insertions[i],
            "<br>Status: ", domains$domain_status[i]),
          
          hoverinfo = "text", showlegend = FALSE
        )
    }
    
    
    # Combine plots
    subplot(insertion_plot, domain_plot, nrows = 2, shareX = TRUE, heights = c(0.6, 0.4)) %>%
      layout(title = paste("Domain-level essentiality:", input$selected_gene),
        
        xaxis2 = list(title = "Genomic Position"),
        yaxis = list(title = "Insertions"),
        yaxis2 = list(title = "Protein Domain")
      )
  })
  
  
  
  # Domain results table
  output$essentialTable <- renderDT({
    
    results <- filtered_domains() %>%
      select(gene_id,domain_name, domain_start, domain_end, domain_length,
        total_insertions, insertion_sites, domain_status)
    
    datatable(results, filter = "top", rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE))
  })
  
  
  
  # Dataset summary
  output$Summary <- renderText({
    tnseq <- tnseq_data()
    domains <- domain_data()
    common <- common_genes()
    
    paste0("Tn-seq rows: ", nrow(tnseq),
      
      "\nUnique Tn-seq genes: ", length(unique(tnseq$gene_id)),
      "\n\nDomain rows: ", nrow(domains),
      "\nUnique domain genes: ", length(unique(domains$gene_id)),
      "\n\nGenes in both datasets: ", length(common),
      "\n\nExample overlapping genes:\n",
      paste(head(common, 20), collapse = ", ")
    )
  })
  
  
  
  # Download results
  output$report <- downloadHandler(
    filename = function() {paste0("domain_essentiality_results_", Sys.Date(), ".csv")},
    content = function(file) {write.csv(filtered_domains(), file, row.names = FALSE)}
  )
}


shinyApp(ui, server)
