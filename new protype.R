library(shiny)
library(TwoSampleMR)
library(dplyr)
library(stringr)
library(tibble)
library(DT)

# Attempt to load available outcomes from TwoSampleMR once when the app starts
ao_raw <- tryCatch(
  TwoSampleMR::available_outcomes(),
  error = function(e) NULL
)
if (is.null(ao_raw) || !inherits(ao_raw, "data.frame")) {
  ao <- tibble()  # Empty tibble if fetching fails
  warning("Failed to load available outcomes. 'ao' is empty.")
} else {
  ao <- ao_raw %>% as_tibble()
}

# UI definition
ui <- fluidPage(
  titlePanel("MR Explorer: Select Exposure and Outcome Studies"),
  sidebarLayout(
    sidebarPanel(
      textInput("exposure_search", "Search Exposures:", value = ""),
      DTOutput("exposure_table"),
      verbatimTextOutput("selected_exposure"),
      hr(),
      textInput("outcome_search", "Search Outcomes:", value = ""),
      DTOutput("outcome_table"),
      verbatimTextOutput("selected_outcome"),
      hr(),
      actionButton("run_mr", "Run MR Analysis", class = "btn-primary"),
      width = 4
    ),
    mainPanel(
      h4("MR Results"),
      DTOutput("mr_results_table"),
      plotOutput("mr_scatter_plot"),
      width = 8
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Reactive subset for exposures based on search string
  exposures <- reactive({
    # If 'ao' is empty, return empty tibble with same columns to avoid errors
    if (nrow(ao) == 0) {
      return(ao)
    }
    if (input$exposure_search == "") {
      return(ao[0, ])  # return empty rows but preserve columns
    }
    ao %>% filter(str_detect(trait, regex(input$exposure_search, ignore_case = TRUE)))
  })
  
  # Render exposures table with selection enabled
  output$exposure_table <- renderDT({
    tbl <- exposures()
    if (nrow(tbl) == 0) {
      datatable(tbl, options = list(dom = 't'), rownames = FALSE)
    } else {
      datatable(
        tbl,
        selection = 'single',
        options = list(pageLength = 5, scrollY = '200px')
      )
    }
  })
  
  # Display selected exposure details (ID and trait)
  output$selected_exposure <- renderPrint({
    sel <- input$exposure_table_rows_selected
    tbl <- exposures()
    if (length(sel) && nrow(tbl) >= sel) {
      row <- tbl[sel, ]
      paste0("Selected Exposure ID: ", row$id, " (", row$trait, ")")
    } else {
      "No exposure selected"
    }
  })
  
  # Reactive subset for outcomes
  outcomes <- reactive({
    if (nrow(ao) == 0) {
      return(ao)
    }
    if (input$outcome_search == "") {
      return(ao[0, ])
    }
    ao %>% filter(str_detect(trait, regex(input$outcome_search, ignore_case = TRUE)))
  })
  
  # Render outcomes table
  output$outcome_table <- renderDT({
    tbl <- outcomes()
    if (nrow(tbl) == 0) {
      datatable(tbl, options = list(dom = 't'), rownames = FALSE)
    } else {
      datatable(
        tbl,
        selection = 'single',
        options = list(pageLength = 5, scrollY = '200px')
      )
    }
  })
  
  # Display selected outcome details
  output$selected_outcome <- renderPrint({
    sel <- input$outcome_table_rows_selected
    tbl <- outcomes()
    if (length(sel) && nrow(tbl) >= sel) {
      row <- tbl[sel, ]
      paste0("Selected Outcome ID: ", row$id, " (", row$trait, ")")
    } else {
      "No outcome selected"
    }
  })
  
  # Event-reactive that runs MR pipeline when button is clicked
  mr_pipeline <- eventReactive(input$run_mr, {
    # Ensure 'ao' is not empty
    validate(
      need(nrow(ao) > 0, "No outcome data available to run MR.")
    )
    # Get selected exposure and outcome IDs
    exp_sel <- input$exposure_table_rows_selected
    out_sel <- input$outcome_table_rows_selected
    validate(
      need(length(exp_sel) == 1, "Please select exactly one exposure."),
      need(length(out_sel) == 1, "Please select exactly one outcome.")
    )
    exposure_row <- exposures()[exp_sel, ]
    outcome_row <- outcomes()[out_sel, ]
    exposure_id <- exposure_row$id
    outcome_id <- outcome_row$id
    
    # Extract genome-wide significant instruments for exposure
    exposure_dat <- TwoSampleMR::extract_instruments(outcomes = exposure_id, p1 = 5e-08)
    # If no instruments are returned, throw validation error
    validate(
      need(nrow(exposure_dat) > 0, "No instruments found for selected exposure." )
    )
    
    # Extract outcome data for those SNPs
    outcome_dat <- TwoSampleMR::extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcome_id)
    validate(
      need(nrow(outcome_dat) > 0, "No outcome data found for selected outcome." )
    )
    
    # Harmonise exposure and outcome datasets
    dat_harmonised <- TwoSampleMR::harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
    # Perform MR analysis
    mr_res <- TwoSampleMR::mr(dat_harmonised)
    # Prepare scatter plot object
    scatter_plots <- TwoSampleMR::mr_scatter_plot(mr_res, dat_harmonised)
    
    list(mr_res = mr_res, scatter = scatter_plots)
  })
  
  # Render MR results table
  output$mr_results_table <- renderDT({
    res <- mr_pipeline()
    req(res)
    datatable(res$mr_res, options = list(pageLength = 5))
  })
  
  # Render MR scatter plot (first plot in the list)
  output$mr_scatter_plot <- renderPlot({
    res <- mr_pipeline()
    req(res)
    # mr_scatter_plot returns a list of ggplot objects; plot the first
    print(res$scatter[[1]])
  })
}

# Launch the Shiny app
shinyApp(ui = ui, server = server)
