library(tidyverse)
library(DRIAD)
library(DT)

MAX_GENE_SETS <- 50
MIN_N <- 5
MAX_N <- 300

run_driad <- function(gene_sets, datasets, comparisons) {
  prediction_tasks %>%
    mutate(
      dataset = paste0(dataset, " - ", brain_region)
    ) %>%
    filter(dataset %in% datasets, comparison %in% comparisons) %>%
    rowwise() %>%
    mutate(
      task = prediction_task_data[[id]]["task"],
      pairs = prediction_task_data[[id]]["pairs"],
      res = evalGeneSets(
        gene_sets, task, pairs
      ) %>%
        list()
    ) %>%
    unnest(res) %>%
    select(-id)
}

#' Server module providing UI logic for DRIAD gene set evaluation
mod_server_driad_prediction <- function(
  input, output, session
) {
  ns <- session$ns

  observe({
    req(input$single_multi_choice)

    switch(
      input$single_multi_choice,
      single = showNavPane(ns("pane_single")),
      multi = showNavPane(ns("pane_multi"))
    )
  })

  observe({
    req(input$results_nav)

    switch(
      input$results_nav,
      plots = showNavPane(ns("pane_plots")),
      table = showNavPane(ns("pane_table"))
    )
  })

  r_gene_sets_raw <- reactive({
    req(input$single_multi_choice)

    switch(
      input$single_multi_choice,
      single = {
        req(input$gene_set_single)
        input$gene_set_single %>%
          trimws() %>%
          str_split("[\\s,;]+") %>%
          map(~.x[.x != ""]) %>%
          set_names("user_gene_set")
      },
      multi = {
        req(input$gene_set_upload)
        file <- input$gene_set_upload
        ext <- tools::file_ext(file$datapath)
        gene_sets <- switch(
          ext,
          "csv" = read_csv(file$datapath) %>%
            as.list(),
          "tsv" = read_tsv(file$datapath) %>%
            as.list(),
          "xlsx" = readxl::read_excel(file$datapath) %>%
            as.list(),
          "gmt" = DRIAD::read_gmt(file$datapath)
        )
        shiny::validate(
          need(
            length(gene_sets) <= MAX_GENE_SETS,
            paste("Number of gene sets must be", MAX_GENE_SETS, "or below.")
          )
        )
        gene_sets %>%
          map(na.omit)
      }
    )
  })

  r_gene_sets_valid <- reactive({
    req(r_gene_sets_raw())

    valid_mask <- r_gene_sets_raw() %>%
      # map(str_to_lower) %>%
      map(`%in%`, valid_gene_symbols)

    list(
      valid = map2(r_gene_sets_raw(), valid_mask, `[`),
      invalid = map2(r_gene_sets_raw(), map(valid_mask, `!`), `[`)
    )
  })

  output$gene_set_info <- renderUI({
    gene_sets <- tryCatch(
      r_gene_sets_valid(),
      error = function(e) NULL
    )
    if (is.null(gene_sets))
      return("No gene sets submitted.")
    gene_set_lengths <- map_int(r_gene_sets_valid()[["valid"]], length)
    invalid_symbols <- r_gene_sets_valid()[["invalid"]] %>%
      reduce(union)
    tagList(
      p(
        length(r_gene_sets_valid()[["valid"]]),
        "gene set(s) with between",
        min(gene_set_lengths), "and",
        max(gene_set_lengths), "genes."
      ),
      if (length(invalid_symbols) > 0)
        p(
          class = "text-warning",
          paste("Invalid gene symbols:", paste(invalid_symbols, collapse = ", "))
        ) %>%
          margin(b = 0)
    )
  })

  r_results <- reactiveVal()

  observeEvent(input$submit, {
    validate(
      need(length(input$datasets) > 0, "Must select at least one dataset.")
    )
    validate(
      need(
        any(map_lgl(r_gene_sets_valid()[["valid"]], ~length(.x) > 0)),
        "Must supply a gene set with at least one valid gene."
      )
    )
    run_driad(
      r_gene_sets_valid()[["valid"]],
      input$datasets,
      input$comparison
    ) %>%
      r_results()
  })

  output$results <- renderDT({
    datatable(
      if (!is.null(r_results()))
        select(r_results(), -task, -pairs),
      style = "bootstrap4",
      selection = "none"
    )
  })
}

#' Server module providing UI elements for DRIAD gene set evaluation
mod_ui_driad_prediction <- function(id) {
  ns <- NS(id)

  tagList(
    deck(
      card(
        header = tagList(
          h4("Gene set"),
          navInput(
            appearance = "tabs",
            id = ns("single_multi_choice"),
            choices = c("Single gene set", "Multiple gene sets"),
            values = c("single", "multi"),
            selected = "single",
            class = "card-header-tabs"
          )
        ),
        footer = uiOutput(outputId = ns("gene_set_info")),
        navContent(
          navPane(
            id = ns("pane_single"),
            fade = FALSE,
            formGroup(
              label = tags$label(style = "display: none;"),
              help = "Enter human gene symbols delimited by line breaks, commas, or tabs.",
              textAreaInput(
                inputId = ns("gene_set_single"),
                label = NULL,
                placeholder = "JAK3, BAX, EGFR, ...",
                rows = 10
              )
            ) %>%
              margin(b = 0)
          ),
          navPane(
            id = ns("pane_multi"),
            fade = FALSE,
            formGroup(
              label = tags$label(style = "display: none;"),
              help = "CSV, TSV or Excel file with one gene set per column. The first row should contain the gene set name.",
              fileInput(
                id = ns("gene_set_upload"),
                placeholder = "Choose gene set file",
                multiple = FALSE,
                accept = c(
                  "text/plain",
                  "text/csv",
                  "text/tab-separated-values",
                  "application/vnd.ms-excel",
                  "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                  ".gmt"
                )
              )
            ) %>%
              margin(b = 0)
          )
        )
      ),
      card(
        header = h4("Settings"),
        formGroup(
          label = "Comparison",
          help = paste(
            "Classification task to run.",
            "Can contrast early (A), intermediate (B) and late (C) disease stages,",
            "as defined by the Braak staging through neuropathological assessment"
          ),
          radiobarInput(
            id = ns("comparison"),
            label = NULL,
            choices = c("AB", "AC", "BC"),
            selected = "AC",
            class = "btn-group-primary"
          )
        ),
        formGroup(
          label = "Datasets",
          help = paste(
            "Datasets to run prediction task on."
          ),
          checkbarInput(
            id = ns("datasets"),
            label = NULL,
            choices = c("ROSMAP", "MSBB<br>BM10", "MSBB<br>BM22", "MSBB<br>BM36", "MSBB<br>BM44") %>%
              map(HTML),
            values = c("ROSMAP - Dorsal prefrontal cortex", "MSBB - BM10", "MSBB - BM22", "MSBB - BM36", "MSBB - BM44"),
            selected = c("ROSMAP - Dorsal prefrontal cortex", "MSBB - BM10", "MSBB - BM22", "MSBB - BM36", "MSBB - BM44"),
            class = "btn-group-primary"
          )
        ) %>%
          margin(b = 0)
      )
    ) %>%
      margin(b = 3),
    div(
      actionButton(
        inputId = ns("submit"),
        label = "Submit",
        class = "btn-primary"
      )
    ) %>%
      margin(b = 3),
    card(
      header = tagList(
        h4("Results"),
        navInput(
          id = ns("results_nav"),
          choices = c("Table", "Plots"),
          values  = c("table", "plots"),
          class = "card-header-tabs",
          appearance = "tabs",
          selected = "table"
        )
      ),
      navContent(
        navPane(
          id = ns("pane_table"),
          dataTableOutput(
            outputId = ns("results")
          )
        ),
        navPane(
          id = ns("pane_plots"),
          p("daff")
        )
      )
    )
  )
}
