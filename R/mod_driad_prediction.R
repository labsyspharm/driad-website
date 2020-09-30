library(tidyverse)

MAX_GENE_SETS <- 50

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

  r_gene_sets_raw <- reactive({
    req(input$single_multi_choice)

    switch(
      input$single_multi_choice,
      single = {
        req(input$gene_set_single)
        input$gene_set_single %>%
          trimws() %>%
          str_split("[\\s,;]+") %>%
          set_names("user_gene_set")
      },
      multi = {
        req(input$gene_set_upload)
        file <- input$gene_set_upload
        ext <- tools::file_ext(file$datapath)
        df <- switch(
          ext,
          "csv" = read_csv(file$datapath),
          "tsv" = read_tsv(file$datapath),
          "xlsx" = readxl::read_excel(file$datapath)
        )
        shiny::validate(
          need(
            ncol(df) <= MAX_GENE_SETS,
            paste("Number of gene sets must be", MAX_GENE_SETS, "or below.")
          )
        )
        df %>%
          as.list() %>%
          map(na.omit)
      }
    )
  })

  output$gene_set_info <- renderUI({
    gene_sets <- tryCatch(
      r_gene_sets_raw(),
      error = function(e) NULL
    )
    if (is.null(gene_sets))
      return("No gene sets submitted.")
    gene_set_lengths <- map_int(r_gene_sets_raw(), length)
    paste(
      length(r_gene_sets_raw()),
      "gene set(s) with between",
      min(gene_set_lengths), "and",
      max(gene_set_lengths), "genes"
    )
  })
}

#' Server module providing UI elements for DRIAD gene set evaluation
mod_ui_driad_prediction <- function(id) {
  ns <- NS(id)

  columns(
    column(
      width = c(lg = 6),
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
                  "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
              )
            ) %>%
              margin(b = 0)
          )
        )
      )
    ),
    column(
      width = c(lg = 6),
      card(
        header = h4("Settings"),
        footer = div(),
        formGroup(
          label = "Comparison",
          help = paste(
            "Classification task to run.\n",
            "Can contrast early (A), intermediate (B) and late (C) disease stages,",
            "as defined by the Braak staging through neuropathological assessment"
          ),
          radiobarInput(
            id = ns("comparison"),
            label = NULL,
            choices = c("AB", "AC", "BC"),
            selected = "AC"
            # inline = TRUE
          )
        )
      )
    )
  )
}
