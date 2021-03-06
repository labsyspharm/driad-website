library(here)
library(tidyverse)
library(DRIAD)
library(DT)
library(promises)
library(future)
library(ggridges)
library(memoise)
library(ggplot2)
library(shinyjs)
library(ipc)
plan(multicore)
# plan(sequential)

MAX_GENE_SETS <- 50
MIN_N <- 5
MAX_N <- 300

BK_GENE_SET_SIZES <- c(
  5:29,
  seq(30, 300, by = 5)
)

DT_DOM <- paste0(
  '<"row justify-content-between"<"col-sm-12 col-md-auto"l><"col-sm-12 col-md-auto"B>',
  '<"col-sm-12 col-md-auto ml-md-auto"f>><"row"<"col-sm-12"t>><"row"<"col-sm-12 col-md-5"i>',
  '<"col-sm-12 col-md-7"p>>'
)

# JAK3, BCL2, TP53, EGFR, TYK2

prediction_tasks <- read_rds(here("data", "prediction_tasks.rda"))
valid_gene_symbols <- read_rds(here("data", "valid_gene_symbols.rda"))
background_gene_sets_auc <- read_rds(here("data", "background_gene_sets_auc.rds"))

prediction_task_data <- memoise(
  function(x) {
    read_rds(here("data", paste0(x, ".rds")))
  }
)

run_driad <- function(gene_sets, datasets, comparisons, progress = NULL) {
  has_progress <- !is.null(progress)
  if (has_progress)
    progress$set(value = 0, message = "Loading data")
  task_data <- prediction_tasks %>%
    filter(
      paste0(dataset, " - ", brain_region) %in% datasets,
      comparison %in% comparisons
    ) %>%
    rowwise() %>%
    mutate(
      task = prediction_task_data(id)["task"],
      pairs = prediction_task_data(id)["pairs"]
    ) %>%
    ungroup()
  if (has_progress)
    progress$set(value = 0.1, message = "Performing prediction tasks")
  # browser()
  tasks <- task_data %>%
    crossing(enframe(gene_sets, "Set", "Feats"))
  n_tasks <- nrow(tasks)
  tasks %>%
    mutate(
      res = pmap(
        .,
        function(Set, Feats, task, pairs, dataset, brain_region, ...) {
          if (has_progress)
            progress$inc(
              0.9/n_tasks,
              detail = paste0(Set, " - ", dataset, " ", brain_region)
            )
          evalGeneSet(Feats, task, pairs, nBK = 0) %>%
            select(-Feats, -Name)
        }
      )
    ) %>%
    unnest(res) %>%
    mutate(
      n_genes = map_int(Feats, length),
      n_genes_bk = BK_GENE_SET_SIZES[findInterval(n_genes, BK_GENE_SET_SIZES)]
    ) %>%
    inner_join(
      rename(background_gene_sets_auc, BK = background_auc),
      by = c("dataset", "brain_region", "comparison", "n_genes_bk" = "gene_set_size")
    ) %>%
    mutate(pval = map2_dbl(AUC, BK, ~`if`(length(.y) == 0, NA, mean(.x <= .y)))) %>%
    select(-id, -task, -pairs)
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
        if (is.null(input$gene_set_single) || input$gene_set_single == "")
          return(list("user_gene_set" = character()))
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
        gene_sets
      }
    ) %>%
      map(na.omit) %>%
      map(unique)
  })

  r_gene_sets_cleaned <- reactive({
    req(r_gene_sets_raw())

    valid_mask <- r_gene_sets_raw() %>%
      # map(str_to_lower) %>%
      map(`%in%`, valid_gene_symbols)

    list(
      valid = map2(r_gene_sets_raw(), valid_mask, `[`),
      invalid = map2(r_gene_sets_raw(), map(valid_mask, `!`), `[`)
    )
  })

  r_gene_sets_valid <- reactive({
    req(r_gene_sets_cleaned(), r_gene_sets_cleaned())

    r_gene_sets_cleaned()[["valid"]] %>%
      keep(~length(.x) >= MIN_N && length(.x) <= MAX_N)
  })

  output$gene_set_info <- renderUI({
    gene_set_lengths <- map_int(r_gene_sets_cleaned()[["valid"]], length)
    invalid_symbols <- r_gene_sets_cleaned()[["invalid"]] %>%
      reduce(union) %>%
      unique()
    tagList(
      if (length(r_input_message()) > 0)
        p(class = "text-danger", paste(r_input_message(), collapse = "<br>") %>% HTML()),
      if (length(r_gene_sets_valid()) > 0)
        p(length(r_gene_sets_valid()), "valid gene set(s)."),
      p(
        length(gene_set_lengths),
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

  r_input_message <- reactive({
    c(
      if (!length(input$datasets) > 0) "Must select at least one dataset.",
      if (!length(input$comparison) > 0) "Must select at least one comparison",
      if (!length(r_gene_sets_valid()) > 0) "No valid gene sets."
    )
  })

  r_results <- reactiveVal()

  observeEvent(input$submit, {
    input <- list(
      valid_gene_sets = r_gene_sets_valid(),
      datasets = input$datasets,
      comparisons = input$comparison
    )
    if (length(r_input_message()) > 0)
      return()
    disable(id = "submit")
    removeCssClass(class = "d-none", selector = paste0("#", ns("submit"), " .spinner-border"))
    p <- AsyncProgress$new()
    p$set(value = 0, message = "Starting...")
    res_future <- future(
      {
        run_driad(
          input[["valid_gene_sets"]],
          input[["datasets"]],
          input[["comparisons"]],
          progress = p
        )
      },
      packages = c("dplyr", "magrittr", "DRIAD"),
      globals = c(
        "prediction_task_data",
        "background_gene_sets_auc",
        "BK_GENE_SET_SIZES"
      ),
      seed = 1
    ) %...>%
      r_results() %>%
      finally(
        function() {
          p$close()
          enable("submit")
          addCssClass(class = "d-none", selector = paste0("#", ns("submit"), " .spinner-border"))
        }
      )
    NULL
  })

  r_plot <- reactive({
    if (is.null(r_results()))
      return(
        ggplot(tibble()) +
          annotation_custom(
            grid::textGrob(
              "No gene set submitted.",
              gp = grid::gpar(fontsize = 16)
            )
          ) +
          theme_void()
      )
    .data <- r_results() %>%
      mutate(dataset = paste0(dataset, " - ", brain_region)) %>%
      mutate(across(c(Set, dataset), as.factor))
    # browser()
    BK <- .data %>%
      select(Set, dataset, AUC = BK) %>%
      unnest(AUC)
    ggplot(BK, aes(x = AUC, y = Set)) +
      facet_wrap(~dataset) +
      theme_ridges() +
      geom_density_ridges2() +
      geom_segment(
        aes(x = AUC, xend = AUC, y = as.numeric(Set), yend = as.numeric(Set) + 0.9),
        data = .data,
        color = "red",
        lwd = 1
      ) +
      coord_cartesian(clip = "off")
  })

  output$plots <- renderPlot({
    req(r_plot())
    r_plot()
  })

  r_results_formatted <- reactive({
    if (is.null(r_results()))
      return(tibble(Set = character()))
    select(
      r_results(),
      gene_set = Set, dataset, brain_region, comparison, auc = AUC, pval,
      n_genes, n_genes_background = n_genes_bk, genes = Feats,
      auc_background = BK
    )
  })

  output$results <- DT::renderDataTable({
    req(r_results_formatted())
    .data <- r_results_formatted()
    DT::datatable(
      .data,
      style = "bootstrap4",
      selection = "none",
      extensions = "Buttons",
      rownames = FALSE,
      options = list(
        dom = DT_DOM,
        scrollX = TRUE,
        buttons = list(
          list(
            extend = "colvis",
            text = "Additional columns",
            className = "btn-outline-primary"
          )
        ),
        language = list(
          emptyTable = "No gene set submitted."
        ),
        columnDefs = list(
          list(
            targets = match(
              c("genes", "auc_background", "n_genes_background"),
              colnames(.data)
            ) %>%
              na.omit() %>%
              magrittr::subtract(1),
            visible = FALSE
          )
        )
      )
    )
  })

  r_zip_dl <- reactive({
    req(r_plot(), r_results_formatted())
    c(
      tempfile("plot", fileext = ".pdf") %T>%
        ggsave(r_plot()),
      tempfile("table", fileext = ".csv") %T>%
        {
          write_csv(
            mutate(r_results_formatted(), across(where(is.list), map_chr, paste, collapse = ";")),
            .
          )
        }
    )
  })

  callModule(mod_server_download_button, "csv", r_data = r_results_formatted, type = "csv")
  callModule(mod_server_download_button, "excel", r_data = r_results_formatted, type = "excel")
  callModule(mod_server_download_button, "zip", r_data = r_zip_dl, type = "zip")
}

#' Server module providing UI elements for DRIAD gene set evaluation
mod_ui_driad_prediction <- function(id) {
  ns <- NS(id)

  tagList(
    deck(
      card(
        header = tagList(
          h4("Gene set"),
          p(
            "Gene sets are evaluated for their ability to predict Alzheimer's",
            "disease progression from gene expression data."
          ),
          p("Gene sets must be between", MIN_N, "and", MAX_N, "genes in length."),
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
              help = "Enter human gene symbols delimited by line breaks, commas, spaces, or tabs.",
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
              help = paste(
                "CSV, TSV or Excel file with one gene set per column.",
                "The first row should contain the gene set name.",
                "GMT files are also accepted."
              ),
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
          checkbarInput(
            id = ns("comparison"),
            label = NULL,
            choices = c("AB", "AC", "BC"),
            selected = "AC",
            class = "btn-group-primary"
          )
        ),
        formGroup(
          label = "Datasets",
          help = tagList(
            p("AD datasets to run prediction task on.") %>%
              margin(b = 2),
            tags$ul(
              tags$li(
                "ROSMAP - The Religious Orders Study and Memory and Aging Project",
                tags$ul(
                  tags$li(
                    "Dorsal prefrontal cortex"
                  )
                )
              ),
              tags$li(
                "MSBB - The Mount Sinai/JJ Peters VA Medical Center Brain Bank",
                tags$ul(
                  c(
                    "BM10 - anterior prefrontal cortext",
                    "BM22 - superior temporal gyrus",
                    "BM36 - perirhinal cortext",
                    "BM44 - pars opercularis"
                  ) %>%
                    map(tags$li)
                )
              )
            )
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
        label = tagList(
          span(
            class = "spinner-border spinner-border-sm d-none",
            role = "status",
            `aria-hidden` = "true"
          ),
          "Submit"
        ),
        class = "btn-primary"
      ),
      textOutput(
        ns("submit_error"),
        inline = TRUE
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
          DT::dataTableOutput(
            outputId = ns("results")
          )
        ),
        navPane(
          id = ns("pane_plots"),
          plotOutput(
            outputId = ns("plots")
          ),
          p("The plots show a distribution of how well random gene sets perform",
            "in the AD prediction task. Random gene sets are chosen to be the same",
            "size as the query gene set. The red line indicates how well the",
            "user-submitted query gene set perfomed. It's position relative to",
            "random gene sets is an indicator of how predictive the query",
            "gene set is for AD progression.")
        )
      ),
      footer = tagList(
        mod_ui_download_button(ns("csv"), "Download CSV"),
        mod_ui_download_button(ns("excel"), "Download Excel"),
        mod_ui_download_button(ns("zip"), "Download table and plot")
      )
    ) %>%
      margin(b = 3)
  )
}
