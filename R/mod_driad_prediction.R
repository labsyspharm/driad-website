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
}

#' Server module providing UI elements for DRIAD gene set evaluation
mod_ui_driad_prediction <- function(id) {
  ns <- NS(id)

  columns(
    column(
      width = c(lg = 6),
      card(
        header = tagList(
          p("asdf"),
          navInput(
            appearance = "tabs",
            id = ns("single_multi_choice"),
            choices = c("Single gene set", "Multiple gene sets"),
            values = c("single", "multi"),
            selected = "single",
            class = "card-header-tabs"
          )
        ),
        navContent(
          navPane(
            id = ns("pane_single"),
            fade = FALSE,
            formGroup(
              label = "Gene set",
              help = "Enter human gene symbols delimited by newlines, commas or tabs.",
              textAreaInput(
                inputId = ns("gene_set_single"),
                label = NULL,
                placeholder = "JAK3, BAX, EGFR, ...",
                rows = 10
              )
            )
          ),
          navPane(
            id = ns("pane_multi"),
            fade = FALSE,
            formGroup(
              label = "Gene set file",
              help = "CSV, TSV or Excel file with one gene set per column. First row should contain the gene set name.",
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
            )
          )
        )
      )
    ),
    column(
      width = c(lg = 6),
      p("asdf")
    )
  )
}
