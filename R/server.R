
app_server <- function(input, output, session) {
  callModule(mod_server_driad_prediction, "main")
  observeEvent(input$explanation_toggle, {
    toggleCollapsePane("explanation_pane")
  })
  observeEvent(input$funding, {
    showModal(modal_funding)
  })
}
