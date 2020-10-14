library(shinyjs)

app_navbar <- navbar(
  # class = "navbar-dark bg-dark",
  a(
    href="http://sorger.med.harvard.edu/",
    img(src = "driad/img/logo.png", width = 32, height = 32),
    class = "navbar-brand mr-auto"
  ),
  div(
    class = "navbar-nav",
    linkInput(
      class = "nav-link",
      id = "about",
      label = "About"
    ) %>%
      margin(l = 2),
    linkInput(
      class = "nav-link",
      id = "funding",
      label = "Funding"
    ) %>%
      margin(l = 2),
    tags$a(
      class = "nav-link",
      href = "https://forms.gle/dSpCJSsbaavTbCkP6",
      target = "_blank",
      icon("comments", class = "fa-lg"),
      " Feedback"
    ) %>%
      margin(l = 2),
    tags$a(
      class = "nav-link",
      href = "https://github.com/labsyspharm/sms-website",
      target = "_blank",
      icon("github", class = "fa-lg")
    ) %>%
      margin(l = 2)
  )
) %>%
  margin(b = 5) %>%
  shadow()

app_main <- container(
  columns(
    column(
      width = 8,
      h1(
        "DRIAD - Drug Repurposing in Alzheimer's Disease"
      ),
      p(
        "The DRIAD tool evaluates gene sets for their ability to predict the",
        "progression of Alzheimer's disease (AD). It utilizes expression data from",
        "post-mortem samples of AD patients collected by the",
        a("AMP-AD consortium.", href = "https://adknowledgeportal.synapse.org", target = "_blank")
      ),
      h2(
        "Gene set AD progression prediction"
      ),
      mod_ui_driad_prediction("main")
    )
  ) %>%
    flex(align = "center", justify = "center")
)

app_ui <- function() {
  tagList(
    htmltools::htmlDependency(
      "font-awesome",
      "5.13.0", "www/shared/fontawesome", package = "shiny",
      stylesheet = c("css/all.min.css", "css/v4-shims.min.css")
    ),
    tags$head(
      tags$title("Drug Repurposing in Alzheimer's Disease"),
      tags$link(href="https://fonts.googleapis.com/css?family=Lato:400,700&display=swap", rel="stylesheet"),
      # tags$link(rel = "stylesheet", type = "text/css", href = "driad/css/driad-flatly.css"),
      tags$link(rel = "stylesheet", type = "text/css", href = "driad/css/main.css"),
      # tags$script(src = "sms/js/main.js"),
      # tags$link(rel = "icon", type = "image/png", href = "sms/assets/img/favicon.png"),
      useShinyjs()
    ),
    webpage(
      nav = app_navbar,
      app_main
    )
  )
}
