library(shinyjs)

modal_funding <- modal(
  id = "modal_funding",
  header = "Funding",
  size = "sm",
  p(
    "We gratefully acknowledge support by NIA grant R01 AG058063: Harnessing",
    "Diverse BioInformatic Approaches to Repurpose Drugs for Alzheimers Disease."
  )
)

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
      id = "funding",
      label = "Funding"
    ) %>%
      margin(l = 2),
    tags$a(
      class = "nav-link",
      href = "https://github.com/labsyspharm/driad-website/issues",
      target = "_blank",
      icon("comments", class = "fa-lg"),
      " Feedback"
    ) %>%
      margin(l = 2),
    tags$a(
      class = "nav-link",
      href = "https://github.com/labsyspharm/driad-website",
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
      p(
        "The background and methodology are described in the preprint",
        a("Rodriguez and Hug et.al.", href = "https://doi.org/10.1101/2020.05.15.098749", target = "_blank")
      ),
      buttonInput(
        id = "explanation_toggle",
        label = "Show/hide details",
        class = "btn-primary"
      ),
      collapsePane(
        id = "explanation_pane",
        tags$figure(
          class = "figure",
          img(
            class = "figure-img img-fluid mx-auto d-block",
            style = "max-width: 40em;",
            src = "driad/img/machine_learning_framework-01.png"
          ),
          tags$figcaption(
            class = "figure-caption text-justify",
            "Overview of the machine learning framework used to establish potential",
            "associations between gene lists and Alzheimers Disease. (i) The",
            "framework accepts as input gene lists derived from experimental data",
            "or extracted from database resources or literature. (ii) Given a gene",
            "expression matrix, the framework subsamples it to a particular gene",
            "list of interest, and (iii) subsequently trains and evaluates through",
            "cross-validation a predictor of Braak stage of disease. (iv) The",
            "process is repeated for randomly-selected gene lists of equal lengths",
            "to determine whether predictor performance associated with the gene",
            "list of interest is significantly higher than whats expected by chance."
          )
        ),
        tags$figure(
          class = "figure",
          img(
            class = "figure-img img-fluid mx-auto d-block",
            style = "max-width: 25em;",
            src = "driad/img/brain_regions-01.png"
          ),
          tags$figcaption(
            class = "figure-caption text-justify",
            "AMP-AD datasets used by the machine learning framework. The three",
            "datasets used to evaluate the predictive power of gene lists are",
            "provided by The Religious Orders Study and Memory and Aging Project",
            "(ROSMAP), The Mayo Clinic Brain Bank (MAYO) and The Mount Sinai/JJ",
            "Peters VA Medical Center Brain Bank (MSBB). The schematic highlights",
            "regions of the brain that are represented in each dataset. The MSBB",
            "dataset spans four distinct regions, which are designated using",
            "Brodmann (BM) area codes."
          )
        ) %>%
          margin(b = 0),
        class = "border border-info rounded"
      ) %>%
        padding(a = 3),
      h2(
        "Gene set AD progression prediction"
      ) %>%
        margin(t = 3),
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
