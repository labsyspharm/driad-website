library(here)

load(here("data", "prediction_tasks.rda"))

valid_gene_symbols <- prediction_tasks %>%
  rowwise() %>%
  mutate(
    task = if (dataset == "MSBB") select(task, -(1:6)) %>% list() else select(task, -(1:7)) %>% list()
  ) %>%
  pull(task) %>%
  map(colnames) %>%
  reduce(intersect)
  # str_to_lower()
