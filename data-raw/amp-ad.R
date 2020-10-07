library(DRIAD)
library(tidyverse)
library(here)

fnROSMAP <- wrangleROSMAP(tempdir())
# fnMSBB <- Sys.glob(here::here("data-raw", "msbb*.tsv.gz"))
fnMSBB <- wrangleMSBB(tempdir())

prediction_tasks_all <- tribble(
  ~brain_region, ~dataset, ~path,
  "Dorsal prefrontal cortex", "ROSMAP", fnROSMAP,
  "BM10", "MSBB", fnMSBB[1],
  "BM22", "MSBB", fnMSBB[2],
  "BM36", "MSBB", fnMSBB[3],
  "BM44", "MSBB", fnMSBB[4],
) %>%
  crossing(
    comparison = c("AB", "AC", "BC")
  ) %>%
  rowwise() %>%
  mutate(
    task = prepareTask(path, comparison) %>%
      list(),
    pairs = preparePairs(task) %>%
      list()
  ) %>%
  select(-path) %>%
  ungroup() %>%
  mutate(
    id = paste0("prediction_task_", 1:n())
  )

pmap(
  prediction_tasks_all,
  function(task, pairs, id, ...) {
    browser()
    x <- list(task = task, pairs = pairs)
    write_rds(
      x,
      here("data", paste0(id, ".rds")),
      compress = "xz"
    )
  }
)

prediction_tasks <- prediction_tasks_all %>%
  select(-task, -pairs)

write_rds(
  prediction_tasks,
  here("data", paste0("prediction_tasks", ".rda")),
  compress = "xz"
)

valid_gene_symbols <- prediction_tasks_all %>%
  rowwise() %>%
  mutate(
    task = if (dataset == "MSBB") select(task, -(1:6)) %>% list() else select(task, -(1:7)) %>% list()
  ) %>%
  pull(task) %>%
  map(colnames) %>%
  reduce(intersect) %>%
  setdiff(
    c("ID", "PMI", "AOD", "CDR",
      "Braak", "Barcode", "Label")
  )

write_rds(
  valid_gene_symbols,
  here("data", paste0("valid_gene_symbols", ".rda")),
  compress = "xz"
)

# gene_set_sizes <- c(
#   5:29,
#   seq(30, 300, by = 5)
# )
#
# prediction_tasks_bk <- prediction_tasks %>%
#   select(-pairs) %>%
#   crossing(
#     gene_set_size = gene_set_sizes
#   ) %>%
#   rowwise() %>%
#   mutate(
#     background_sets = DRIAD:::genBK(
#       gene_symbols[1:gene_set_size],
#       task,
#       1000
#     ) %>%
#       list()
#   ) %>%
#   ungroup() %>%
#   select(-task) %>%
#   mutate(
#     id = paste0("prediction_task_background_", 1:n())
#   )
#
# pmap(
#   prediction_tasks_bk,
#   function(background_sets, id, ...) {
#     # browser()
#     assign(id, background_sets)
#     save(list = id, file = here("data", paste0(id, ".rda")), version = 3, compress = "xz")
#   }
# )
