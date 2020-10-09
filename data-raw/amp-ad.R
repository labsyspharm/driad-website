library(DRIAD)
library(tidyverse)
library(here)
library(furrr)
library(carrier)

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

write_rds(
  prediction_tasks_all,
  here("data", paste0("prediction_tasks_all", ".rda")),
  compress = "xz"
)

pmap(
  prediction_tasks_all,
  function(task, pairs, id, ...) {
    # browser()
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

gene_set_sizes <- c(
  5:29,
  seq(30, 300, by = 5)
)

prediction_tasks_bk_gene_sets <- prediction_tasks_all %>%
  mutate(
    data = tibble(gene_set_size = gene_set_sizes) %>%
      list()
  ) %>%
  mutate(
    data = pmap(
      list(task, data),
      function(task, df) {
        df %>%
          rowwise() %>%
          mutate(
            background_sets = DRIAD:::genBK(
              valid_gene_symbols[1:gene_set_size],
              task,
              1000
            ) %>%
              set_names(., paste0("BK_", seq_along(.))) %>%
              list()
          ) %>%
          ungroup()
      }
    )
  )
# %>%
  # mutate(
  #   id = paste0("prediction_task_background_", 1:n())
  # )

plan(sequential(workers = 1))

eval_gene_sets <- crate(
  function(task, pairs, df) {
    df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        background_auc = DRIAD::evalGeneSets(
          background_sets,
          task,
          pairs
        ) %>%
          list()
      ) %>%
      dplyr::ungroup()
  },
  `%>%` = `%>%`
)

plan(multisession(workers = 8))
# plan(sequential)
prediction_tasks_bk_auc <- prediction_tasks_bk_gene_sets %>%
  mutate(
    data = pmap(
      list(
        task,
        pairs,
        data
      ),
      eval_gene_sets
    ),
    # .progress = TRUE,
    # options = future_options(globals = FALSE, scheduling = FALSE)
  )
#
# pmap(
#   prediction_tasks_bk,
#   function(background_sets, id, ...) {
#     # browser()
#     assign(id, background_sets)
#     save(list = id, file = here("data", paste0(id, ".rda")), version = 3, compress = "xz")
#   }
# )
