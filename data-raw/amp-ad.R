library(DRIAD)

fnROSMAP <- wrangleROSMAP(tempdir())
fnMSBB <- Sys.glob(here::here("data-raw", "msbb*.tsv.gz"))

prediction_tasks <- tribble(
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
  ungroup()


gene_symbols <- prediction_tasks %>%
  rowwise() %>%
  mutate(
    task = if (dataset == "MSBB") select(task, -(1:6)) %>% list() else select(task, -(1:7)) %>% list()
  ) %>%
  pull(task) %>%
  map(colnames) %>%
  reduce(union)

usethis::use_data(prediction_tasks, overwrite = TRUE, internal = FALSE, compress = "xz", version = 3)
