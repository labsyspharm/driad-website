# Run on O2

library(DRIAD)
library(tidyverse)
library(here)
library(batchtools)

wd <- file.path("/n", "scratch3", "users", "c", "ch305", "driad")
dir.create(wd)

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
  here("data", paste0("prediction_tasks_all", ".rds"))
)

# prediction_tasks_all <- read_rds(here("data", paste0("prediction_tasks_all", ".rds")))

pwalk(
  prediction_tasks_all,
  function(task, pairs, id, ...) {
    # browser()
    x <- list(task = task, pairs = pairs)
    write_rds(
      x,
      here("data", paste0(id, ".rds"))
      # compress = "xz"
    )
  }
)

prediction_tasks <- prediction_tasks_all %>%
  select(-task, -pairs)

write_rds(
  prediction_tasks,
  here("data", paste0("prediction_tasks", ".rds")),
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

prediction_tasks_gene_sets <- prediction_tasks_all %>%
  crossing(
    gene_set_size = gene_set_sizes
  ) %>%
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

write_rds(
  prediction_tasks_gene_sets,
  here("data", paste0("prediction_tasks_gene_sets", ".rds"))
)

# prediction_tasks_gene_sets <- read_rds(here("data", paste0("prediction_tasks_gene_sets", ".rds")))

prediction_tasks_chunks <- prediction_tasks_gene_sets %>%
  mutate(
    task_file = normalizePath(here("data", paste0(id, ".rds")))
  ) %>%
  select(-task, -pairs) %>%
  split(rep_len(1:500, nrow(.))) %>%
  enframe("chunk", "data") %>%
  # slice(1) %>%
  mutate(
    input_file = file.path(wd, paste0("prediction_task_", chunk, ".rds"))
  )

pwalk(
  prediction_tasks_chunks,
  function(input_file, data, ...)
    write_rds(data, input_file)
)


# Set up jobs
reg <- makeRegistry(
  file.dir = file.path(wd, paste0("registry_", gsub(" ", "_", Sys.time()))),
  seed = 1,
  conf.file = here("data-raw", "batchtools-conf.R")
)
#reg$cluster.functions <- makeClusterFunctionsSlurm(template = "slurm-simple")

run_bk_job <- function(input_file, n_workers = 4) {
  library(tidyverse)
  library(DRIAD)
  library(furrr)

  message("Reading input...")
  df <- read_rds(input_file)
  message("Evaluating gene sets...")
  out <- df %>%
    rowwise() %>%
    mutate(
      task = read_rds(task_file) %>%
        list(),
      background_auc = DRIAD::evalGeneSets(
        background_sets,
        task[["task"]],
        task[["pairs"]]
      ) %>%
        list()
    ) %>%
    ungroup()
  message("Writing results...")
  write_rds(
    out %>%
      select(-task),
    file.path(
      dirname(input_file),
      paste0(tools::file_path_sans_ext(basename(input_file)), "_output.rds")
    )
  )
  message("Done...")
  input_file
}

batchMap(
  fun = run_bk_job,
  input_file = prediction_tasks_chunks[["input_file"]]
)

job_table <- findJobs() %>%
  # Chunk jobs into a single array job
  mutate(chunk = 1)

submitJobs(
  job_table,
  resources = list(
    memory = "8gb",
    ncpus = 1L,
    partition = "short",
    walltime = 5*60*60,
    chunks.as.arrayjobs = TRUE,
    # For some reason these nodes fail to execute R because of an "illegal instruction"
    exclude = "compute-f-17-[09-25]"
  )
)

prediction_tasks_outputs <- prediction_tasks_chunks %>%
  # rowwise() %>%
  mutate(
    result = map(
      paste0(tools::file_path_sans_ext(basename(input_file)), "_output.rds"),
      possibly(read_rds, otherwise = NULL)
    )
  )
  # ungroup()

background_auc_df <- prediction_tasks_outputs %>%
  # filter(map_lgl(result, Negate(is.null))) %>%
  pull(result) %>%
  bind_rows() %>%
  arrange(dataset, brain_region, comparison, gene_set_size) %>%
  select(dataset, brain_region, comparison, gene_set_size, background_auc) %>%
  mutate(
    background_auc = map(background_auc, pull, AUC)
  )

write_rds(
  background_auc_df,
  file.path(wd, "background_gene_sets_auc.rds")
)
