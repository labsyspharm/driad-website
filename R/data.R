library(here)

prediction_tasks <- read_rds(here("data", "prediction_tasks.rda"))
valid_gene_symbols <- read_rds(here("data", "valid_gene_symbols.rda"))

# https://towardsdatascience.com/lazy-loading-data-in-r-2b100acb63fc
lazy.loader <- function(inputs, FUN) {
  obj <- new.env()  # define an new environment
  obj$storage <- list()
  obj$inputs <- inputs  # copy the inputs to the environment
  obj$FUN <- FUN  # user defines the function
  class(obj) <- "lazyloader"  # slap class label for method dispatch
  return(obj)  # return
}

`[[.lazyloader` <-  function(obj, idx, ...) {
  if (!idx %in% obj$inputs)
    stop("Invalid index", idx)
  if (!idx %in% names(obj$storage)) {
    obj$storage[[idx]] <- obj$FUN(idx, ...)
  }
  obj$storage[[idx]]
}

`names.lazyloader` <- function(obj) {
  obj$inputs
}

prediction_task_data <- lazy.loader(
  prediction_tasks[["id"]],
  function(x) {
    read_rds(here("data", paste0(x, ".rds")))
  }
)
