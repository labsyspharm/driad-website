library(here)

cluster.functions <- makeClusterFunctionsSlurm(template = here("data-raw", "slurm-o2-template.tmpl"))
