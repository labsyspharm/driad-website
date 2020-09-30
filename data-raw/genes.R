library(tidyverse)

tmp_file <- tempfile()
download.file(
  "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2020-09-01.txt",
  tmp_file
)

genes <- read_tsv(
  tmp_file,
  col_types = cols_only(
    hgnc_id = "c",
    symbol = "c",
    name = "c",
    alias_symbol = "c",
    prev_symbol = "c",
    entrez_id = "i",
    ensembl_gene_id = "c"
  )
)

all_symbols <- genes %>%
  select(hgnc_id, symbol, alias_symbol, prev_symbol) %>%
  pivot_longer(-hgnc_id, names_to = "type", values_to = "symbol") %>%
  drop_na(symbol) %>%
  mutate(
    symbol = symbol %>%
      str_split(fixed("|"))
  ) %>%
  unchop(symbol)

usethis::use_data(genes, overwrite = TRUE, compress = "xz", version = 3)
