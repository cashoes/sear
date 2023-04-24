library(DBI)
library(RSQLite)
library(dplyr)

# downloaded from: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = 'data-raw/msigdb/msigdb_v2023.1.Hs.db')
DBI::dbListTables(con)

# define tables we want to combine
geneset_db <- dplyr::tbl(con, 'gene_set')                                              # standard_name, collection_name
details_db <- dplyr::tbl(con, 'gene_set_details')                                      # description_brief, description_full
geneset_genesymbol_db <- dplyr::tbl(con, 'gene_set_gene_symbol')                       # meat and potatoes
genesymbol_db <- dplyr::tbl(con, 'gene_symbol')                                        # mapping from ids to gene symbols
collection_db <- dplyr::tbl(con, 'collection') %>% dplyr::select(collection_name, full_name)  # collection metadata

# join tables
msigdb <- geneset_db %>%
  dplyr::left_join(details_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(collection_db, by = 'collection_name') %>%
  dplyr::left_join(geneset_genesymbol_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(genesymbol_db, by = c('gene_symbol_id' = 'id'))

# cleanup and make tibble
msigdb <- msigdb %>%
  dplyr::select(collection = collection_name, subcollection = full_name, geneset = standard_name, description = description_brief, symbol) %>%
  dplyr::as_tibble() %>%
  tidyr::nest(members_mrna = symbol) %>%
  dplyr::mutate(members_mrna = purrr::map(members_mrna, 'symbol'))


attributes(msigdb)$version <- 'v2023.1.Hs'

# clean up
DBI::dbDisconnect(con)
rm(con, geneset_db, details_db, collection_db, geneset_genesymbol_db, genesymbol_db)
