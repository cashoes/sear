#' Add additional geneset collections to MSigDB data.frame
#' @author C.Shannon
#'
#'   Currently only blood transcriptome modules:
#'
#'   Chaussabel D, Baldwin N. Democratizing systems immunology with modular
#'   transcriptional repertoire analyses. Nat Rev Immunol. 2014;14: 271–280.
#'   doi:10.1038/nri3642
#'
#'   Li S, Rouphael N, Duraisingham S, Romero-Steiner S, Presnell S, Davis C, et
#'   al. Molecular signatures of antibody responses derived from a systems
#'   biology study of five human vaccines. Nature Immunology. 2013;15: 195–204.
#'   doi:10.1038/ni.2789
#'
#'   and tissue specific genesets derived from:
#'
#'   Benita Y, Cao Z, Giallourakis C, Li C, Gardet A, Xavier RJ. Gene
#'   enrichment profiles reveal T-cell development, differentiation, and
#'   lineage-specific transcription factors including ZBTB25 as a novel NF-AT
#'   repressor. Blood. 2010;115: 5376–5384. doi:10.1182/blood-2010-01-263855

# helper: parse gmt files ----
parse_gmt <- function(gmt){
  geneset <- readLines(gmt)
  geneset <- sapply(geneset, function(x) strsplit(x, split = "\t"))
  names(geneset) <- lapply(geneset, function(x) head(x, 1))
  geneset <- lapply(geneset, function(x) tail(x, -2))
}


# BTM as gmt ----
btm <- parse_gmt('data-raw/misgdb/btm.all.v1.0.symbols.gmt')
data.frame(collection = 'BTM',
           subcollection = '',
           geneset = names(btm),
           description = '',
           members_mrna = I(btm),
           stringsAsFactors = F) %>%
  dplyr::tbl_df() -> btms

rm(btm)

# Benita et al. ----
tissues <- readRDS("data-raw/Tissues_tbl_df.rds")
tissues_map <- readRDS("data-raw/Tissues_map_tbl_df.rds")
all(tissues_map$celltype == tissues$geneset)

tissues <- tissues %>%
  dplyr::mutate(subcollection = toupper(tissues_map$category),
                description = '') %>%
  dplyr::select(collection, subcollection, geneset, description, members_mrna = members)
rm(tissues_map)

# Combine ----
collections <- rbind(msigdb, tissues, btms)
collections <- collections %>% arrange(collection, subcollection, geneset)
