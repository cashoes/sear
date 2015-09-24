#' Simple enrichment analysis in R using the MSigDB collections.
#'
#' Carry out enrichment analysis against some or all of the MSigDB collections
#' (1), Blood Transcriptome geness (2) and Tissue Enrichment Sets (3).
#'
#' A list of genes is compared to each annotated gene set in turn by performing
#' a hypergeometric test of the overlap. The size of the input gene list, gene
#' set, intersection and resulting p-value are returned. P-values can be
#' adjusted over each collection (optimistic) or globally (pessimistic), using
#' any of the methods available to `p/adjust()`.
#'
#' @param genes A character vector of HUGO gene naming commitee (HGNC) gene
#'   symbols.
#' @param genesets A data.frame of annotated gene sets: defaults to a prepared
#'   aggregate of the collections described above.
#' @export
#' @examples
#' sear(c("ACTB", "B2M", "SDHA", "LTBR", "HBB"), genesets)
sear <- function(genes) {
  genesets %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n_genes = length(genes),
           n_geneset = length(members),
           intersect = length(intersect(genes, members)),
           p_value = 1 - phyper(intersect - 1,
                                n_genes,
                                UNIVERSE_MRNA - n_genes,
                                n_geneset)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-members)
}

UNIVERSE_MRNA <- 45956
UNIVERSE_MIRNA <- 1881
