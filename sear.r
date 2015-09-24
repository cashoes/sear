#' Get a glimpse of your data.
#'
#' Carry out enrichment analysis against some or all of the MSigDB collections
#' (1), Blood Transcriptome Modules (2) and Tissue Enrichment Sets (3).
#'
#' A list of genes is compared to each annotated gene set in turn and a
#' hypergeometric test of the overlap is performed. The size of the input gene
#' list, gene set, intersection and resulting p-value are returned. P-values can
#' be adjusted over each collection (optimistic) or globally (pessimistic),
#' using any of the methods available to `p/adjust()`.
#'
#' @param genes A character vector of HUGO gene naming commitee (HGNC) gene symbols.
#' @param references A data.frame of annotated gene sets: defaults to a prepared
#'   aggregate of the collections described above.
#' @export
#' @examples
#' sear(c("ACTB", "B2M", "SDHA", "LTBR", "HBB"), references)
sear <- function(module, references) {
  references %>%
    rowwise() %>%
    mutate(n_module = length(module),
           n_geneset = length(members),
           intersect = length(intersect(module, members)),
           p_value = 1 - phyper(intersect - 1,
                                n_module,
                                UNIVERSE_MRNA - n_module,
                                n_geneset)) %>%
    ungroup() %>%
    select(-members)
}
