#' Simple enrichment analysis in R using the MSigDB collections.
#'
#' Carry out enrichment analysis against some or all of the MSigDB collections
#' (1), Blood Transcriptome genes (2) and Tissue Enrichment Sets (3).
#'
#' A list of genes is compared to each annotated gene set in turn by performing
#' a hypergeometric test of the overlap. The size of the input gene list, gene
#' set, intersection and resulting p-value are returned. P-values can be
#' adjusted over each collection (optimistic) or globally (pessimistic), using
#' any of the methods available to `p/adjust()`.
#'
#' @param input A character vector of HUGO gene naming commitee (HGNC) gene
#'   symbols.
#' @param type A character string. The type of symbols in the input: 'mrna' for
#'   gene symbols, 'mirna' for miRNAs. The reference gene sets used differ for
#'   either types. See: Godard, P., and van Eyll, J. (2015). Pathway analysis
#'   from lists of microRNAs: common pitfalls and alternative strategy. Nucl.
#'   Acids Res. 43, 3490-3497.
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#'
#' library(dplyr)
#' data('genesets')
#' input <- genesets %>%
#'   filter(collection == 'tissues') %>%
#'   nth(4) %>%
#'   unlist() %>%
#'   unique() %>%
#'   sample(100)
#' output <- sear(input, type = 'mrna') %>%
#'   arrange(fdr) %>%
#'   slice(1:100)
sear <- function(input, type = c("mrna", "mirna")) {
  type <- match.arg(type)
  tbl <- switch(type,
                mrna   = dplyr::select(genesets, collection:geneset,
                                       members = members_mrna),
                mirna  = dplyr::select(genesets, collection:geneset,
                                       members = members_mirna))
  uni <- tbl$members %>% unlist() %>% unique()

  # check type of input
  if (sum(input %in% uni)/length(input) < 0.25)
    stop(sprintf("Type = '%s', but many features are not recognized.", type))

  # only keep valid symbols
  input <- input[input %in% uni]

  # check size of input
  if (length(unique(input)) < 10)
    warning("Submitted <10 valid symbols. Results may not be meaningful.")

  tbl %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(n_input   = length(input),
                  n_geneset = length(members),
                  intersect = length(intersect(input, members)),
                  p_value   = phyper(intersect - 1,
                                     n_input,
                                     length(uni) - n_input,
                                     n_geneset, lower.tail = F)) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(collection) %>%
    dplyr::mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup(.)
}
