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
#' @param return_members A logical. Whether to include a list-column in the
#'   return object containing the gene set member genes.
#' @param return_intersect A logical. Whether to include a list-column in the
#'   return object containing the genes found in common between the gene set and
#'   the input.
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#'
#' library(dplyr)
#' data('collections')
#' input <- collections$members_mrna %>% unlist() %>% unique() %>% sample(100)
#' output <- sear(input, type = 'mrna') %>%
#'   arrange(fdr) %>%
#'   slice(1:100)
sear <- function(input, type = c("mrna", "mirna"), return_members = F, return_intersect = T) {
  data("collections", envir = environment())

  type <- match.arg(type)
  tbl <- switch(
    type,
    mrna   = dplyr::select(collections, collection:geneset, members = members_mrna),
    mirna  = dplyr::select(collections, collection:geneset, members = members_mirna)
  )
  uni <- tbl$members %>% unlist() %>% unique()

  # warn on potential input issues
  recognized <- input[input %in% uni]
  if (length(recognized) == 0) {
    stop("No symbols were found in gene set collection.")
  }
  if (length(recognized) < length(input)) {
    warning(sprintf("Submitted %s symbols, but only %s are recognized.", length(input), length(recognized)))
  }
  input <- recognized

  tbl <- tbl %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
      n_input   = length(input),
      n_geneset = length(members),
      intersect = list(intersect(input, members)),
      p_value   = phyper(
        length(intersect) - 1,
        n_geneset,
        length(uni) - n_geneset,
        n_input, lower.tail = F
      )
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(collection) %>%
    dplyr::mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup(.)

  # remove members
  if(!return_members)
    tbl <- tbl %>% dplyr::select(-members)

  # remove intersect
  if(!return_intersect)
    tbl <- tbl %>% dplyr::select(-intersect)

  return(tbl)
}
