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
#' @param genes A character vector of HUGO gene naming commitee (HGNC) gene
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
#' input <- c("TTLL4", "CTSK", "NNMT", "TGFBR2", "MAGED2", "ASB13", "CCDC80", "APBB2", "RABEP1", "FBP1")
#' sear(input, type = 'mrna')
sear <- function(genes, type = c("mrna", "mirna")) {
  type <- match.arg(type)

  tbl <- switch(type,
                mrna  = genesets %>% dplyr::select(collection, subcollection, geneset, members = members_mrna),
                mirna = genesets %>% dplyr::select(collection, subcollection, geneset, members = members_mirna))
  uni <- tbl$members %>% unlist() %>% unique()

  # check type of input
  if (sum(genes %in% uni)/length(genes) < 0.25)
    stop(sprintf("You selected type = '%s', but many of your features are not recognized.", type))

  # only keep valid symbols
  genes <- genes[genes %in% uni]

  # check size of input
  if (length(unique(genes)) < 10)
    warning("You submitted <10 valid symbols. Results may not be meaningful with so few genes.")

  tbl %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(n_genes = length(genes),
                  n_geneset = length(members),
                  intersect = length(intersect(genes, members)),
                  p_value = 1 - phyper(intersect - 1,
                                       n_genes,
                                       length(uni) - n_genes,
                                       n_geneset)) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-members) %>%
    dplyr::group_by(collection) %>%
    dplyr::mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup(.)
}

.jaccard  <- function(s1, s2) {
  length(intersect(s1, s2))/length(union(s1, s2))
}

.process_nodes <- function(nodes) {
  nodes %>%
    dplyr::add_rownames('rowid') %>%
    dplyr::mutate(rowid = as.numeric(rowid) - 1,
                  size = unlist(purrr::map(members, length))) %>%
    dplyr::select(rowid, geneset, collection, subcollection, size, members)
}

.process_links <- function(nodes) {
  ref <- nodes$members
  names(ref) <- nodes$rowid
  t(combn(as.numeric(names(ref)), 2)) %>%
    as.data.frame(.) %>%
    dplyr::tbl_df(.) %>%
    dplyr::rename(source = V1, target = V2) %>%
    dplyr::mutate(jaccard  = unlist(map2(ref[as.character(source)], ref[as.character(target)], .jaccard))) %>%
    dplyr::mutate(source = match(source, nodes$rowid) - 1,
                  target = match(target, nodes$rowid) - 1)
}

clear <- function(leading_edge, type = c('mrna', 'mirna'), cutoff = 0.25) {
  type <- match.arg(type)
  tbl <- switch(type,
                mrna  = genesets %>% dplyr::select(collection, subcollection, geneset, members = members_mrna),
                mirna = genesets %>% dplyr::select(collection, subcollection, geneset, members = members_mirna))

  nodes <- .process_nodes(tbl) %>%
    filter(geneset %in% leading_edge)

  links <- .process_links(nodes) %>%
    filter(jaccard >= cutoff)

  networkD3::forceNetwork(Links = links, #links_small,
                          Nodes = nodes, #nodes_small,
                          NodeID = 'geneset', Nodesize = 'size', Group = 'subcollection',
                          Source = 'source', Target = 'target', Value = 'jaccard',
                          linkDistance = networkD3::JS("function(d) { return d.value * 100; }"),
                          fontSize = 24, fontFamily = 'sans-serif', opacity = 0.75,
                          zoom = F, legend = T, bounded = T, opacityNoHover = 0.25,
                          charge = -250)
}
