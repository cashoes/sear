#' Clear
#'
#' @param leading_edge The output of running sear on an input list of mRNAs or
#'   miRNA symbols.
#' @param cutoff A numeric value between 0 and 1. The minimum value of the
#'   jaccard coefficient to show links between genesets.
#' @param trim A logical value. Whether to perform trimming of non-connected
#'   nodes.
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
#' clear(output, cutoff = 0.25, trim = TRUE)
clear <- function(leading_edge, cutoff = 0.25, trim = FALSE) {

  nodes <- leading_edge %>%
    tibble::rownames_to_column('rowid') %>%
    dplyr::mutate(rowid = as.numeric(rowid) - 1,
                  group = -log10(fdr + .Machine$double.xmin))

  links <- rebase_links(nodes) %>%
    filter(jaccard >= cutoff)

  if (trim) {
    selection <- c(links$source, links$target) %>%
      unique() %>%
      match(., (1:nrow(nodes) - 1), nomatch = NA) %>%
      sort()
    nodes <- nodes %>%
      dplyr::slice(selection) %>%
      dplyr::mutate(rowid = (0:(n() - 1)))
    links <- rebase_links(nodes) %>% dplyr::filter(jaccard >= cutoff)
  }

  networkD3::forceNetwork(
    Links = links,
    Nodes = nodes,
    NodeID = 'geneset', Nodesize = 'n_geneset', Group = 'group',
    Source = 'source', Target = 'target', Value = 'jaccard',
    linkDistance = networkD3::JS("function(d) { return d.value * 100; }"),
    colourScale = create_colorscale(nodes, 'BuPu'),
    fontSize = 16, fontFamily = 'sans-serif', opacity = 0.85,
    zoom = F, legend = T, bounded = T, opacityNoHover = 0.00,
    charge = -300)
}

jaccard  <- function(s1, s2) {
  length(intersect(s1, s2))/length(union(s1, s2))
}

rebase_links <- function(nodes) {
  ref <- nodes$members
  names(ref) <- nodes$rowid
  t(combn(as.numeric(names(ref)), 2)) %>%
    as.data.frame(.) %>%
    dplyr::tbl_df(.) %>%
    dplyr::rename(source = V1, target = V2) %>%
    dplyr::mutate(jaccard  = unlist(purrr::map2(ref[as.character(source)],
                                                ref[as.character(target)],
                                                jaccard))) %>%
    dplyr::mutate(source = match(source, nodes$rowid) - 1,
                  target = match(target, nodes$rowid) - 1)
}

create_colorscale <- function(nodes, palette) {
  cols <- RColorBrewer::brewer.pal(9, palette)[-c(1:2)]
  # cols <- c('blue', 'purple', 'red')
  cols <- paste0("'", paste(cols, collapse = "', '"), "'")
  # range <- c(0, 5, 10, 50, 100, 200, -log10(.Machine$double.xmin) %>% ceiling())
  max <- -log10(nodes$fdr) %>% max() %>% ceiling()
  range <- quantile(0:max) %>% ceiling()
  range <- paste(range, collapse = ", ")
  networkD3::JS(paste0('d3.scale.linear().domain([', range, ']).range([', cols, '])'))
}
