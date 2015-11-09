.jaccard  <- function(s1, s2) {
  length(intersect(s1, s2))/length(union(s1, s2))
}

.rebase_links <- function(nodes) {
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

.trim_links <- function(nodes, links, cutoff) {
    selection <- links %>%
      dplyr::group_by(source) %>%
      dplyr::select(node = source, jaccard) %>%
      dplyr::summarise(n = sum(jaccard >= cutoff)) %>%
      dplyr::arrange(desc(n)) %>%
      dplyr::filter(n >= 1) %>%
      dplyr::first(.)
    nodes <- nodes %>% dplyr::slice(selection + 1)
}

.create_colorscale <- function(nodes, palette) {
  cols <- RColorBrewer::brewer.pal(9, palette)[-c(1:2)]
  # cols <- c('blue', 'purple', 'red')
  cols <- paste0("'", paste(cols, collapse = "', '"), "'")
  range <- c(0, 5, 10, 50, 100, 200, 300)
  range <- paste(range, collapse = ", ")
  networkD3::JS(paste0('d3.scale.linear().domain([', range, ']).range([', cols, '])'))
}

clear <- function(leading_edge, cutoff = 0.25, trim = FALSE) {

  nodes <- leading_edge %>%
    dplyr::add_rownames('rowid') %>%
    dplyr::mutate(rowid = as.numeric(rowid) - 1,
                  group = -log10(fdr))

  links <- .rebase_links(nodes) %>%
    filter(jaccard >= cutoff)

  if (trim) {
    selection <- links %>%
      dplyr::group_by(source) %>%
      dplyr::select(node = source, jaccard) %>%
      dplyr::summarise(n = sum(jaccard >= cutoff)) %>%
      dplyr::arrange(desc(n)) %>%
      dplyr::filter(n >= 1) %>%
      dplyr::first(.)
    nodes <- nodes %>% dplyr::slice(selection + 1)
    links <- .rebase_links(nodes) %>% dplyr::filter(jaccard >= cutoff)
  }

  networkD3::forceNetwork(Links = links,
                          Nodes = nodes,
                          NodeID = 'geneset', Nodesize = 'n_geneset', Group = 'group',
                          Source = 'source', Target = 'target', Value = 'jaccard',
                          linkDistance = networkD3::JS("function(d) { return d.value * 100; }"),
                          # colourScale = JS("d3.scale.category20c()"),
                          colourScale = create_colorscale(nodes, 'BuPu'),
                          # colourScale = networkD3::JS("d3.scale.ordinal().domain(['foo', 'bar', 'baz']).range(colorbrewer.RdBu[9])"),
                          fontSize = 16, fontFamily = 'sans-serif', opacity = 0.75,
                          zoom = F, legend = T, bounded = T, opacityNoHover = 0.25,
                          charge = -300)
}
