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

clear <- function(leading_edge, cutoff = 0.25, trim = TRUE) {

  nodes <- leading_edge %>%
    dplyr::add_rownames('rowid') %>%
    dplyr::mutate(rowid = as.numeric(rowid) - 1)

  links <- .rebase_links(nodes) %>%
    filter(jaccard >= cutoff)

  if (trim) {
    selection <- links %>%
      dplyr::group_by(source) %>%
      dplyr::select(node = source, jaccard) %>%
      dplyr::summarise(n = sum(jaccard >= cutoff)) %>%
      dplyr::arrange(desc(n)) %>%
      dplyr::filter(n >= 5) %>%
      dplyr::first(.)
    nodes <- nodes %>% dplyr::slice(selection + 1)
    links <- .rebase_links(nodes) %>% dplyr::filter(jaccard >= cutoff)
  }

  networkD3::forceNetwork(Links = links,
                          Nodes = nodes,
                          NodeID = 'geneset', Nodesize = 'size', Group = 'subcollection',
                          Source = 'source', Target = 'target', Value = 'jaccard',
                          linkDistance = networkD3::JS("function(d) { return d.value * 100; }"),
                          fontSize = 24, fontFamily = 'sans-serif', opacity = 0.75,
                          zoom = F, legend = T, bounded = T, opacityNoHover = 0.25,
                          charge = -250)
}
