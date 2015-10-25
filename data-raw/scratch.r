# scratch.r

# write a function to compute Jaccard given a pair of 'genesets'
jaccard  <- function(s1, s2) length(intersect(s1, s2))/length(union(s1, s2))
# jaccard <- function(s1, s2) length(intersect(s1, s2))/min(length(s1), length(s2))

# visualize
# install.packages('networkD3')
library(networkD3)

# networkD3 takes 2 data.frames to render network:
#
# Links: a data frame object
# with the links between the nodes. It should include the Source and Target for
# each link. These should be numbered starting from 0. An optional Value
# variable can be included to specify how close the nodes are to one another.
#
# Nodes: a data frame containing the node id and properties of the nodes. If no
# ID is specified then the nodes must be in the same order as the Source
# variable column in the Links data frame. Currently only a grouping variable is
# allowed.

# # selection from SGCCDA manuscript
# selection <- readClipboard()

# using purrr and map2 instead - about half the time...
library(purrr)

# prepare nodes and links data.frames for use by networkD3 ---------------------
data('genesets')

# NodeID = 'geneset', Nodesize = 'size', Group = 'collection',
nodes <- genesets %>%
  dplyr::add_rownames('rowid') %>%
  dplyr::mutate(size = unlist(purrr::map(members, length))) %>%
  dplyr::select(rowid, geneset, collection, subcollection, size, members)

rm(genesets)

# list for lookup - numbered version
foo <- nodes$members
names(foo) <- 0:(nrow(nodes) - 1)

# Source = 'source', Target = 'target', Value = 'jaccard',
system.time({
  links <- data.table::CJ(as.numeric(names(foo)), as.numeric(names(foo))) %>%
    dplyr::tbl_df(.) %>%
    dplyr::rename(source = V1, target = V2) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(jaccard  = jaccard(unlist(foo[source + 1]), unlist(foo[target + 1]))) %>%
    # dplyr::mutate(jaccard  = unlist(map2(foo[source + 1], foo[target + 1], jaccard))) %>%
    dplyr::filter(source != target) # remove self-references
})

rm(foo)
save(nodes, links, file = 'data-raw/genesets_adjacency.rda')

# come up with a way to filter from the whole table of links and
load('data-raw/genesets_adjacency.rda')
library(networkD3)

# given a filtered nodes data.frames subset links and return numbered version
# suitable for use with networkD3
# rebase <- function(nodes, links) {
#   # subset
#   # index <- which(links$source %in% nodes$rowid & links$target %in% nodes$rowid)
#   links %>%
#     dplyr::filter(source %in% nodes$rowid, target %in% nodes$rowid) %>%
#     # dplyr::slice(index) %>%
#     dplyr::mutate(source = rep(1:nrow(nodes), each = nrow(nodes)),
#            target = rep(1:nrow(nodes), nrow(nodes)))
# }

rebase <- function(nodes) {
  foo <- nodes$members
  names(foo) <- 0:(nrow(nodes) - 1)

  data.table::CJ(as.numeric(names(foo)), as.numeric(names(foo))) %>%
    dplyr::tbl_df(.) %>%
    dplyr::rename(source = V1, target = V2) %>%
    dplyr::mutate(jaccard  = unlist(map2(foo[source + 1], foo[target + 1], jaccard))) %>%
    dplyr::filter(source != target) # remove self-references
}

nodes_small <- nodes %>% filter(collection == 'TISSUES')
system.time(links_small <- rebase(nodes_small) %>% filter(jaccard > 0.15))
forceNetwork(Links = links_small,
             Nodes = nodes_small,
             NodeID = 'geneset', Nodesize = 'size', Group = 'subcollection',
             Source = 'source', Target = 'target', Value = 'jaccard',
             fontSize = 10, fontFamily = 'sans', opacity = 0.75,
             zoom = F, legend = T, bounded = F, opacityNoHover = 0,
             charge = -120, linkDistance = 50)

# # try ggnet2: omg so slow... -------------------------------------------------
# # devtools::install_github("briatte/ggnet")
# # install.packages(c('sna', 'network'))
# library(ggplot2)
# library(grid)
# library(scales)
#
# library(sna)
# library(network)
# library(ggnet)
#
# # example ggnet2
# # node information
# ids = read.csv("https://github.com/briatte/ggnet/raw/master/inst/extdata/nodes.tsv", sep = "\t")
# # edge list
# df = read.csv("https://github.com/briatte/ggnet/raw/master/inst/extdata/network.tsv", sep = "\t")
# # network object
# net = network(df)
#
# # party affiliation
# x = data.frame(Twitter = network.vertex.names(net))
# x = merge(x, ids, by = "Twitter", sort = FALSE)$Groupe
# net %v% "party" = as.character(x)
#
# # colour palette
# library(RColorBrewer)
# y = brewer.pal(9, "Set1")[ c(3, 1, 9, 6, 8, 5, 2) ]
# names(y) = levels(x)
#
# # network plot
# ggnet2(net, color = "party", palette = y, alpha = 0.75, size = 4, edge.alpha = 0.5)
#
#
#
# ids_2 <- nodes %>% select(-members) %>% as.data.frame()
# df_2 <- links %>% select(-jaccard) %>% as.data.frame()
# net_2 <- network(df_2)
#
# network.vertex.names(net_2)
# subs <- data.frame(geneset = network.vertex.names(net_2), stringsAsFactors = F) %>%
#   inner_join(ids_2, by = 'geneset') %>%
#   nth(3)
#
# net_2 %v% 'subcollection' <- subs
#
# colors <- colorRampPalette(RColorBrewer::brewer.pal(10, name = "Spectral"))(length(unique(subs)))
# names(colors) <- unique(subs)
#
# ggnet2(net_2, color = 'subcollection', palette = colors, alpha = 0.75, size = 4, edge.alpha = 0.5)

