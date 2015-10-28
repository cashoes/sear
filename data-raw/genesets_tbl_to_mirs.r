# creat mirna version of genesets ----------------------------------------------
# using mirWalk2.0 data:
# http://www.umm.uni-heidelberg.de/apps/zmf/mirwalk/holistic.html

# read in as list of lists
files <- list.files(path = "data-raw/mirwalk_2/", pattern = "*.gmt$", full.names = T)
names(files) <- c('validated', '3UTR', '5UTR', 'Body', 'Promoter')

# read files in
tbl_mirs <- lapply(files, function(gmt){
  gmt <- readLines(gmt)
  gmt <- sapply(gmt, function(x) strsplit(x, split = "\t"))
  names(gmt) <- lapply(gmt, function(x) head(x, 1))
  gmt <- lapply(gmt, function(x) tail(x, -3))
})
# convert lists to list-columns and bind into data.frame
tbl_mirs <- lapply(tbl_mirs, function(x) data.frame(geneset = names(x), members = I(x))) %>%
  dplyr::bind_rows(.id = "subcollection") %>%
  dplyr::mutate(collection = 'mirwalk') %>%
  dplyr::select(collection, subcollection, geneset, members)

getmirs <- function(geneset, mirs_annot) {
  mirs_annot %>%
    dplyr::last() %>%
    purrr::map_int(function(x) length(intersect(x, geneset))) %>%
    `>`(0) %>%
    which() %>%
    dplyr::slice(mirs_annot, .) %>%
    dplyr::nth(3)
}

# system.time({
#   mirs_annot %>%
#     dplyr::last() %>%
#     purrr::map_int(function(x) length(intersect(x, geneset))) %>%
#     `>`(0) %>%
#     which() %>%
#     dplyr::slice(mirs_annot, .) %>%
#     dplyr::nth(3)
# })
#
# # this is slower and also returns partial matches... don't use
# system.time({
#   mirs_annot$geneset[grep(paste(geneset, collapse = "|"), mirs_annot$members)] %>%
#     strsplit(",") %>%
#     unlist() %>%
#     unique()
# })


# caution, slow ~6min with 'validated' subcollection only
mirs_annot <- tbl_mirs %>% filter(subcollection == 'validated')
genesets_mirs <- genesets %>%
  dplyr::rowwise() %>%
  dplyr::mutate(members = I(list(getmirs(members, mirs_annot))))

save(genesets_mirs, file = 'data/genesets_mirs.rda')

# system.time({
#   foo <- genesets %>%
#   slice(1:100) %>%
#   last() %>%
#   map(getmirs, mirs_annot)
# })
#
# system.time({
#   genesets %>%
#     dplyr::slice(1:100) %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(members = I(list(getmirs(members, mirs_annot))))
# })