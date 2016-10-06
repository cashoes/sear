#' creat mirna version of genesets
#' @author C.Shannon
#'
#' using mirWalk2.0 data:
#' http://www.umm.uni-heidelberg.de/apps/zmf/mirwalk/holistic.html

# read in as list of lists
files <- list.files(path = "data-raw/mirwalk2/", pattern = "*.gmt$", full.names = T)
names(files) <- c('validated', '3UTR', '5UTR', 'Body', 'Promoter')

# only validated
files <- files[1]

# read files in
tbl_mirs <- lapply(files, function(gmt){
  gmt <- readLines(gmt)
  gmt <- sapply(gmt, function(x) strsplit(x, split = "\t"))
  names(gmt) <- lapply(gmt, function(x) head(x, 1))
  gmt <- lapply(gmt, function(x) tail(x, -3))
})
# convert lists to list-columns and bind into data.frame
tbl_mirs <- lapply(tbl_mirs, function(x) dplyr::data_frame(geneset = names(x), members = I(x))) %>%
  dplyr::bind_rows(.id = "subcollection") %>%
  dplyr::mutate(collection = 'mirwalk') %>%
  dplyr::select(collection, subcollection, geneset, members) %>%
  dplyr::tbl_df()

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

system.time({
  collections <- collections %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(members_mirna = list(getmirs(members_mrna, mirs_annot)))
})

library(ggplot2)
collections %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(n = length(members_mirna)) %>%
  ggplot(aes(n, fill = collection)) +
  geom_histogram() +
  scale_x_log10() +
  facet_wrap(~collection, scales = 'free')

# clean up
rm(files, tbl_mirs, mirs_annot, btms, tissues, msigdb)

# finally, save object
collections <- collections %>% dplyr::ungroup()
devtools::use_data(collections, overwrite = T)
