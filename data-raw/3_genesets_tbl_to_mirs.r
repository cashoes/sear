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
tbl_mirs <- lapply(tbl_mirs, function(x) dplyr::tibble(geneset = names(x), members = x)) %>%
  dplyr::bind_rows(.id = "subcollection") %>%
  dplyr::mutate(collection = 'mirwalk') %>%
  dplyr::select(collection, subcollection, geneset, members) %>%
  dplyr::as_tibble()

# caution, slow ~6min with 'validated' subcollection only
mirs_annot <- tbl_mirs %>% filter(subcollection == 'validated')

getmirs <- function(geneset, mirs_annot) {
  mirs_annot %>%
    tidyr::unnest(members) %>%
    dplyr::filter(members %in% !!geneset) %>%
    dplyr::distinct(geneset) %>%
    purrr::flatten_chr()
}

system.time({
  library(future)
  plan(multisession, workers = 10)
  collections$members_mirna <- collections$members_mrna %>%
    furrr::future_map(
      .options = furrr::furrr_options(globals = c('mirs_annot', 'getmirs')),
      .progress = T,
      ~ getmirs(.x, mirs_annot)
    )
  plan(sequential)
})

library(ggplot2)
collections %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(n = length(members_mirna)) %>%
  ggplot2::ggplot(aes(n, fill = collection)) +
  ggplot2::geom_histogram() +
  ggplot2::scale_x_log10() +
  ggplot2::facet_wrap(~collection, scales = 'free')

# clean up
rm(files, tbl_mirs, mirs_annot, btms, tissues, msigdb)

# finally, save object
usethis::use_data(collections, overwrite = T)
