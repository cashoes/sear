#' Parse all MSigDB genesets from single XML release
#' @author C.Shannon
#'
#'   Vast improvement over previous version, as only requires single download
#'   and category, sub-category is obtained directly from file. Also added
#'   description from metadat for eacg geneset.

library(xml2)
library(dplyr)
library(purrr)

m <- read_xml('data-raw/msigdb/msigdb_v5.1.xml')
header <- xml_attrs(m)
sets <- xml_children(m)
# construct table
data.frame(collection    = xml_attr(sets, 'CATEGORY_CODE'),
           subcollection = xml_attr(sets, 'SUB_CATEGORY_CODE'),
           geneset       = xml_attr(sets, 'STANDARD_NAME'),
           description   = xml_attr(sets, 'DESCRIPTION_BRIEF'),
           members_mrna  = xml_attr(sets, 'MEMBERS_SYMBOLIZED'),
           stringsAsFactors = F) %>% tbl_df() -> msigdb

msigdb <- msigdb %>%
  mutate(members_mrna = map(members_mrna, function(x) strsplit(x, ',') %>% unlist()))

attributes(msigdb)$version <- header

# clean up
rm(m, header, sets)
