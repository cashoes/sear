#' Source code for parsing MSigDB GMT files obtained from the Broad Institute.
parse_msigdb_gmts <- function(path) {
  files <- list.files(path = path, pattern = ".gmt$", full.names = T)
  names(files) <- c("BTM", "C1_positional", "C2_canonical", "C3_computational",
                    "C4_computational", "C5_ontology", "C6_oncological",
                    "C7_immunological", "H_hallmark")

  tbl_genesets <- lapply(files, function(gmt){
    geneset <- readLines(gmt)
    geneset <- sapply(geneset, function(x) strsplit(x, split = "\t"))
    names(geneset) <- lapply(geneset, function(x) head(x, 1))
    geneset <- lapply(geneset, function(x) tail(x, -2))
  })

  tbl_genesets <- lapply(tbl_genesets, function(x) {
    data.frame(geneset = names(x), members = I(x))
  }) %>%
    bind_rows(.id = "collection")
}
