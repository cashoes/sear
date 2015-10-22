#' Code for parsing MSigDB GMT files obtained from the Broad Institute

# creat gene symbol version of genesets ----------------------------------------
files <- list.files(path = "data-raw/misgdb/", pattern = ".gmt$", full.names = T)
names(files) <- c("BTM_ALL",
                  "C1_positional_ALL",
                  "C2_canonical_CGP",
                  "C2_canonical_CP_BIOCARTA",
                  "C2_canonical_CP_KEGG",
                  "C2_canonical_CP_REACTOME",
                  "C2_canonical_CP_PID",
                  "C3_motif_MIR",
                  "C3_motif_TFT",
                  "C4_computational_CGN",
                  "C4_computational_CN",
                  "C5_ontology_BP",
                  "C5_ontology_CC",
                  "C5_ontology_MF",
                  "C6_oncogenic_ALL",
                  "C7_immunologic_ALL",
                  "H_hallmark_ALL")

# read in as list of lists
tbl_genesets <- lapply(files, function(gmt){
  geneset <- readLines(gmt)
  geneset <- sapply(geneset, function(x) strsplit(x, split = "\t"))
  names(geneset) <- lapply(geneset, function(x) head(x, 1))
  geneset <- lapply(geneset, function(x) tail(x, -2))
})

# clean up C2 collection - CP includes all other CP_x + CP_PID
gsub("^(BIOCARTA|KEGG|REACTOME|PID).+", "\\1", names(tbl_genesets$C2_canonical_CP_PID)) %>%
  table()

# check to make sure all all included in respective sub-collections
all(grep("KEGG", names(tbl_genesets$C2_canonical_CP_PID), value = T) %in% names(tbl_genesets$C2_canonical_CP_KEGG))
all(grep("BIOCARTA", names(tbl_genesets$C2_canonical_CP_PID), value = T) %in% names(tbl_genesets$C2_canonical_CP_BIOCARTA))
all(grep("REACTOME", names(tbl_genesets$C2_canonical_CP_PID), value = T) %in% names(tbl_genesets$C2_canonical_CP_REACTOME))

# move others to new sub-collection
tbl_genesets$C2_canonical_CP_OTHER <- tbl_genesets$C2_canonical_CP_PID[grep("BIOCARTA|KEGG|REACTOME|PID", names(tbl_genesets$C2_canonical_CP_PID), invert = T)]

# keep only PID in PID sub-collection
tbl_genesets$C2_canonical_CP_PID <- tbl_genesets$C2_canonical_CP_PID[grep("PID", names(tbl_genesets$C2_canonical_CP_PID))]

# sort list by names
tbl_genesets <- tbl_genesets[order(names(tbl_genesets))]

# convert to list of data.frames
tbl_genesets <- lapply(tbl_genesets, function(x) data.frame(geneset = names(x), members = I(x)))

# convert to data.frame
tbl_genesets <- dplyr::bind_rows(tbl_genesets, .id = "collection")

# add sub-collections
tbl_genesets$subcollection <- gsub("BTM_(ALL)|C[0-9]_[a-z]+_(.+)|H_[a-z]+_(ALL)", "\\1\\2\\3", tbl_genesets$collection)
tbl_genesets$collection <- gsub("(BTM)_ALL|(C[0-9]_[a-z]+)_.+|(H_[a-z]+)_ALL", "\\1\\2\\3", tbl_genesets$collection)

# reorder
tbl_genesets <- dplyr::select(tbl_genesets, collection, subcollection, geneset, members)

# add tissues
tissues <- readRDS("data-raw/Tissues_tbl_df.rds")
tissues_map <- readRDS("data-raw/Tissues_map_tbl_df.rds")
all(tissues_map$celltype == tissues$geneset)
tissues$subcollection <- tissues_map$category
tissues <- dplyr::select(tissues, collection, subcollection, geneset, members)

# final object
genesets <- rbind(tbl_genesets, tissues)

# cleanup ' /// ' separators across all genesets
genesets$members <- map(genesets$members, function(x) unlist(strsplit(x, split = ' /// ')))
