#' 10,786 annotated gene sets
#'
#' A dataset containing nearly 11,000 annotated gene sets from the following
#' sources:
#'
#' \itemize{
#'   \item MSigDB: Subramanian, Aravind, et al. "Gene set enrichment
#' analysis: a knowledge-based approach for interpreting genome-wide expression
#' profiles." Proceedings of the National Academy of Sciences of the United
#' States of America 102.43 (2005): 15545-15550.
#'   \item Blood Transciptome
#' Modules (BTM): Li, Shuzhao, et al. "Molecular signatures of antibody
#' responses derived from a systems biology study of five human vaccines."
#' Nature immunology 15.2 (2014): 195-204
#'   \item Tissue Enrichment: Benita, Yair,
#' et al. "Gene enrichment profiles reveal T-cell development, differentiation,
#' and lineage-specific transcription factors including ZBTB25 as a novel NF-AT
#' repressor." Blood 115.26 (2010): 5376-5384.
#' }
#'
#' Gene sets were also mapped to equivalent validated targetting miRNAs to
#' enable unbiased enrichment analysis of miRNA lists, as described in:
#' Godard, P. & van Eyll, J. Pathway analysis from lists of microRNAs: common
#' pitfalls and alternative strategy. Nucl. Acids Res. 43, 3490–3497 (2015).
#'
#' Variables are as follows:
#'
#' @format A data frame with 10786 rows and 4 variables:
#' \itemize{
#'   \item collection: origin of the gene set (MSigDB: h_hallmark, c1_position -
#'   c7_immunologic collections, btm or tissues)
#'   \item subcollection: collections are sometimes further subdivided (e.g.
#'   c2_canonical: cp_biocarta, cp_kegg, cp_pid, cp_reactome, etc.)
#'   \item geneset: name of the gene set in the collection
#'   \item members_mrna: gene members (using HUGO Gene Nomenclature Committee
#'   [HGNC] gene symbols)
#'   #'   \item members_mirna: miRNA members (using naming convention used by
#'   mirBase; described in:
#'   Victor Ambros, Bonnie Bartel, et al. A uniform system for microRNA
#'   annotation. RNA 2003 9(3):277-279.)
#' }
"genesets"
