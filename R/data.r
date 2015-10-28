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
#' Variables are as follows:
#'
#' @format A data frame with 10786 rows and 4 variables:
#' \itemize{
#'   \item collection: origin of the gene set (MSigDB: H_hallmark, C1_position -
#'   C7_immunologic collections, BTM or TISSUES)
#'   \item subcollection: collections are sometimes further subdivided (e.g.
#'   C2_canonical: CP_BIOCARTA, CP_KEGG, CP_PID, CP_REACTOME, etc.)
#'   \item geneset: name of the gene set in the collection
#'   \item members: gene members (using HUGO Gene Nomenclature Committee [HGNC]
#'   gene symbols)
#' }
"genesets"

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
#' Gene sets were mapped to equivalent validated targetting miRNAs using the
#' approach proposed in: Godard, P. & van Eyll, J. Pathway analysis from lists
#' of microRNAs: common pitfalls and alternative strategy. Nucl. Acids Res. 43,
#' 3490–3497 (2015).
#'
#' Variables are as follows:
#'
#' @format A data frame with 10786 rows and 4 variables:
#' \itemize{
#'   \item collection: origin of the gene set (MSigDB: H_hallmark, C1_position -
#'   C7_immunologic collections, BTM or TISSUES)
#'   \item subcollection: collections are sometimes further subdivided (e.g.
#'   C2_canonical: CP_BIOCARTA, CP_KEGG, CP_PID, CP_REACTOME, etc.)
#'   \item geneset: name of the gene set in the collection
#'   \item members: miRNA members (using naming convention used by mirBase; described:
#'   Victor Ambros, Bonnie Bartel, David P. Bartel, Christopher B. Burge, James
#'   C. Carrington, Xuemei Chen, Gideon Dreyfuss, Sean R. Eddy, Sam
#'   Griffiths-Jones, Mhairi Marshall, Marjori Matzke, Gary Ruvkun, and Thomas
#'   Tuschl. A uniform system for microRNA annotation. RNA 2003 9(3):277-279.)
#' }
"genesets_mirs"

#' Adjacency of 10,786 annotated gene sets as a dyadic table
#'
#' A dataset containing the Jaccard coefficient for all pairwise comparisons
#' between all genesets. Self reflective edges and zero edges are removed. The
#' variables are as follows:
#'
#' @format A data frame with 22,081,622 rows and 3 variables:
#' \itemize{
#'   \item source: ...
#'   \item target: ...
#'   \item jaccard: ...

#' }
"genesets_links"
