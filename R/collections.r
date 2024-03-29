#' Collection of 34,010 annotated gene sets
#'
#' A dataset containing >34,000 annotated gene sets from the following sources:
#'
#' \itemize{
#'   \item MSigDB (33,591): Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America. 2005;102: 15545-15550. doi:10.1073/pnas.0506580102
#'   \item Blood Transciptome Modules (BTM; 346): Chaussabel D, Baldwin N. Democratizing systems immunology with modular transcriptional repertoire analyses. Nat Rev Immunol. 2014;14: 271-280. doi:10.1038/nri3642
#'   \item Tissue Enrichment (73): Benita Y, Cao Z, Giallourakis C, Li C, Gardet A, Xavier RJ. Gene enrichment profiles reveal T-cell development, differentiation, and lineage-specific transcription factors including ZBTB25 as a novel NF-AT repressor. Blood. 2010;115: 5376-5384. doi:10.1182/blood-2010-01-263855
#' }
#'
#' Gene set members were also mapped to equivalent validated targetting miRNAs
#' to enable unbiased enrichment analysis of miRNA lists, as described in:
#'
#' \itemize{
#'   \item Godard, P. & van Eyll, J. Pathway analysis from lists of microRNAs: common pitfalls and alternative strategy. Nucl. Acids Res. 43, 3490-3497 (2015).
#' }
#'
#' @format A data frame with 34,010 rows and 4 variables:
#' \itemize{
#'   \item collection: origin of the gene set (MSigDB: H - hallmark,
#'   C1 - positional, C2 - curated, C3 - motif, C4 - computational,
#'   C5 - GO gene ontology, C6 - oncogenic, C7 - immunologic; BTM or TISSUES)
#'   \item subcollection: collections are sometimes further subdivided (e.g.
#'   C2: CP:BIOCARTA, CP:KEGG, CP:REACTOME, etc.)
#'   \item geneset: name of the gene set in the collection
#'   \item description: a short description of the gene set (only available for
#'   some MSigDB collections)
#'   \item members_mrna: gene members (using HUGO Gene Nomenclature Committee
#'   [HGNC] gene symbols)
#'   \item members_mirna: miRNA members (using naming convention used by
#'   mirBase; described in:
#'   Victor Ambros, Bonnie Bartel, et al. A uniform system for microRNA
#'   annotation. RNA 2003 9(3):277-279.)
#' }
"collections"
