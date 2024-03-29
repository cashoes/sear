% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sear.r
\name{sear}
\alias{sear}
\title{Simple enrichment analysis in R using the MSigDB collections.}
\usage{
sear(
  input,
  type = c("mrna", "mirna"),
  return_members = F,
  return_intersect = T
)
}
\arguments{
\item{input}{A character vector of HUGO gene naming commitee (HGNC) gene
symbols.}

\item{type}{A character string. The type of symbols in the input: 'mrna' for
gene symbols, 'mirna' for miRNAs. The reference gene sets used differ for
either types. See: Godard, P., and van Eyll, J. (2015). Pathway analysis
from lists of microRNAs: common pitfalls and alternative strategy. Nucl.
Acids Res. 43, 3490-3497.}

\item{return_members}{A logical. Whether to include a list-column in the
return object containing the gene set member genes.}

\item{return_intersect}{A logical. Whether to include a list-column in the
return object containing the genes found in common between the gene set and
the input.}
}
\description{
Carry out enrichment analysis against some or all of the MSigDB collections
(1), Blood Transcriptome genes (2) and Tissue Enrichment Sets (3).
}
\details{
A list of genes is compared to each annotated gene set in turn by performing
a hypergeometric test of the overlap. The size of the input gene list, gene
set, intersection and resulting p-value are returned. P-values can be
adjusted over each collection (optimistic) or globally (pessimistic), using
any of the methods available to `p/adjust()`.
}
\examples{

library(dplyr)
data('collections')
input <- collections$members_mrna \%>\% unlist() \%>\% unique() \%>\% sample(100)
output <- sear(input, type = 'mrna') \%>\%
  arrange(fdr) \%>\%
  slice(1:100)
}
