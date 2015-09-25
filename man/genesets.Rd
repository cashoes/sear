\name{genesets}
\alias{genesets}
\docType{data}
\title{
Annotated gene sets from various sources.
}
\description{
A dataset containing annotated gene sets from the Broad Institute's MSigDB collection (1), Blood Transcriptome Modules (BTM) gene sets from (2) and tissue enrichment gene sets from (3).
}
\usage{data("genesets")}
\format{
  A data frame with 10767 observations on the following 3 variables.
  \describe{
    \item{\code{collection}}{name of the collection from which the gene set was obtained}
    \item{\code{geneset}}{name of the gene set - usually a functional annotation or provenance}
    \item{\code{members}}{a list - the genes that make up the gene set}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
\url{http://www.broadinstitute.org/gsea/msigdb/index.jsp}
}
\references{
1. Subramanian, Aravind, et al. "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." Proceedings of the National Academy of Sciences of the United States of America 102.43 (2005): 15545-15550.

2. Li, Shuzhao, et al. "Molecular signatures of antibody responses derived from a systems biology study of five human vaccines." Nature immunology 15.2 (2014): 195-204

3. Benita, Yair, et al. "Gene enrichment profiles reveal T-cell development, differentiation, and lineage-specific transcription factors including ZBTB25 as a novel NF-AT repressor." Blood 115.26 (2010): 5376-5384.
}
\examples{
data(genesets)
}
\keyword{datasets}
