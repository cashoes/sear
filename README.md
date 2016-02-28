# sear
##Simple (effin') Enrichment Analysis in R

A package for simple pathway over-representation analysis.

```
install.packages('devtools')
library(devtools)
install_github('cashoes/sear')

library(sear)
candidates <- c('ACTB', 'B2M', 'SERPINC1', 'UBC')
sear(candidates)
```

Annotated genesets are aggregated from the following ressources:

  1.MSigDB: Subramanian, Aravind, et al. "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." Proceedings of the National Academy of Sciences of the United States of America 102.43 (2005): 15545-15550.
  
  2.Blood Transciptome Modules (BTM): Li, Shuzhao, et al. "Molecular signatures of antibody responses derived from a systems biology study of five human vaccines." Nature immunology 15.2 (2014): 195-204
  
  3. Chaussabel D, Baldwin N. Democratizing systems immunology with modular transcriptional repertoire analyses. Nat Rev Immunol. 2014;14: 271–280. doi:10.1038/nri3642
  
  4.Tissue Enrichment: Benita, Yair, et al. "Gene enrichment profiles reveal T-cell development, differentiation, and lineage-specific transcription factors including ZBTB25 as a novel NF-AT repressor." Blood 115.26 (2010): 5376-5384.
  
In addition, gensets are mapped to facilitate analysis of miRNA candidate lists using the strategy described in:

  1. Godard P, van Eyll J. Pathway analysis from lists of microRNAs: common pitfalls and alternative strategy. Nucl Acids Res. 2015;43: 3490–3497. doi:10.1093/nar/gkv249
