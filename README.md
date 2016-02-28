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

1. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America. 2005;102: 15545 –15550. doi:10.1073/pnas.0506580102
2. Benita Y, Cao Z, Giallourakis C, Li C, Gardet A, Xavier RJ. Gene enrichment profiles reveal T-cell development, differentiation, and lineage-specific transcription factors including ZBTB25 as a novel NF-AT repressor. Blood. 2010;115: 5376–5384. doi:10.1182/blood-2010-01-263855
3. Li S, Rouphael N, Duraisingham S, Romero-Steiner S, Presnell S, Davis C, et al. Molecular signatures of antibody responses derived from a systems biology study of five human vaccines. Nature Immunology. 2013;15: 195–204. doi:10.1038/ni.2789
4. Chaussabel D, Baldwin N. Democratizing systems immunology with modular transcriptional repertoire analyses. Nat Rev Immunol. 2014;14: 271–280. doi:10.1038/nri3642
  
In addition, genesets are mapped to facilitate analysis of miRNA candidate lists using the strategy described in:

1. Godard P, van Eyll J. Pathway analysis from lists of microRNAs: common pitfalls and alternative strategy. Nucl Acids Res. 2015;43: 3490–3497. doi:10.1093/nar/gkv249
