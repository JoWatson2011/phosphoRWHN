# RWHN_Phosphoproteomics
A method to infer the function of individual phosphorylation sites based on phosphoproteomics data, using random walks on heterogeneous networks (RWHN).

# To install package

Run the following in an R session:
`devtools::install_github("JoWatson2011/phosphoRWHN",  build_vignettes = T)`

The [Bioconductor](https://www.Bioconductor.org) packages `org.Hs.eg.db`, `GOSemSim`, `GO.db`, `AnnotationDbi` are also required and may need installing separately.

# To view vignette
Run the following in an R session:
```
library(phosphoRWHN)
vignette("phosphoRWHN_MAPKERK")
```


![Figure 1](https://raw.githubusercontent.com/JoWatson2011/RWHN_Phosphoproteomics/master/figure/Figure1.png?raw=true "Title")
 
Manuscript can be viewed [here](https://doi.org/10.1021/acs.jproteome.1c00150):
Watson, Joanne, Jean-Marc Schwartz, and Chiara Francavilla. 
“Using Multilayer Heterogeneous Networks to Infer Functions of Phosphorylated Sites.” 
Journal of Proteome Research 20, no. 7 (July 2, 2021): 3532–48. https://doi.org/10.1021/acs.jproteome.1c00150.

