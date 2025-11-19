![Scripting](https://img.shields.io/badge/Language-R-yellow.svg) ![Copyright](https://img.shields.io/badge/Copyright-(c)_2025_Max\_Stammnitz\_@CRG\_Barcelona-green.svg)

## Analysis scripts accompanying Stammnitz and Lehner, 2025

_Lehner Laboratory (Genetic Systems), Centre for Genomic Regulation (2022-2025)_

This repository contains custom R scripts which - in conjunction with the associated supplementary files - can be used to replicate the main and supplementary figures presented in: 

**[The genetic architecture of an allosteric hormone receptor (bioRxiv, 2025)](https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1)**

![example](/aux/cover.png)

The scripts are written to run on **[R](https://www.r-project.org/)** version 4.4.1 or later. You should be able to check your current version of R by running the command below:

```
R --version
```

The current R version is also shown when opening RStudio or the R Console.

Scripts require the following R packages: [**`readxl`**](https://cran.r-project.org/web/packages/readxl/index.html), [**`stringr`**](https://cran.r-project.org/web/packages/stringr/index.html), [**`scales`**](https://cran.r-project.org/web/packages/scales/index.html), [**`bio3d`**](https://cran.r-project.org/web/packages/bio3d/index.html), [**`drc`**](https://cran.r-project.org/web/packages/drc/index.html), [**`growthrates`**](https://cran.r-project.org/web/packages/growthrates/index.html), [**`reshape`**](https://cran.r-project.org/web/packages/reshape/index.html), [**`ggplot2`**](https://cran.r-project.org/web/packages/ggplot2/index.html), [**`ggtext`**](https://cran.r-project.org/web/packages/ggtext/index.html), [**`ggrepel`**](https://cran.r-project.org/web/packages/ggrepel/index.html), [**`GGally`**](https://cran.r-project.org/web/packages/GGally/index.html), [**`cowplot`**](https://cran.r-project.org/web/packages/cowplot/index.html), [**`pheatmap`**](https://cran.r-project.org/web/packages/pheatmap/index.html), [**`beeswarm`**](https://cran.r-project.org/web/packages/beeswarm/index.html), [**`viridis`**](https://cran.r-project.org/web/packages/viridis/index.html), [**`wesanderson`**](https://cran.r-project.org/web/packages/wesanderson/index.html).

Although care has been taken to make the code distribution-independent, it is possible that some of the scripts only work on Unix/MacOS systems, and may need to be modified in order to run on Windows systems.

Sequencing data for reproducing the raw [variant fitness estimates](/data/DiMSum) can be found under the European Nucleotide Archive BioProject accession no. [PRJEB89674](https://www.ebi.ac.uk/ena/browser/view/PRJEB89674).

For any feedback or requests, please get in touch with me via email: maximilian.stammnitz@crg.eu
