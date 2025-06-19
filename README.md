# IsoPairFinder <img src="man/figures/logo_250612.png" align="right" alt="[IsoPairFinder]" width="150" />


[![](https://www.r-pkg.org/badges/version/masstools?color=green)](https://cran.r-project.org/package=IsoPairFinder)

- Author: Zhiwei Zhou (zhouzw@stanford.edu)
- Created: 06/19/2025
- Last modified: 06/19/2025

## Overview
**IsoPairFinder** is an R package designed to identify potential intermediates from Stable Isotope Tracing (SIT) metabolomics data for accelerating the elucidation of the metabolism pathway in the gut microbes. It provides the first end-to-end workflow serving this objective, including (1) identifying differential ion signals introduced by the gene mutation; (2) merging of the redundant LC-MS signals (isotopes, adducts, and in-source fragments); (3) pairing 12C/13C features to determine potential intermediates. It is also compatible with different metabolomics data processing tools. 



## Installation
This workflow is based on R, which requires installing some dependent packages first. 

```r
# intall public packages
if (!require(devtools)){
    install.packages("devtools")
}

if (!require(BiocManager)){
    install.packages("BiocManager")
}

# Required packages
required_pkgs <- c("dplyr","tidyr","readr", "stringr", "tibble", "purrr",
"ggplot2", "igraph", "pbapply", "Rdisop", "randomForest", "pryr", "magrittr", "rmarkdown", "caret")
BiocManager::install(required_pkgs)
devtools::install_github("ZhuMetLab/SpectraTools")
devtools::install_github("tidymass/masstools")

# Install DoddLabTracer
devtools::install_github("DoddLab/DoddLabTracer")
```

## Usage
This is a basic example of HyuA mutants to get you started. The Study design could be found [here](xxx)

```r
library(tidyverse)
library(DoddLabTracer)

# analysis of HyuA 
find_intemidates(peak_table_unlabel = 'hyuA_UA_48h_area.txt',
                 peak_table_label = 'hyuA_13CUA_48h_area.txt',
                 path = '~/Project/00_Uric_Acid_project/Data/20240319_isotope_tracing_analysis/hyuA/',
                 control_group = c("WT_UA", "WT_13CUA"),
                 case_group = c('hyuA_UA', 'hyuA_13CUA'),
                 polarity = 'positive',
                 mz_tol = 10,
                 rt_tol = 0.05,
                 p_value_cutoff = 0.05,
                 fold_change_cutoff = 10,
                 p_adjust = FALSE,
                 is_recognize_adducts = TRUE)

```

For more detailed examples and vignettes, please refer to the package documentation:
```r
# View all vignettes
browseVignettes("IonPairFinder")

```


## Citation
If you use `IsoPairFinder` in your research, please consider citing it. 

- Zhiwei Zhou, IsoPairFinder: An R package for identifying potential intermediates from Stable Isotope Tracing metabolomics data, In preparing

You can also find citation information by running:
```r
citation("IsoPairFinder")
```


## Contact
For questions or issues, please open an issue on GitHub or contact Zhiwei Zhou (zhouzw@stanford.edu).


## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a> 
This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)


