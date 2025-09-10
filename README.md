scPS: power and sample size calculation in differential expression analysis of scRNAseq data
================
Chih-Yuan Hsu

Sept/10/2025

Chih-Yuan Hsu, Qi Liu, and Yu Shyr (2024). A distribution-free and analytic method for power and sample size calculation in single-cell differential expression. *Bioinformatics*. Volume 40, Issue 9, https://doi.org/10.1093/bioinformatics/btae540

## Installation

Download nifts_0.5.2.tar.gz and locally install it, or execute the following code:
``` r
library(devtools)
install_github("cyhsuTN/scPS")
```

## Usage

``` r
library(scPS)
```

### scPS Guidance for Experimental Design
#### Types of comparisons

- [Comparison between independent two groups](scPS_indep.md)
- [Comparison between paired groups](scPS_paired.md)

### Web Version of scPS
To simplify parameter settings, the web version supports power and sample size calculation for only one specific cell type of interest. However, the total number of required cells per subject can be estimated by dividing the number of required cells of interest by their proportion in the cell population. For example, if 100 cells of interest are needed and their proportion is 20%, the total required cells would be 500 (=100/0.2).

- Independent two-group comparison:
  <https://cyhsutn.shinyapps.io/scPS_shiny_Indep/>
- Paired-group comparison:
  <https://cyhsutn.shinyapps.io/scPS_shiny_Paired/>




