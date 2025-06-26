scPS: power and sample size calculation in differential expression analysis of scRNAseq data
================
Chih-Yuan Hsu

June/26/2025

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

- Independent two-group comparison:
  <https://cyhsutn.shinyapps.io/scPS_shiny_Indep/>
- Paired-group comparison:
  <https://cyhsutn.shinyapps.io/scPS_shiny_Paired/>

