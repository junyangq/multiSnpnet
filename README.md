# Fit a sparse reduced rank regression model on large-scale SNP data and multivariate responses

License: GPL-2

### Reference: 
  - Qian, Junyang, Wenfei Du, Yosuke Tanigawa, Matthew Aguirre, Robert Tibshirani, Manuel A. Rivas, and Trevor Hastie. "A Fast and Flexible Algorithm for Solving the Lasso in Large-scale and Ultrahigh-dimensional Problems." bioRxiv (2019): https://www.biorxiv.org/content/10.1101/630079v1
  - Ruilin Li, Christopher Chang, Johanne Marie Justesen, Yosuke Tanigawa, Junyang Qian, Trevor Hastie, Manuel A. Rivas, Robert Tibshirani. "Fast Lasso method for Large-scale and Ultrahigh-dimensional Cox Model with applications to UK Biobank." bioRxiv (2020): https://www.biorxiv.org/content/10.1101/2020.01.20.913194v1.full.pdf

### Installation:
Most of the requirements of snpnet are available from CRAN. It also depends on the `pgenlibr` and `glmnet/glmnetPlus` packages. One can install them by running the following commands in R. Notice that the installation of `pgenlibr` requires [zstd(>=1.4.4)](https://github.com/facebook/zstd). It can be built from source or simply available from [conda](https://anaconda.org/conda-forge/zstd), [pip](https://pypi.org/project/zstd/) or [brew](https://formulae.brew.sh/formula/zstd).

```r
library(devtools)
install_github("junyangq/glmnetPlus")
install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
```
We assume the users already have PLINK 2.0. Otherwise it can be installed from https://www.cog-genomics.org/plink/2.0/.

### Sample Usage:
library(multisnpnet)

genotype_file <- file.path(system.file("extdata", package = "multisnpnet"), "sample")
phenotype_file <- system.file("extdata", "sample.phe", package = "multisnpnet")

phenotype_names = c("QPHE", "BPHE")
covariate_names = c("age", "sex", paste0("PC", 1:10))

out <- multisnpnet(
  genotype_file = genotype_file, 
  phenotype_file = phenotype_file, 
  phenotype_names = phenotype_names, 
  covariate_names = covariate_names, 
  rank = 2, 
  nlambda = 10, 
  validation = TRUE,
  split_col = "split",
  batch_size = 100, 
  standardize_response = TRUE, 
  save = TRUE
)
