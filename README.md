# Fast Multi-Phenotype SRRR on Genetic Data

Fit a sparse reduced rank regression model on large-scale genetic data and multivariate responses with batch variable screening and alternating minimization. It computes a full solution path on a grid of penalty values. Can deal with larger-than-memory SNP data, missing values, and confounding covariates.

License: GPL-2

### Installation:
Most of the requirements of snpnet are available from CRAN. It also depends on `snpnet`, `pgenlibr` and `glmnet/glmnetPlus` packages. One can install them by running the following commands in R. Notice that the installation of `pgenlibr` requires [zstd(>=1.4.4)](https://github.com/facebook/zstd). It can be built from source or simply available from [conda](https://anaconda.org/conda-forge/zstd), [pip](https://pypi.org/project/zstd/) or [brew](https://formulae.brew.sh/formula/zstd).

```r
library(devtools)
install_github("junyangq/snpnet")
install_github("junyangq/glmnetPlus")
install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
```
We assume the users already have PLINK 2.0. Otherwise it can be installed from https://www.cog-genomics.org/plink/2.0/.

### Sample Usage:

For description of the arguments, see the function documentation (e.g. by `?multisnpnet`).

```r
library(multiSnpnet)

genotype_file <- file.path(system.file("extdata", package = "multiSnpnet"), "sample")
phenotype_file <- system.file("extdata", "sample.phe", package = "multiSnpnet")

phenotype_names <- c("QPHE", "BPHE")
covariate_names <- c("age", "sex", paste0("PC", 1:10))

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
```

### References:

[1] Junyang Qian, Yosuke Tanigawa, Ruilin Li, Robert Tibshirani, Manuel A. Rivas, and Trevor Hastie. "Large-Scale Sparse Regression for Multiple Responses with Applications to UK Biobank." bioRxiv, 2020. doi: 10.1101/2020.05.30.125252

[2] Lisha Chen, and Jianhua Huang. "Sparse reduced-rank regression for simultaneous dimension reduction and variable selection." Journal of the American Statistical Association 107.500 (2012): 1533-1545.

