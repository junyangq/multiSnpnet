# "Rscript ../job_mr.R ${rank[$idx]} -1 ${weighted[$idx]} ${use_plink2[$idx]} ${mmem[$idx]} 16 2000 ../phe_long_721.csv
library(devtools)
#devtools::load_all('/oak/stanford/groups/mrivas/software/snpnet/snpnet_v.0.3.4/')
#snpnet must be from install.packages("remotes")
#remotes::install_github("junyangq/snpnet")
library(snpnet)
library(BEDMatrixPlus)
library(glmnetPlus)
library(reshape2)
library(ggplot2)
library(data.table)
library(magrittr)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

# provide 8 arguments: desired rank, prev_iter (if none, set 0), weighted (set 0 for now), use_plink2 (TRUE, FALSE), memory, ncores, batchsize, pheno_list

if (length(args) < 3) {
  stop("At least two arguments must be supplied (rank, prev_iter, weighted).", call.=FALSE)
} else {
  # default output file
  rank <- as.integer(args[1])
  prev_iter <- as.integer(args[2])
  weighted <- as.integer(args[3])
  if (length(args) > 3) {
    use_plink2 <- as.integer(args[4])
  } else {
    use_plink2 <- FALSE
  }
}

source("../multnet.R")

#####----------- Configs -----------#####

genotype_file <- "/scratch/groups/mrivas/ukbb24983/array_combined/pgen/split/train.bed"   # training bed
genotype_file_val <- "/scratch/groups/mrivas/ukbb24983/array_combined/pgen/split/val.bed"   # validation bed

genotype_p2file <- "/scratch/groups/mrivas/ukbb24983/array_combined/pgen/split/train"   # training bed
genotype_p2file_val <- "//scratch/groups/mrivas/ukbb24983/array_combined/pgen/split/val"   # validation bed

phe_list_in <- args[8]
if(phe_list_in == "../phe_long_721.csv"){
phenotype_file <- "/oak/stanford/groups/mrivas/projects/biomarkers/snpnet/biomarkers/biomarkers_covar.phe"   # path to the phenotype file
} else{
phenotype_file <- "/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/master_phe/master.20190509.phe"
}

phe_list <- fread(phe_list_in, sep = ",")  # read list of phenotypes to analyze, two columns: (GBE_ID, phenotype)
results_dir <- paste0(getwd(), "/results_rank_",rank)  # parent results directory to save intermediate results
# covariate_names <- c("age", "sex", paste0("PC", 1:10))
covariate_names <- c()

standardize_response <- TRUE  # should we standardize response beforehand (rescale when prediction for sure)
save <- TRUE  # should the program save intermediate results?

nlambda <- 100
lambda.min.ratio <- 0.01
batch_size <- as.integer(args[7])  # size of candidate variant batch in the screening
max.iter <- 50  # maximum BAY iteration
thresh <- 1e-7  # convergence threshold
validation <- TRUE  # is validation set provided
early_stopping <- TRUE  # should we adopt early stopping

# other computational configurations
configs <- list(
  missing.rate = 0.2,  # variants above this missing rate are discarded
  MAF.thresh = 0.001,  # MAF threshold
  nCores = as.integer(args[6]),  # number of cores to be used
  bufferSize = 50000,  # number of COLUMNS (Variants) the memory can hold at a time
  standardize.variant = FALSE,  # standardize predictors or not
  results.dir = paste0("results_rank_", rank, "/"),  # subdirectory for each rank
  meta.dir = "meta/",
  use_plink2 = use_plink2,
  mem = as.integer(args[5]),  # memeory guidance for PLINK2
  KKT.verbose = FALSE
)

######------------------------------########

r <- rank
phenotype_names <- phe_list[["GBE_ID"]]
weight <- NULL

fall <- SRRR_path(genotype_file = genotype_file, 
                  phenotype_file = phenotype_file, 
                  phenotype_names = phenotype_names, 
                  covariate_names = covariate_names, 
                  results_dir = results_dir, 
                  r = r, 
                  batch_size = batch_size, 
                  max.iter = max.iter, 
                  thresh = thresh, 
                  configs = configs,
                  lambda.min.ratio = lambda.min.ratio,
                  nlambda = nlambda,
                  standardize_response = standardize_response,
                  save = save, 
                  validation = validation, 
                  genotype_file_val = genotype_file_val,
                  prev_iter = prev_iter,
                  early_stopping = early_stopping,
                  glmnet_thresh = thresh,
                  weight = weight,
                  use_plink2 = use_plink2,
                  genotype_p2file = genotype_p2file,
                  genotype_p2file_val = genotype_p2file_val)

