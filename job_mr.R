library(snpnet)
library(glmnetPlus)
library(reshape2)
library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

# provide 3 arguments: desired rank, prev_iter (if none, set 0), weighted (set 0 for now)

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

source("multnet.R")

#####----------- Configs -----------#####

genotype_file <- "/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/biobank-methods-dev/private_data/data-split/train.bed"   # training bed
genotype_file_val <- "/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/biobank-methods-dev/private_data/data-split/val.bed"   # validation bed

genotype_p2file <- "/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/biobank-methods-dev/private_data/data-split/train"   # training bed
genotype_p2file_val <- "/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/biobank-methods-dev/private_data/data-split/val"   # validation bed

phenotype_file <- "/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/master_phe/master.20190509.phe"   # path to the phenotype file

phe_list <- fread("phe_long_3.csv")  # read list of phenotypes to analyze, two columns: (GBE_ID, phenotype)
results_dir <- "/scratch/users/junyangq/multiresponse/results_mr_long_3_renew/"  # parent results directory to save intermediate results
covariate_names <- c("age", "sex", paste0("PC", 1:10))


standardize_response <- TRUE  # should we standardize response beforehand (rescale when prediction for sure)
save <- TRUE  # should the program save intermediate results?

nlambda <- 100
lambda.min.ratio <- 0.01
batch_size <- 1000  # size of candidate variant batch in the screening
max.iter <- 50  # maximum BAY iteration
thresh <- 1e-7  # convergence threshold
validation <- TRUE  # is validation set provided
early_stopping <- TRUE  # should we adopt early stopping

# other computational configurations
configs <- list(
  missing.rate = 0.1,  # variants above this missing rate are discarded
  MAF.thresh = 0.001,  # MAF threshold
  nCores = 8,  # number of cores to be used
  bufferSize = 10000,  # number of COLUMNS (Variants) the memory can hold at a time
  standardize.variant = FALSE,  # standardize predictors or not
  results.dir = paste0("results_rank_", rank, "/"),  # subdirectory for each rank
  meta.dir = "meta/",
  use_plink2 = use_plink2,
  mem = 128000  # memeory guidance for PLINK2
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

