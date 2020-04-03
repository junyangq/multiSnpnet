library(data.table)
library(devtools)
library(snpnet)
library(magrittr)
library(reshape2)
library(ggplot2)
library(data.table)
library(magrittr)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)

source("/oak/stanford/groups/mrivas/users/mrivas/repos/multiresponse-ukbb/multnet.R")

#####----------- Configs -----------#####
full_rank <- 10
rank <- 2:10
cores <- 8
mem <- 120000
buffer_size <- 10000
use_plink2 <- T

folder_name <- file.path("/oak/stanford/groups/mrivas/users/mrivas/repos/multiresponse-ukbb/", paste0("results_rank_", full_rank))  # parent directory
results_dir <- paste0("/scratch/users/junyangq/multiresponse/results_mr_lazy_rank_", full_rank)
dir.create(results_dir)

genotype_file_train <- "/scratch/users/ytanigaw/tmp/snpnet/geno/array_combined/train.bed"
genotype_file_val <- "/scratch/users/ytanigaw/tmp/snpnet/geno/array_combined/val.bed"
genotype_p2file <- "/scratch/users/ytanigaw/tmp/snpnet/geno/array_combined/train"   # training bed
genotype_p2file_val <- "/scratch/users/ytanigaw/tmp/snpnet/geno/array_combined/val"   # validation bed
phenotype_file <- "/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/master_phe/master.20190509.phe"

prefix_result_file <- "output_lambda_"
suffix_result_file <- ".RData"
train_name <- "metric_train"
val_name <- "metric_val"
files_in_dir <- list.files(folder_name)
result_files <- files_in_dir[startsWith(files_in_dir, prefix_result_file)]
max_iter <- max(as.numeric(gsub(suffix_result_file, "", gsub(pattern = prefix_result_file, "", result_files))))
last <- max_iter
latest_result <- file.path(folder_name, paste0(prefix_result_file, max_iter, suffix_result_file))

e <- new.env()
load(latest_result, envir = e)

phenotype_names <- colnames(e[[train_name]])
binary_names <- colnames(e[["AUC_train"]])
covariate_names <- rownames(as.matrix(e[["W_init"]]))

save <- T
validation <- T
standardize_response <- T
# 
# if (weighted == 1) {
#   weight <- c(rep(1, 7), 7)
# } else {
#   weight <- rep(1, length(phenotype_names))
# }

weight <- rep(1, length(phenotype_names))

configs <- list(
  missing.rate = e$configs$missing.rate,
  MAF.thresh = e$configs$MAF.thresh,
  nCores = cores,
  bufferSize = buffer_size,
  standardize.variant = FALSE,
  meta.dir = "meta/",
  use_plink2 = use_plink2,
  results.dir = results_dir,
  mem = mem
)

######------------------------------########

configs <- setup_configs_directories(configs, covariate_names, save, results_dir, validation)

names(weight) <- phenotype_names
weight <- weight / sum(weight) * length(phenotype_names)

phe_master <- fread(phenotype_file, colClasses = c("FID" = "character", "IID" = "character"), select = c("FID", "IID", covariate_names, phenotype_names))
fill_missing(phe_master, phenotype_names, -9, NA) # replace -9 with NA
chr_train <- BEDMatrixPlus(genotype_file_train)
n_chr_train <- nrow(chr_train)
if (validation) {
  chr_val <- BEDMatrixPlus(genotype_file_val)
  n_chr_val <- nrow(chr_val)
}
q_train <- length(phenotype_names)
cat_ids <- paste(phe_master$FID, phe_master$IID, sep = "_")
ids_valid_phe <- cat_ids[apply(as.matrix(phe_master[, phenotype_names, with = F]), 1, function(x) any(!is.na(x)))]
ids_valid_gen <- rownames(chr_train)
ids_valid <- intersect(ids_valid_phe, ids_valid_gen)
rowIdx_subset_gen_train <- match(ids_valid, ids_valid_gen)
if (validation) {
  ids_valid_gen_val <- rownames(chr_val)
  ids_valid_val <- intersect(ids_valid_phe, ids_valid_gen_val)
  rowIdx_subset_gen_val <- match(ids_valid_val, ids_valid_gen_val)
}

phe_train <- phe_master[match(ids_valid, cat_ids), ]
if (validation) phe_val <- phe_master[match(ids_valid_val, cat_ids), ]

## summary statistics: missing rate, mean, standard deviation (if needed) ##
if (!use_plink2) {
  stats <- snpnet:::computeStats(chr_train, rowIdx_subset_gen_train, c("pnas", "means", "sds"),
                                 path = file.path(results_dir, "meta"), save = save, configs = configs, verbose = TRUE, buffer.verbose = TRUE)
} else {
  stats <- computeStats_P2(genotype_p2file, ids_valid, configs = configs)
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0('zstdcat ', paste0(genotype_p2file, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
}

response_train_0 <- as.matrix(phe_train[, phenotype_names, with = F])

if (standardize_response) {
  std_obj <- y_standardization(phe_train, phenotype_names, weight)
  phe_train <- std_obj$response
}

if (length(covariate_names) > 0) {
  covariates_train <- phe_train[, covariate_names, with = F]
  if (validation) covariates_val <- phe_val[, covariate_names, with = F]
} else {
  covariates_train <- covariates_val <- NULL
}

n_subset_train <- nrow(phe_train)
response_train <- as.matrix(phe_train[, phenotype_names, with = F])
missing_response_train <- is.na(response_train)
if (validation) {
  n_subset_val <- nrow(phe_val)
  response_val <- as.matrix(phe_val[, phenotype_names, with = F])
  missing_response_val <- is.na(response_val)
}

feature_names <- e$active

features_train <- snpnet:::prepareFeatures(chr_train, feature_names, stats, rowIdx_subset_gen_train)
if (validation) {
  features_val <- snpnet:::prepareFeatures(chr_val, feature_names, stats, rowIdx_subset_gen_val)
}

MSE_train_all <- array(NA, dim = c(length(rank), length(phenotype_names), last),
                       dimnames = list(paste0("rank-", rank), phenotype_names, 1:last))
MSE_val_all <- array(NA, dim = c(length(rank), length(phenotype_names), last),
                     dimnames = list(paste0("rank-", rank), phenotype_names, 1:last))

AUC_train_all <- array(NA, dim = c(length(rank), length(binary_names), last),
                       dimnames = list(paste0("rank-", rank), binary_names, 1:last))
AUC_val_all <- array(NA, dim = c(length(rank), length(binary_names), last),
                     dimnames = list(paste0("rank-", rank), binary_names, 1:last))

nactive <- matrix(NA, length(rank), last, dimnames = list(paste0("rank-", rank), 1:last))

var_train <- apply(response_train_0, 2, function(x) mean((x - mean(x, na.rm = T))^2, na.rm = T))
var_val <- apply(response_val, 2, function(x) mean((x - mean(x, na.rm = T))^2, na.rm = T))

if (!is.null(covariates_train)) covariates_train_matrix <- as.matrix(covariates_train)
if (!is.null(covariates_val)) covariates_val_matrix <- as.matrix(covariates_val)

for (j in 1:max_iter) {
  print(j)
  fname <- file.path(folder_name, paste0("output_lambda_", j, ".RData"))
  if (file.exists(fname)) {
    e <- new.env()
    load(fname, envir = e)
    activenames <- e$active
    activeC <- as.matrix(e$fit$C)[activenames, , drop = F]
    features_train_current <- features_train[, activenames, with = F]
    features_train_current_matrix <- as.matrix(features_train_current)
    
    features_val_current <- features_val[, activenames, with = F]
    features_val_current_matrix <- as.matrix(features_val_current)
    if (length(activenames) > 0) svdC <- svd(activeC)
    for (ir in 1:length(rank)) {
      r <- rank[ir]
      if (length(activenames) < r) {
        C_approx <- activeC
      } else {
        C_approx <- tcrossprod(svdC$u[, 1:r, drop = F] %*% diag(svdC$d[1:r], r), svdC$v[, 1:r, drop = F])
        colnames(C_approx) <- colnames(activeC)
      }
      pred_cov_train <- if (!is.null(covariates_train)) covariates_train_matrix %*% as.matrix(e$fit$W) else 0
      pred_cov_val <- if (!is.null(covariates_val)) covariates_val_matrix %*% as.matrix(e$fit$W) else 0
      pred_feat_train <- if (length(activenames) > 0) features_train_current_matrix %*% C_approx else matrix(0, nrow(features_train), ncol(response_train_0))
      pred_feat_val <- if (length(activenames) > 0) features_val_current_matrix %*% C_approx else matrix(0, nrow(features_val), ncol(response_val))
      pred_train <- sweep(pred_cov_train + pred_feat_train, 2, e$fit$a0, FUN = "+")
      pred_val <- sweep(pred_cov_val + pred_feat_val, 2, e$fit$a0, FUN = "+")
      if (standardize_response) {
        pred_train <- y_de_standardization(pred_train, std_obj$means, std_obj$sds, weight)
        pred_val <- y_de_standardization(pred_val, std_obj$means, std_obj$sds, weight)
      }
      MSE_train <- apply(pred_train - response_train_0, 2, function(x) mean((x-mean(x, na.rm=T))^2, na.rm=T))
      MSE_train_all[ir, , j] <- 1 - MSE_train / var_train
      
      for (ik in 1:length(binary_names)) {
        response_train_binary <- response_train_0[, binary_names[ik]] - 1
        response_val_binary <- response_val[, binary_names[ik]] - 1
        data_train_binary <- data.frame(response = response_train_binary, covariates_train_matrix, score = pred_train[, binary_names[ik]])
        fit_binary <- glm(response ~ ., family = binomial(), data = data_train_binary)
        prob_train <- predict(fit_binary, newdata = data_train_binary, type = "response")
        not_missing_train <- !is.na(response_train_binary)
        AUC_train_all[ir, ik, j] <- snpnet:::computeMetric(as.matrix(prob_train[not_missing_train], ncol = 1), response_train_binary[not_missing_train], "binomial")
        data_val_binary <- data.frame(covariates_val_matrix, score = pred_val[, binary_names[ik]])
        prob_val <- predict(fit_binary, newdata = data_val_binary, type = "response")
        not_missing_val <- !is.na(response_val_binary)
        AUC_val_all[ir, ik, j] <- snpnet:::computeMetric(as.matrix(prob_val[not_missing_val], ncol = 1), response_val_binary[not_missing_val], "binomial")
      }
      
      MSE_val <- apply(pred_val - response_val, 2, function(x) mean((x-mean(x, na.rm=T))^2, na.rm=T))
      R2_val <- 1 - MSE_val / var_val
      MSE_val_all[ir, , j] <- R2_val
      nactive[ir, j] <- length(e$current_active)
    }
    save(MSE_train_all, MSE_val_all, AUC_train_all, AUC_val_all, nactive, file = file.path(results_dir, "summary_results.RData"))
  } else {
    break
  }
}

for (ir in 1:length(rank)) {
  r <- rank[ir]
  metric_train <- t(MSE_train_all[ir, , ])
  metric_val <- t(MSE_val_all[ir, , ])
  AUC_train <- t(as.matrix(AUC_train_all[ir, , ], nrow = length(binary_names)))
  AUC_val <- t(as.matrix(AUC_val_all[ir, , ], nrow = length(binary_names)))
  dir.create(file.path(results_dir, paste0("results_rank_", rank[ir])))
  save(metric_train, metric_val, AUC_train, AUC_val, file = file.path(results_dir, paste0("results_rank_", rank[ir]), "output_lambda_0.RData"))
}
