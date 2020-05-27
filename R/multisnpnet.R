#' Fast Multi-Phenotype SRRR on SNP Data
#'
#' Fit a sparse reduced rank regression model on large-scale SNP data and multivariate responses
#' with batch variable screening and alternating minimization. It computes a full solution path on a
#' grid of penalty values. Can deal with larger-than-memory SNP data, missing values and adjustment
#' covariates.
#'
#' @usage multisnpnet(genotype_file, phenotype_file, phenotype_names, binary_phenotypes = NULL,
#'   covariate_names, rank, nlambda = 100, lambda.min.ratio = 0.01, standardize_response = TRUE,
#'   weight = NULL, validation = FALSE, split_col = NULL, mem = NULL,
#'   batch_size = 100, prev_iter = 0, max.iter = 10, configs = NULL, save = TRUE,
#'   early_stopping = FALSE)
#'
#' @param genotype_file Path to the suite of genotype files. genotype_file.{pgen, psam, pvar.zst}
#'   must exist.
#' @param phenotype_file Path to the phenotype. The header must include FID, IID, covariate_names
#'   and phenotype_names.
#' @param binary_phenotypes Names of the binary phenotypes. AUC will be evaluated for binary
#'   phenotypes.
#' @param covariate_names Character vector of the names of the adjustment covariates.
#' @param rank Target rank of the model.
#' @param nlambda Number of penalty values.
#' @param lambda.min.ratio Ratio of the minimum penalty to the maximum penalty.
#' @param standardize_response Boolean. Whether to standardize the responses before fitting to deal
#'   with potential different units of the responses.
#' @param p.factor Named vector of separate penalty factors applied to each coefficient. This is a
#'   number that multiplies \code{lambda} to allow different shrinkage. Default is 1 for all
#'   variables. Can specify partially and the rest will be set to 1. Must be positive.
#' @param weight Numberic vector that specifies the (importance) weights for the responses.
#' @param validation Boolean. Whether to evaluate on validation set.
#' @param split_col Name of the column in the phenotype file that specifies whether each sample
#'   belongs to the training split or the validation split. The values are either "train" or "val".
#' @param mem Memory available for the program. It tells PLINK 2.0 the amount of memory it can
#'   harness for the computation. IMPORTANT if using a job scheduler.
#' @param batch_size Number of variants used in batch screening.
#' @param prev_iter Index of the iteration to start from (e.g. to resume a previously interrupted
#'   computation).
#' @param max.iter Maximum number of iterations allowed for alternating minimization.
#' @param configs List of additional configuration parameters. It can include:
#'                \describe{
#'                \item{nCores}{number of cores for the PLINK computation (default: 1)}
#'                \item{thresh}{convergence threshold for alternating minimization (default: 1E-7)}
#'                \item{glmnet.thresh}{convergence threshold for glmnet(Plus) (default: 1E-7)}
#'                \item{plink2.path}{path to the PLINK2.0 program, if not on the system path}
#'                \item{zstdcat.path}{path to the zstdcat program, if not on the system path}
#'                }
#' @param save Boolean. Whether to save intermediate results.
#' @param early_stopping. Whether to stop the process early if validation metric starts to fall.
#'
#' @importFrom data.table ':='
#' @importFrom magrittr '%>%'
#'
#' @export
multisnpnet <- function(genotype_file, phenotype_file, phenotype_names, binary_phenotypes = NULL, covariate_names,
                        rank, nlambda = 100, lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04), standardize_response = TRUE,
                        weight = NULL, p.factor = NULL, validation = FALSE, split_col = NULL, mem = NULL,
                        batch_size = 100, prev_iter = 0, max.iter = 10, configs = NULL, save = TRUE,
                        early_stopping = FALSE) {

  configs <- setupMultiConfigs(configs, genotype_file, phenotype_file, phenotype_names, covariate_names,
                               nlambda, mem, standardize_response, max.iter, rank, prev_iter, batch_size)

  start_all <- Sys.time()

  if (rank > length(phenotype_names)) stop("The specified rank (", rank, ") should not be greater than the number of responses (", length(phenotype_names), ").")

  cat("Start Sparse Reduced Rank Regression for ", paste(phenotype_names, collapse = ", "), ".\n", sep = "")

  ### ---- pre-process phenotypes and genotype data ---- ###
  ids <- list()
  ids[["psam"]] <- snpnet:::readIDsFromPsam(paste0(genotype_file, '.psam'))

  ctype <- c("FID" = "character", "IID" = "character")
  if (!is.null(split_col)) ctype[split_col] <- "character"
  phe_master <- data.table::fread(phenotype_file, colClasses = ctype, select = c("FID", "IID", split_col, covariate_names, phenotype_names))
  if (length(covariate_names) > 0) {  # remove individuals with missing covariate values
    cov_master <- as.matrix(phe_master[, covariate_names, with = F])
    cov_no_missing <- apply(cov_master, 1, function(x) all(!is.na(x)))
    phe_master <- phe_master[cov_no_missing, ]
  }
  phe_master[["ID"]] <- paste(phe_master[["FID"]], phe_master[["IID"]], sep = "_")
  phe_master <- phe_master %>%
    dplyr::left_join(data.frame(ID = ids[["psam"]], sort_order = seq_along(ids[["psam"]]), stringsAsFactors = FALSE), by = "ID") %>%
    dplyr::arrange(sort_order) %>% dplyr::select(-sort_order) %>%
    data.table::as.data.table()
  fill_missing(phe_master, phenotype_names, -9, NA) # replace -9 with NA

  if (is.null(split_col)) {
    splits <- c('train')
    ids[['train']] <- phe_master$ID
  } else {
    splits <- c('train', 'val')
    for (s in splits) {
      ids[[s]] <- phe_master$ID[phe_master[[split_col]] == s]
    }
  }

  q_train <- length(phenotype_names)
  # exclude samples that have missingness in all phenotypes, could change to some lower bound of phenotype missingness
  cat_ids <- phe_master[["ID"]]
  ids_valid_phe <- cat_ids[apply(as.matrix(phe_master[, phenotype_names, with = F]), 1, function(x) any(!is.na(x)))]
  ids_valid_gen <- ids[["train"]]
  ids_valid <- intersect(ids_valid_phe, ids_valid_gen)
  rowIdx_subset_gen <- match(ids_valid, ids_valid_gen)
  if (validation) {
    ids_valid_gen_val <- ids[["val"]]
    ids_valid_val <- intersect(ids_valid_phe, ids_valid_gen_val)
    rowIdx_subset_gen_val <- match(ids_valid_val, ids_valid_gen_val)
  }

  ## summary statistics: missing rate, mean, standard deviation (if needed) ##
  phe_train <- phe_master[match(ids_valid, cat_ids), ]
  if (validation) phe_val <- phe_master[match(ids_valid_val, cat_ids), ]

  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ', paste0(genotype_file, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  pvar <- pgenlibr::NewPvar(paste0(genotype_file, '.pvar.zst'))
  chr_train <- pgenlibr::NewPgen(paste0(genotype_file, '.pgen'), pvar = pvar, sample_subset = match(ids_valid, ids[['psam']]))
  n_chr_train <- length(ids[['train']])
  if (validation) {
    chr_val <- pgenlibr::NewPgen(paste0(genotype_file, '.pgen'), pvar = pvar, sample_subset = match(ids_valid_val, ids[['psam']]))
    n_chr_val <- length(ids[['val']])
  }
  pgenlibr::ClosePvar(pvar)

  stats <- snpnet:::computeStats(genotype_file, ids_valid, configs = configs)
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0('zstdcat ', paste0(genotype_file, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID

  if (is.null(p.factor)) {
    p.factor <- rep(1, length(vars))
    names(p.factor) <- vars
  } else {
    if (!all(vars %in% names(p.factor))) {
      warning("p.factor does not cover all variants. The missing penalties are set to 1.\n")
      p.factor[setdiff(vars, names(p.factor))] <- 1
    }
  }
  p.factor[covariate_names] <- 0

  if (is.null(binary_phenotypes)) {
    binary_phenotypes <- phenotype_names[apply(as.matrix(phe_train[, phenotype_names, with = F]), 2, function(x) all(unique(x[!is.na(x)]) %in% c(1, 2)))]
    print(binary_phenotypes)
  }

  if (length(binary_phenotypes) > 0) {
    phe_train[, binary_phenotypes] <- phe_train[, binary_phenotypes, with = F] - 1
    if (validation) phe_val[, binary_phenotypes] <- phe_val[, binary_phenotypes, with = F] - 1
  }

  response_train_0 <- as.matrix(phe_train[, phenotype_names, with = F])
  var_train <- apply(response_train_0, 2, function(x) mean((x - mean(x, na.rm = T))^2, na.rm = T))

  if (save) {
    dir.create(file.path(configs[["results.dir"]], configs[["meta.dir"]]), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(configs[["results.dir"]], "train"), showWarnings = FALSE)
    dir.create(file.path(configs[["results.dir"]], "val"), showWarnings = FALSE)

    saveRDS(phe_train, file = file.path(configs[["results.dir"]], configs[["meta.dir"]], "phe_train.rds"))
    if (validation) saveRDS(phe_val, file = file.path(configs[["results.dir"]], configs[["meta.dir"]], "phe_val.rds"))
  }

  if (is.null(weight)) {
    weight <- rep(1, length(phenotype_names))
    names(weight) <- phenotype_names
  } else {
    if (is.null(names(weight))) {
      names(weight) <- phenotype_names
    }
    weight <- weight / sum(weight) * length(phenotype_names)
  }
  configs[["weight"]] <- weight

  if (standardize_response) {
    std_obj <- y_standardization(phe_train, phenotype_names, weight)
    phe_train <- std_obj$response
  }

  ## find covariates and response ##
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
    response_val <- as.matrix(phe_val[, phenotype_names, with = F])  ## validation responses are not standardized
    missing_response_val <- is.na(response_val)
    var_val <- apply(response_val, 2, function(x) mean((x - mean(x, na.rm = T))^2, na.rm = T))
  }

  ### --- Step 0: initial Y imputation based on the covariates --- ###
  fit_init <- initial_Y_imputation(response_train, covariates_train, missing_response_train)
  response_train <- fit_init$response
  object0 <- mean((fit_init$residual)^2)/2*ncol(response_train)

  rownames(fit_init$residual) <- ids_valid
  colnames(fit_init$residual) <- colnames(response_train)
  prod_full <- snpnet:::computeProduct(fit_init$residual, genotype_file, vars, stats, configs, iter=0) / length(rowIdx_subset_gen)
  # }
  score <- row_norm2(prod_full)
  score <- score / p.factor[names(score)]
  nobs <- n_subset_train
  nvars <- length(vars)-length(stats[["excludeSNP"]])
  configs[["lambda.min.ratio"]] <- lambda.min.ratio

  ## compute lambda sequence ##
  full_lams <- snpnet:::computeLambdas(score, nlambda, lambda.min.ratio)
  configs[["lambda"]] <- full_lams

  ### --- Start fitting --- ###
  fit_list <- vector("list", nlambda)

  active <- c()
  features_train <- NULL
  B_init <- NULL
  W_init <- NULL
  A_init <- NULL

  Z1 <- matrix(1, nrow(response_train), 1, dimnames = list(NULL, c("intercept")))
  if (!is.null(covariates_train)) Z1 <- cbind(Z1, as.matrix(covariates_train))
  PZ <- solve(crossprod(Z1), t(Z1))
  metric_train <- matrix(NA, nlambda, q_train)
  metric_val <- matrix(NA, nlambda, q_train)
  colnames(metric_train) <- colnames(metric_val) <- phenotype_names
  AUC_train <- matrix(NA, nlambda, length(binary_phenotypes))
  AUC_val <- matrix(NA, nlambda, length(binary_phenotypes))
  colnames(AUC_train) <- colnames(AUC_val) <- binary_phenotypes
  nactive <- rep(NA, 100)

  if (rank == ncol(response_train) && !is.null(covariates_train)) {
    features_train <- covariates_train
    if (validation) features_val <- covariates_val
  }

  if (prev_iter < 0) {
    prev_iter <- 0
    for (idx in 1:nlambda) {
      fname <- file.path(configs[["results.dir"]], paste0("output_lambda_", idx, ".RData"))
      if (file.exists(fname)) {
        prev_iter <- idx
      }
    }
  }

  if (prev_iter != 0) {
    cat("Recover iteration ", prev_iter, ". Now time: ", as.character(Sys.time()), "\n", sep = "")
    load_start <- Sys.time()
    new_configs <- configs
    new_pfactor <- p.factor
    load(file.path(configs[["results.dir"]], paste0("output_lambda_", prev_iter, ".RData")))
    check_configs_diff(configs, new_configs)
    if (!identical(new_pfactor, p.factor)) {
      warning("New p.factor and the saved p.factor are not the same. The new p.factor will be used.\n")  # to allow for running from old results; backward compatibility
    }
    p.factor <- new_pfactor
    configs <- new_configs
    start_lambda <- ilam + 1
    if (rank == ncol(response_train) && !is.null(covariates_train)) {
      features_train[, (feature_names) := snpnet:::prepareFeatures(chr_train, vars, feature_names, stats)]
      if (validation) features_val[, (feature_names) := snpnet:::prepareFeatures(chr_val, vars, feature_names, stats)]
    } else {
      features_train <- snpnet:::prepareFeatures(chr_train, vars, feature_names, stats)
      if (validation) {
        features_val <- snpnet:::prepareFeatures(chr_val, vars, feature_names, stats)
      }
    }

    if (rank == ncol(response_train)) {
      pred_train_0 <- sweep(safe_product(as.matrix(features_train), fit$CC), 2, fit$a0, FUN = "+")
    } else {
      if (!is.null(covariates_train)) {
        pred_train_0 <- sweep(safe_product(as.matrix(covariates_train), fit$W) + safe_product(as.matrix(features_train), fit$C), 2, fit$a0, FUN = "+")
      } else {
        pred_train_0 <- sweep(safe_product(as.matrix(features_train), fit$C), 2, fit$a0, FUN = "+")
      }
    }
    response_train[missing_response_train] <- pred_train_0[missing_response_train]

    prev_lambda <- nrow(metric_train)
    if (prev_lambda < nlambda) {
      metric_train <- rbind(metric_train, matrix(NA, nlambda-prev_lambda, q_train))
      metric_val <- rbind(metric_val, matrix(NA, nlambda-prev_lambda, q_train))
      AUC_train <- rbind(AUC_train, matrix(NA, nlambda-prev_lambda, length(binary_phenotypes)))
      AUC_val <- rbind(AUC_val, matrix(NA, nlambda-prev_lambda, length(binary_phenotypes)))
      nactive <- c(nactive, rep(NA, nlambda-prev_lambda))
    }
    load_end <- Sys.time()
    cat("Time elapsed on loading back features:", time_diff(load_start, load_end), "\n")
  } else {
    start_lambda <- 1
  }

  for (ilam in start_lambda:nlambda) {  # consider batch-type algorithm later
    cat("Current lambda:", ilam, "\n")
    lam <- full_lams[ilam]
    discard <- setdiff(colnames(features_train), c(active, covariate_names))
    if (ilam == 1) {
      norm_prod <- score
    }
    if (!is.null(features_train) && length(discard) > 0) {
      features_train[, (discard) := NULL]
      if (validation) features_val[, (discard) := NULL]
    }
    var_to_ignore <- active
    norm_prod[var_to_ignore] <- NA
    var_violate <- c()
    check <- FALSE
    while (!check) {
      if (length(var_violate) == 0) {
        order_norm <- order(norm_prod, decreasing = TRUE, na.last = NA)
        var_strong <- head(names(norm_prod[order_norm]), batch_size)
        feature_names_add <- var_strong
      } else {
        feature_names_add <- var_violate
      }
      if (length(feature_names_add) == 0) stop("Empty list of features to be added.")
      if (!is.null(features_train) && ncol(features_train) != 0) {
        features_train[, (feature_names_add) := snpnet:::prepareFeatures(chr_train, vars, feature_names_add, stats)]
        if (validation) features_val[, (feature_names_add) := snpnet:::prepareFeatures(chr_val, vars, feature_names_add, stats)]
      } else {
        features_train <- snpnet:::prepareFeatures(chr_train, vars, feature_names_add, stats)
        if (validation) features_val <- snpnet:::prepareFeatures(chr_val, vars, feature_names_add, stats)
      }
      if (rank == ncol(response_train) || ilam < 3) {
        if (rank == ncol(response_train)) {
          # penalty_factor <- rep(1, ncol(features_train))
          penalty_factor <- p.factor[colnames(features_train)]
        } else {
          features_train_combined <- cbind(covariates_train, features_train)
          # penalty_factor <- rep(1, ncol(features_train_combined))
          penalty_factor <- p.factor[colnames(features_train_combined)]
        }
        penalty_factor[covariate_names] <- 0
        lam_adjusted <- full_lams[ilam] * sum(penalty_factor) / length(penalty_factor)  # adjustment to counteract automatic normalization in glmnet
        if (rank == ncol(response_train)) {
          fit <- alternate_Y_glmnet(features_train, response_train, missing_response_train,
                                    lam_adjusted, penalty_factor, configs,
                                    num_covariates = length(covariate_names), r = rank, thresh = configs[["thresh"]], object0 = object0,
                                    W_init = W_init, B_init = B_init, A_init = A_init, glmnet_thresh = configs[["glmnet.thresh"]],
                                    max.iter = max.iter)
        } else {
          fit <- alternate_Y_glmnet(features_train_combined, response_train, missing_response_train,
                                    lam_adjusted, penalty_factor, configs,
                                    num_covariates = length(covariate_names), r = rank, thresh = configs[["thresh"]], object0 = object0,
                                    W_init = W_init, B_init = B_init, A_init = A_init, glmnet_thresh = configs[["glmnet.thresh"]],
                                    max.iter = max.iter)
        }
        W_init <- fit$W
        A_init <- fit$A
      } else {
        fit <- SRRR_iterative_missing_covariates(as.matrix(features_train), response_train,
                                                 missing_response_train, Z1, PZ, lam,
                                                 rank, max.iter, B_init, configs[["thresh"]], object0,
                                                 configs[["is.warm.start"]], configs[["is.A.converge"]],
                                                 glmnet_thresh = configs[["glmnet.thresh"]])
      }
      response_train <- fit$response
      residuals <- as.matrix(fit$residuals)
      rownames(residuals) <- ids_valid
      colnames(residuals) <- colnames(response_train)
      start_KKT <- Sys.time()
      cat("Start checking KKT condition ...\n")
      prod_resid <- snpnet:::computeProduct(residuals, genotype_file, vars, stats, configs, iter=0) / length(rowIdx_subset_gen)
      norm_prod <- row_norm2(prod_resid)
      norm_prod <- norm_prod / p.factor[names(norm_prod)]
      norm_prod[configs[["excludeSNP"]]] <- NA
      norm_prod_inner <- norm_prod
      norm_prod_inner[setdiff(colnames(features_train), covariate_names)] <- NA
      current_active <- which_row_active(fit$B)
      lam_active <- ifelse(length(current_active) > 0, max(norm_prod_inner[current_active]), lam)
      check <- all(norm_prod_inner <= lam_active, na.rm = T)
      var_violate <- names(norm_prod_inner)[which(norm_prod_inner > lam_active)]
      var_danger <- names(norm_prod_inner)[which(norm_prod_inner > lam & norm_prod_inner <= lam_active)]
      B_init <- fit$B
      end_KKT <- Sys.time()
      cat("Finish checking KKT condition: ", ifelse(check, "TRUE", "FALSE"), ". Number of violation variables: ",
          length(var_violate), ". Number of dangered variables: ", length(var_danger), ". Time elapsed: ",
          time_diff(start_KKT, end_KKT), "\n", sep = "")
    }

    pred_train <- response_train - residuals
    if (standardize_response) pred_train <- y_de_standardization(pred_train, std_obj$means, std_obj$sds, weight)
    MSE_train <- apply(pred_train - response_train_0, 2, function(x) mean((x-mean(x, na.rm=T))^2, na.rm = T))
    R2_train <- 1 - MSE_train / var_train
    cat("R2_train:\n")
    print(R2_train)
    metric_train[ilam, ] <- R2_train
    saveRDS(pred_train, file = file.path(configs[["results.dir"]], "train", paste0("pred_score_", ilam, ".rds")))

    if (validation) {
      if (rank == ncol(response_train)) {
        pred_val <- sweep(safe_product(as.matrix(features_val), fit$CC), 2, fit$a0, FUN = "+")
      } else {
        if (!is.null(covariates_val)) {
          pred_val <- sweep(safe_product(as.matrix(covariates_val), fit$W) + safe_product(as.matrix(features_val), fit$C), 2, fit$a0, FUN = "+")
        } else {
          pred_val <- sweep(safe_product(as.matrix(features_val), fit$C), 2, fit$a0, FUN = "+")
        }
      }
      if (standardize_response) pred_val <- y_de_standardization(pred_val, std_obj$means, std_obj$sds, weight)
      MSE_val <- apply(pred_val - response_val, 2, function(x) mean((x-mean(x, na.rm=T))^2, na.rm=T))
      R2_val <- 1 - MSE_val / var_val
      cat("R2_val:\n")
      print(R2_val)
      metric_val[ilam, ] <- R2_val
      saveRDS(pred_val, file = file.path(configs[["results.dir"]], "val", paste0("pred_score_", ilam, ".rds")))
    }

    if (length(binary_phenotypes) > 0) {
      for (bphe in binary_phenotypes) {
        data_logistic_train <- data.frame(response = response_train_0[, bphe], covariates_train, score = pred_train[, bphe])
        logitfit <- glm(response ~ ., data = data_logistic_train, family = binomial())
        pred_prob_train <- predict(logitfit, newdata = data_logistic_train, type = "response")
        not_missing_train <- !missing_response_train[, bphe]
        pred_obj <- ROCR::prediction(pred_prob_train[not_missing_train], response_train_0[not_missing_train, bphe])
        auc_obj <- ROCR::performance(pred_obj, measure = 'auc')
        AUC_train[ilam, bphe] <- auc_obj@y.values[[1]]
        if (validation) {
          data_logistic_val <- data.frame(covariates_val, score = pred_val[, bphe])
          pred_prob_val <- predict(logitfit, newdata = data_logistic_val, type = "response")
          not_missing_val <- !missing_response_val[, bphe]
          pred_obj <- ROCR::prediction(pred_prob_val[not_missing_val], response_val[not_missing_val, bphe])
          auc_obj <- ROCR::performance(pred_obj, measure = 'auc')
          AUC_val[ilam, bphe] <- auc_obj@y.values[[1]]
        }
      }
      cat("AUC_train:\n")
      print(AUC_train[ilam, ])
      if (validation) {
        cat("AUC_val:\n")
        print(AUC_val[ilam, ])
      }
    }

    fit[["stats"]] <- stats
    fit[["std_obj"]] <- std_obj
    fit[["weight"]] <- weight
    fit_list[[ilam]] <- fit
    current_active <- which_row_active(fit$B)
    nactive[ilam] <- length(current_active)
    cat("Number of active variables:", length(current_active), "\n")
    active <- unique(c(active, current_active))  ## ever-active set
    cat("Number of ever-active variables:", length(active), "\n")
    cat("Time since start of the procedure: ", time_diff(start_all, Sys.time()), "\n\n")
    B_init <- fit$B[active, , drop = F]
    norm_prod[active] <- NA

    if (save) {
      feature_names <- setdiff(colnames(features_train), covariate_names)
      save(fit, ilam, current_active, active, feature_names, norm_prod, B_init, W_init, A_init,
           metric_train, metric_val, AUC_train, AUC_val, nactive, weight, configs, p.factor,
           file = file.path(configs[["results.dir"]], paste0("output_lambda_", ilam, ".RData")))
      saveRDS(fit_list, file = file.path(configs[["results.dir"]], "fit_list.rds"))
    }

    if (early_stopping && ilam > 2 && all(metric_val[ilam, ] < metric_val[ilam-1, ]) &&
        all(metric_val[ilam-1, ] < metric_val[ilam-2, ])) {
      cat("None of the phenotype metrics is improving anymore. DONE. \n")
      break
    }

  }
  class(fit_list) <- "multisnpnet"
  fit_list
}
