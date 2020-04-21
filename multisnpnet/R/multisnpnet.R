multisnpnet <- function(genotype_file, phenotype_file, phenotype_names, covariate_names, results_dir,
                      r, nlambda = 100, batch_size = 100, lambda.min.ratio = 0.01,
                      max.iter = 10, is.warm.start = TRUE, is.A.converge = TRUE, thresh = 1e-7, glmnet_thresh = 1e-7,
                      standardize_response = FALSE,
                      configs, save = TRUE, validation = FALSE, genotype_file_val = NULL, early_stopping = FALSE,
                      prev_iter = 0, weight = NULL, binary_phenotypes = NULL,
                      use_plink2 = FALSE, genotype_p2file, genotype_p2file_val = NULL) {

  configs <- setup_configs_directories(configs, covariate_names, save, results_dir, validation, nlambda, lambda.min.ratio, standardize_response)

  start_all <- Sys.time()

  if (r > length(phenotype_names)) stop("The specified rank (", r, ") should not be greater than the number of responses (", length(phenotype_names), ").")

  cat("Start Sparse Reduced Rank Regression for ", paste(phenotype_names, collapse = ", "), ".\n", sep = "")

  ### ---- pre-process phenotypes and genotype data ---- ###
  phe_master <- fread(phenotype_file, colClasses = c("FID" = "character", "IID" = "character"), select = c("FID", "IID", covariate_names, phenotype_names))
  fill_missing(phe_master, phenotype_names, -9, NA) # replace -9 with NA
  # phe_master[, binary_phenotypes] <- phe_master[, binary_phenotypes]  ## should -1 !!!
  chr_train <- BEDMatrixPlus(genotype_file)
  n_chr_train <- nrow(chr_train)
  if (validation) {
    chr_val <- BEDMatrixPlus(genotype_file_val)
    n_chr_val <- nrow(chr_val)
  }
  q_train <- length(phenotype_names)
  # exclude samples that have missingness in all phenotypes, could change to some lower bound of phenotype missingness
  # ids_valid_phe <- phe_master$ID[apply(as.matrix(phe_master[, phenotype_names, with = F]), 1, function(x) any(!is.na(x)))]
  cat_ids <- paste(phe_master$FID, phe_master$IID, sep = "_")
  ids_valid_phe <- cat_ids[apply(as.matrix(phe_master[, phenotype_names, with = F]), 1, function(x) any(!is.na(x)))]
  # ids_valid_gen <- sapply(strsplit(rownames(chr_train), "_"), function(x) x[1])
  ids_valid_gen <- rownames(chr_train)
  ids_valid <- intersect(ids_valid_phe, ids_valid_gen)
  rowIdx_subset_gen <- match(ids_valid, ids_valid_gen)
  if (validation) {
    # ids_valid_gen_val <- sapply(strsplit(rownames(chr_val), "_"), function(x) x[1])
    ids_valid_gen_val <- rownames(chr_val)
    ids_valid_val <- intersect(ids_valid_phe, ids_valid_gen_val)
    rowIdx_subset_gen_val <- match(ids_valid_val, ids_valid_gen_val)
  }

  ## summary statistics: missing rate, mean, standard deviation (if needed) ##
  phe_train <- phe_master[match(ids_valid, cat_ids), ]
  if (validation) phe_val <- phe_master[match(ids_valid_val, cat_ids), ]

  if (!use_plink2) {
    stats <- snpnet:::computeStats(chr_train, rowIdx_subset_gen, c("pnas", "means", "sds"),
                                   path = file.path(results_dir, "meta"), save = save, configs = configs, verbose = TRUE, buffer.verbose = TRUE)
  } else {
    stats <- computeStats_P2(genotype_p2file, ids_valid, configs = configs)
    vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0('zstdcat ', paste0(genotype_p2file, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  }

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
    saveRDS(phe_train, file = file.path(results_dir, configs[["meta.dir"]], "phe_train.rds"))
    if (validation) saveRDS(phe_val, file = file.path(results_dir, configs[["meta.dir"]], "phe_val.rds"))
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

  # response_train_0[, binary_phenotypes] <- response_train_0[, binary_phenotypes] - 1

  if (standardize_response) {
    std_obj <- y_standardization(phe_train, phenotype_names, weight)
    phe_train <- std_obj$response
  }

  # browser()

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
  }
  # var_response <- var(as.numeric(response_train[!missing_response_train]))

  var_val <- apply(response_val, 2, function(x) mean((x - mean(x, na.rm = T))^2, na.rm = T))

  ### --- Step 0: initial Y imputation based on the covariates --- ###
  fit_init <- initial_Y_imputation(response_train, covariates_train, missing_response_train)
  response_train <- fit_init$response
  object0 <- mean((fit_init$residual)^2)/2*ncol(response_train)

  if (!use_plink2) {
    prod_full <- snpnet:::computeProduct(fit_init$residual, chr_train, rowIdx_subset_gen, stats,
                                         configs, verbose = TRUE, path = genotype_file)
  } else {
    rownames(fit_init$residual) <- ids_valid
    colnames(fit_init$residual) <- colnames(response_train)
    prod_full <- computeProduct_P2(fit_init$residual, genotype_p2file, vars, stats, configs, iter=0) / length(rowIdx_subset_gen)
  }
  score <- row_norm2(prod_full)
  if (is.null(lambda.min.ratio)) {
    lambda.min.ratio <- ifelse(n_subset_train < ncol(chr_train)-length(stats[["excludeSNP"]]), 0.01, 0.0001)
  }
  ## compute lambda sequence ##
  configs[["lambda.min.ratio"]] <- lambda.min.ratio
  configs[["nlambda"]] <- nlambda
  full_lams <- snpnet:::computeLambdas(score, nlambda, lambda.min.ratio)

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

  if (r == ncol(response_train) && !is.null(covariates_train)) {
    features_train <- covariates_train
    if (validation) features_val <- covariates_val
  }

  if (prev_iter < 0) {
    prev_iter <- 0
    for (idx in 1:nlambda) {
      fname <- file.path(results_dir, configs[["results.dir"]], paste0("output_lambda_", idx, ".RData"))
      if (file.exists(fname)) {
        prev_iter <- idx
      } else {
        break
      }
    }
  }

  if (prev_iter != 0) {
    cat("Recover iteration ", prev_iter, ". Now time: ", as.character(Sys.time()), "\n", sep = "")
    load_start <- Sys.time()
    load(file.path(results_dir, configs[["results.dir"]], paste0("output_lambda_", prev_iter, ".RData")))
    response_train <- fit$response
    start_lambda <- ilam + 1
    if (r == ncol(response_train) && !is.null(covariates_train)) {
      features_train[, (feature_names) := snpnet:::prepareFeatures(chr_train, feature_names, stats, rowIdx_subset_gen)]
      if (validation) features_val[, (feature_names) := snpnet:::prepareFeatures(chr_val, feature_names, stats, rowIdx_subset_gen_val)]
    } else {
      features_train <- snpnet:::prepareFeatures(chr_train, feature_names, stats, rowIdx_subset_gen)
      if (validation) {
        features_val <- snpnet:::prepareFeatures(chr_val, feature_names, stats, rowIdx_subset_gen_val)
      }
    }
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
  # browser()


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
      # B_init <- B_init[active, , drop = FALSE]
    }
    # var_strong <- c()
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
        features_train[, (feature_names_add) := snpnet:::prepareFeatures(chr_train, feature_names_add, stats, rowIdx_subset_gen)]
        if (validation) features_val[, (feature_names_add) := snpnet:::prepareFeatures(chr_val, feature_names_add, stats, rowIdx_subset_gen_val)]
      } else {
        features_train <- snpnet:::prepareFeatures(chr_train, feature_names_add, stats, rowIdx_subset_gen)
        if (validation) features_val <- snpnet:::prepareFeatures(chr_val, feature_names_add, stats, rowIdx_subset_gen_val)
      }
      if (r == ncol(response_train) || ilam < 3) {
        if (r == ncol(response_train)) {
          penalty_factor <- rep(1, ncol(features_train))
        } else {
          features_train_combined <- cbind(covariates_train, features_train)
          penalty_factor <- rep(1, ncol(features_train_combined))
        }
        penalty_factor[seq_len(length(covariate_names))] <- 0
        lam_adjusted <- full_lams[ilam] * sum(penalty_factor) / length(penalty_factor)  # adjustment to counteract automatic normalization in glmnet
        if (r == ncol(response_train)) {
          fit <- alternate_Y_glmnet(features_train, response_train, missing_response_train,
                                    lam_adjusted, penalty_factor, configs,
                                    num_covariates = length(covariate_names), r = r, thresh = thresh, object0 = object0,
                                    W_init = W_init, B_init = B_init, A_init = A_init, glmnet_thresh = glmnet_thresh)
        } else {
          fit <- alternate_Y_glmnet(features_train_combined, response_train, missing_response_train,
                                    lam_adjusted, penalty_factor, configs,
                                    num_covariates = length(covariate_names), r = r, thresh = thresh, object0 = object0,
                                    W_init = W_init, B_init = B_init, A_init = A_init, glmnet_thresh = glmnet_thresh)
        }
        W_init <- fit$W
        A_init <- fit$A

        # fit$C_full <- as.matrix(do.call(cbind, fit$glmnet.fit$beta))
        # idx_geno_C <- which(!(rownames(fit$C) %in% covariate_names))
        # fit$C <- fit$C_full[idx_geno_C, , drop = F]
        # svd_obj <- svd(fit$C)
        # fit$B <- svd_obj$u[, 1:r] %*% diag(svd_obj$d[1:r], r)
        # fit$A <- svd_obj$v[, 1:r]
        # fit$W <- fit$C_full[which(rownames(fit$C) %in% covariate_names), , drop = F]
        # fit$a0 <- fit$glmnet.fit$a0
      } else {
        # browser()
        fit <- SRRR_iterative_missing_covariates(as.matrix(features_train), response_train,
                                                 missing_response_train, Z1, PZ, lam,
                                                 r, max.iter, B_init,
                                                 thresh, object0, glmnet_thresh = glmnet_thresh)
      }
      response_train <- fit$response
      residuals <- as.matrix(fit$residuals)
      rownames(residuals) <- ids_valid
      colnames(residuals) <- colnames(response_train)
      start_KKT <- Sys.time()
      cat("Start checking KKT condition ...\n")
      if (!use_plink2) {
        prod_resid <- snpnet:::computeProduct(residuals, chr_train, rowIdx_subset_gen, stats, configs, path = "", verbose = TRUE)
      } else {
        prod_resid <- computeProduct_P2(residuals, genotype_p2file, vars, stats, configs, iter=0) / length(rowIdx_subset_gen)
      }
      # check.obj <- KKT_mcheck(residual.full, chr.train, rowIdx.subset.train, current.lams[start.lams:num.lams], ifelse(family == "gaussian" && use.glmnetPlus, 1, lambda.idx),
      # stats, glmfit, configs, verbose, KKT.verbose, path = file.path(genotype.dir, "train.bed"))
      # prod_resid <- crossprod(X, Y %*% fit$A - features_train %*% fit$B)
      norm_prod <- row_norm2(prod_resid)
      ### It's fine to overwrite norm_prod because if KKT check fails,
      ### only the violated variables will be added, the previous norm_prod is not useful anymore
      norm_prod[configs[["excludeSNP"]]] <- NA
      norm_prod_inner <- norm_prod
      norm_prod_inner[setdiff(colnames(features_train), covariate_names)] <- NA
      # norm_prod_with_NA <- norm_prod
      # norm_prod_with_NA[c(var_to_ignore, idx_strong)] <- NA
      current_active <- which_row_active(fit$B)
      lam_active <- ifelse(length(current_active) > 0, max(norm_prod_inner[current_active]), lam)
      check <- all(norm_prod_inner <= lam_active, na.rm = T)
      var_violate <- names(norm_prod_inner)[which(norm_prod_inner > lam_active)]
      # var_violate <- which(norm_prod_with_NA > max(n*lam, max_norm_strong))
      # num_violations[ilam] <- num_violations[ilam] + 1
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
    saveRDS(pred_train, file = file.path(results_dir, configs[["results.dir"]], "train", paste0("pred_score_", ilam, ".rds")))

    if (validation) {
      if (r == ncol(response_train)) {
        pred_val <- sweep(as.matrix(features_val) %*% fit$CC, 2, fit$a0, FUN = "+")
      } else {
        if (!is.null(covariates_val)) {
          pred_val <- sweep(as.matrix(covariates_val) %*% fit$W + as.matrix(features_val) %*% fit$C, 2, fit$a0, FUN = "+")
        } else {
          pred_val <- sweep(as.matrix(features_val) %*% fit$C, 2, fit$a0, FUN = "+")
        }
      }
      # browser()
      if (standardize_response) pred_val <- y_de_standardization(pred_val, std_obj$means, std_obj$sds, weight)
      # metric_val[ilam] <- mean((pred_val[!missing_response_val] - response_val[!missing_response_val])^2)
      MSE_val <- apply(pred_val - response_val, 2, function(x) mean((x-mean(x, na.rm=T))^2, na.rm=T))
      R2_val <- 1 - MSE_val / var_val
      cat("R2_val:\n")
      print(R2_val)
      metric_val[ilam, ] <- R2_val
      saveRDS(pred_val, file = file.path(results_dir, configs[["results.dir"]], "val", paste0("pred_score_", ilam, ".rds")))
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
           metric_train, metric_val, AUC_train, AUC_val, nactive, weight, configs,
           file = file.path(results_dir, configs[["results.dir"]], paste0("output_lambda_", ilam, ".RData")))
    }

    if (early_stopping && ilam > 2 && all(metric_val[ilam, ] < metric_val[ilam-1, ]) &&
        all(metric_val[ilam-1, ] < metric_val[ilam-2, ])) {
      cat("None of the phenotype metrics is improving anymore. DONE. \n")
      break
    }

  }
  out <- list(fit = fit_list)
  out
}
