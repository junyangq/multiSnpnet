fill_missing <- function(data, colnames, key, values) {
  if (length(values) == 1) values <- rep(values, length(colnames))
  if (is.na(key)) {
    for (j in 1:length(colnames)) {
      set(data, i=which(is.na(data[[colnames[j]]])), j=colnames[j], value=values[j])
    }
  } else {
    for (j in 1:length(colnames)) {
      set(data, i=which(data[[colnames[j]]] == key), j=colnames[j], value=values[j])
    }
  }
}

initial_Y_imputation <- function(response, covariates, missing_response) {
  residual <- matrix(NA, nrow(response), ncol(response))
  colnames(residual) <- colnames(response)
  if (is.null(covariates)) covariates <- data.frame(intercept = rep(1, nrow(response)))
  for (k in 1:ncol(response)) {
    fit <- lm(response[, k] ~ ., data = covariates)
    pred <- predict(fit, newdata = covariates)
    response[missing_response[, k], k] <- pred[missing_response[, k]]
    residual[, k] <- response[, k] - pred
  }
  out <- list(response = response, residual = residual)
  out
}


time_diff <- function(start_time, end_time) {
  paste(round(end_time-start_time, 4), units(end_time-start_time))
}


alternate_Y_glmnet <- function(features, response, missing_response, lambda, penalty_factor, configs, 
                               num_covariates, r, thresh = 1E-7, object0, W_init, B_init, A_init, glmnet_thresh = 1e-7) {
  # browser()
  converge <- FALSE
  features_matrix <- as.matrix(features)
  niter <- 0
  obj_values <- c()
  cat("    Start Y-C (glmnet) iteration ...\n")
  message <- "Terminated"
  
  if (is.null(W_init) || is.null(B_init) || is.null(A_init)) {
    CC <- NULL
  } else {
    CC <- matrix(0, ncol(features_matrix), ncol(response))
    rownames(CC) <- colnames(features_matrix)
    # browser()
    CC[rownames(B_init), ] <- tcrossprod(as.matrix(B_init), A_init)
    CC[rownames(W_init), ] <- as.matrix(W_init)
  }
  
  while (!converge) {
    niter <- niter + 1
    fit <- glmnetPlus::glmnet(features_matrix, response, family = "mgaussian", lambda = lambda, penalty.factor = penalty_factor,
                              standardize = configs[["standardize.variant"]], standardize.response = FALSE, beta0 = CC, thresh = glmnet_thresh)
    CC <- do.call(cbind, fit$beta)
    pred <- predict(fit, newx = features_matrix, type = "response")[, , 1]
    response[missing_response] <- pred[missing_response]
    obj_values[niter] <- sum((response-pred)^2) / 2 / nrow(response) + lambda * sum(apply(CC, 1, function(x) sqrt(sum(x^2))))
    cat("         Objective (", niter, "): ", round(obj_values[niter], digits = 8), "\n", sep = "")
    if (niter > 1) {
      delta <- obj_values[niter] - obj_values[niter-1]
    }
    if (niter > 1 && delta < thresh*object0) {
      message <- "Converged"
      converge <- TRUE
    }
  }
  cat("    Finish Y-C (glmnet) iteration:", message, "in", niter, "iterations.\n")
  colnames(CC) <- colnames(response)
  W <- CC[1:num_covariates, , drop = F]
  C <- CC[(num_covariates+1):ncol(features), , drop = F]
  if (r < ncol(response)) {  ## but assumes that the number of selected variables is less than the target rank
    svd_obj <- svd(C)
    B <- svd_obj$u[, 1:r, drop = F] %*% diag(svd_obj$d[1:r], r)
    rownames(B) <- rownames(C)
    A <- svd_obj$v[, 1:r, drop = F]
    rownames(A) <- colnames(C)
  } else {
    B <- C
    A <- diag(1, nrow = r)
    rownames(A) <- colnames(C)
  }
  a0 <- fit$a0
  # residuals <- response - predict(fit, newx = features_matrix)[, , 1]
  residuals <- response - pred
  out <- list(response = response, a0 = a0, W = W, C = C, CC = CC, B = B, A = A, residuals = residuals, obj_values = obj_values)
  out
}


SRRR_iterative_missing_covariates <- function(X, Y, Y_missing, Z, PZ, lambda, r, niter, B0, thresh = 1e-7, object0, 
                                              is.warm.start = FALSE, is.A.converge = TRUE, glmnet_thresh = 1e-7) {
  # browser()
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  # browser()
  # B <- B0
  B <- matrix(0, ncol(X), r)
  rownames(B) <- colnames(X)
  # browser()
  B[rownames(B0), ] <- as.matrix(B0)
  
  obj_values <- rep(NA, niter)
  A_niter <- rep(NA, niter) 
  message <- "Terminated"
  
  count <- 0
  
  start_BAY <- Sys.time()
  
  for (k in 1:niter) {
    count <- count + 1
    A_niter[k] <- 0
    # fix B, solve A
    score <- as.matrix(X %*% B)
    # crossmat <- crossprod(crossprod(X, Y), B)
    # is_Y_converge <- FALSE
    impute_iter_count <- 0
    # fit_score <- 
    projected_score_Z <- Z %*% (PZ %*% score)
    RS <- score - projected_score_Z
    # RS <- residuals(fit_score)
    cat("    Start Y-A iteration ...\n")
    B_norm <- sum(row_norm2(B))
    start_Y_A <- Sys.time()
    while (TRUE) {
      A_niter[k] <- A_niter[k] + 1
      impute_iter_count <- impute_iter_count + 1
      # print(impute_iter_count)
      # fit_Y <- lm(Y ~ Z)
      # RY <- residuals(fit_Y)
      projected_Y_Z <- Z %*% (PZ %*% Y)
      RY <- Y - projected_Y_Z
      crossmat <- crossprod(RY, score)
      svd_cross <- svd(crossmat)
      # if (impute_iter_count > 1) A_old <- A
      if (impute_iter_count > 1) obj_old <- obj
      A <- tcrossprod(svd_cross$u, svd_cross$v)
      # Y_new <- tcrossprod(score, A)
      # Y_new <- Y - predict(fit_Y, )
      Y_new <- projected_Y_Z + tcrossprod(RS, A)  # implicit W
      # print(A[1, ])
      # if (impute_iter_count > 1) {
      #   delta <- mean((A - A_old)^2)  ## until A is stable
      # } else {
      #   delta <- 100
      # }
      # print(delta)
      # delta <- mean((Y[Y_missing] - Y_new[Y_missing])^2)
      # print(delta)
      # browser()
      Y[Y_missing] <- Y_new[Y_missing]
      
      obj <- mean((Y - Y_new)^2)/2*ncol(Y) + lambda * B_norm
      cat("         Objective (", impute_iter_count, "): ", round(obj, digits = 8), "\n", sep = "")
      if (impute_iter_count > 1) {
        delta <- obj_old - obj
      } else {
        delta <- 100
      }
      
      if (delta < thresh*object0 || (!is.A.converge && A_niter[k] > 0)) {
        break
      }
    }
    end_Y_A <- Sys.time()
    cat("    Finish Y-A iteration: ", impute_iter_count, ". Time elapsed: ", 
        time_diff(start_Y_A, end_Y_A), "\n", sep = "")
    
    cat("    Start solving for B ...\n")
    
    # browser()
    # fix A, solve B
    # fit_Y <- lm(Y ~ Z)
    # coef_Y_Z <- PZ %*% Y
    # projected_Y_Z <- Z %*% coef_Y_Z
    # RY <- residuals(fit_Y)
    # RY <- Y - projected_Y_Z  # more recent
    # YA <- RY %*% A + projected_score_Z
    # YA <- Y %*% A
    ZW <- Y_new - tcrossprod(score, A)
    # YA <- tcrossprod(score, A)
    # ZW <- projected_Y_Z - tcrossprod(projected_score_Z, A)
    YA <- (Y - ZW) %*% A
    
    # YA <- Y %*% A
    # keep a strong set, update only when KKT is violated ...
    # if (!is.null(lambda_seq)) {
    #   lambda_seq <- sort(unique(c(lambda_seq, lambda)), decreasing = T)
    #   lambda_seq <- lambda_seq[lambda_seq >= lambda]
    # }
    if (k == 1 && !is.warm.start) {
      mfit <- glmnetPlus::glmnet(x = X, y = YA, family = "mgaussian", standardize = F, intercept = F, lambda = lambda, thresh = glmnet_thresh)
      # print(mfit$npasses)
    } else {
      mfit <- glmnetPlus::glmnet(x = X, y = YA, family = "mgaussian", standardize = F, intercept = F, lambda = lambda, beta0 = B, thresh = glmnet_thresh)
      # print(mfit$npasses)
    }
    beta_single <- coef(mfit, s = lambda, x = X, y = YA)
    B <- do.call(cbind, beta_single)[-1, ]
    
    end_B <- Sys.time()
    cat("    Finish solving for B. ", "Time elapsed: ", 
        time_diff(end_Y_A, end_B), "\n", sep = "")
    
    
    if (k > 1) C_old <- C
    C <- tcrossprod(B, A)
    
    MAXLEN <- 2^31 - 1  # deal with long vector
    ncol.chunk <- floor(MAXLEN / as.double(nrow(X)) / 4)
    numChunks <- ceiling(ncol(X) / ncol.chunk)
    
    for (jc in 1:numChunks) {
      idx <- ((jc-1)*ncol.chunk+1):min(jc*ncol.chunk, ncol(X))
      if (jc == 1) {
        score <- as.matrix(X[, idx] %*% B[idx, ])
      } else {
        score <- score + as.matrix(X[, idx] %*% B[idx, ])
      }
    }
    
    # score <- X %*% B
    
    Y_new <- ZW + tcrossprod(score, A)
    Y[Y_missing] <- Y_new[Y_missing]
    residuals <- Y - Y_new
    # obj_values[k] <- 1/(2*n) * sum((residuals[!Y_missing])^2) +
    #   lambda * sum(apply(B, 1, function(x) sqrt(sum(x^2))))
    obj_values[k] <- 1/(2*n) * sum((residuals)^2) +
      lambda * sum(row_norm2(B))
    
    # if (k > 1) print(abs(obj_values[k] - obj_values[k-1])/object0)
    if (k > 1 && abs(obj_values[k] - obj_values[k-1]) < thresh*object0) {
    # if (k > 1 && mean((C - C_old)^2) < thresh) {  ## until C is stable
      message <- "Converged"
      obj_values <- obj_values[1:k]
      A_niter <- A_niter[1:k]
      break
    }
  }
  
  end_BAY <- Sys.time()
  cat("Finish B-A-Y iteration: ", message, " after ", k, " iterations. ", "Time elapsed: ", 
      time_diff(start_BAY, end_BAY), "\n", sep = "")
  
  coef_score_Z <- PZ %*% score
  coef_Y_Z <- PZ %*% Y
  W_full <- coef_Y_Z - tcrossprod(coef_score_Z, A) ## can be different from W on exit of inner loop
  a0 <- W_full[1, ]
  W <- W_full[-1, , drop = F]
  
  # browser()
  out <- list(B = B, A = A, C = C, a0 = a0, W = W, obj_values = obj_values, message = message, 
              niter = count, response = Y, A_niter = A_niter, residuals = residuals)
  out
}


row_norm2 <- function(X) {
  out <- apply(X, 1, function(x) sqrt(sum(x^2)))
  out
}

which_row_active <- function(X) {
  out <- names(which(apply(X, 1, function(x) any(x != 0))))
  out
}

MSE_col <- function(X) {
  MSE <- apply(X, 2, function(x) mean((x - mean(x, na.rm = T))^2, na.rm = T))
  MSE
}

# We assume that each column of Y has been standardized

MSE_coef <- function(C_hat, C) {
  C_full <- matrix(0, nrow(C), ncol(C))
  rownames(C_full) <- rownames(C)
  C_full[rownames(C_hat), ] <- C_hat
}

y_standardization <- function(response, phenotype_names, weight) {
  y_means <- c()
  y_sds <- c()
  for (name in phenotype_names) {
    y_means[name] <- mean(response[[name]], na.rm = T)
    y_sds[name] <- sd(response[[name]], na.rm = T)
    response[[name]] <- (response[[name]] - y_means[name]) / y_sds[name] * sqrt(weight[name])
  }
  out <- list(response = response, means = y_means, sds = y_sds)
  out
}

y_de_standardization <- function(response, means, sds, weight) {
  for (name in colnames(response)) {
    response[, name] <- sds[name] * response[, name] / sqrt(weight[name]) + means[name]
  }
  response
}

setup_configs_directories <- function(configs, covariates, save, results.dir, validation) {
  configs[["covariates"]] <- covariates
  default_settings <- c(missing.rate = 0.1, MAF.thresh = 0.001, nCores = 1,
                        nlams.init = 10, nlams.delta = 5)
  for (name in names(default_settings)) {
    if (!(name %in% names(configs))) configs[[name]] <- as.numeric(default_settings[name])
  }
  if (!("bufferSize" %in% names(configs)))
    stop("bufferSize should be provided to guide the memory capacity.")
  if (!("chunkSize" %in% names(configs)))
    configs[["chunkSize"]] <- configs[["bufferSize"]] / configs[["nCores"]]
  if (save) {
    if (is.null(configs[["meta.dir"]])) configs[["meta.dir"]] <- "meta/"
    if (is.null(configs[["results.dir"]])) configs[["results.dir"]] <- "results/"
    dir.create(file.path(results.dir, configs[["meta.dir"]]), showWarnings = FALSE, recursive = T)
    dir.create(file.path(results.dir, configs[["results.dir"]], "train"), showWarnings = FALSE, recursive = T)
    if (validation) dir.create(file.path(results.dir, configs[["results.dir"]], "val"), showWarnings = FALSE, recursive = T)
  }
  configs
}

SRRR_path <- function(genotype_file, phenotype_file, phenotype_names, covariate_names, results_dir, 
                      r, nlambda = 100, batch_size = 100, lambda.min.ratio = 0.01, 
                      max.iter = 10, is.warm.start = TRUE, is.A.converge = TRUE, thresh = 1e-7, glmnet_thresh = 1e-7,
                      standardize_response = FALSE, 
                      configs, save = TRUE, validation = FALSE, genotype_file_val = NULL, early_stopping = FALSE,
                      prev_iter = 0, weight = NULL, binary_phenotypes = NULL) {
  
  configs <- setup_configs_directories(configs, covariate_names, save, results_dir, validation)
  
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
  stats <- snpnet:::computeStats(chr_train, rowIdx_subset_gen, c("pnas", "means", "sds"),
                                 path = file.path(results_dir, "meta"), save = save, configs = configs, verbose = TRUE, buffer.verbose = TRUE)
  
  phe_train <- phe_master[match(ids_valid, cat_ids), ]
  if (validation) phe_val <- phe_master[match(ids_valid_val, cat_ids), ]
  
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
  prod_full <- snpnet:::computeProduct(fit_init$residual, chr_train, rowIdx_subset_gen, stats, 
                                       configs, verbose = TRUE, path = genotype_file)
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
  
  Z1 <- cbind(intercept = 1, as.matrix(covariates_train))
  PZ <- solve(crossprod(Z1), t(Z1))
  metric_train <- matrix(NA, nlambda, q_train)
  metric_val <- matrix(NA, nlambda, q_train)
  colnames(metric_train) <- colnames(metric_val) <- phenotype_names
  AUC_train <- matrix(NA, nlambda, length(binary_phenotypes))
  AUC_val <- matrix(NA, nlambda, length(binary_phenotypes))
  colnames(AUC_train) <- colnames(AUC_val) <- binary_phenotypes
  nactive <- rep(NA, 100)
  
  if (r == ncol(response_train)) {
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
    if (r == ncol(response_train)) {
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
                                                 missing_response_train, cbind(intercept = 1, as.matrix(covariates_train)), PZ, lam, 
                                                 r, max.iter, B_init, 
                                                 thresh, object0, glmnet_thresh = glmnet_thresh)
      }
      response_train <- fit$response
      residuals <- as.matrix(fit$residuals)
      start_KKT <- Sys.time()
      cat("Start checking KKT condition ...\n")
      prod_resid <- snpnet:::computeProduct(residuals, chr_train, rowIdx_subset_gen, stats, configs, path = "", verbose = TRUE)
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
        pred_val <- sweep(as.matrix(covariates_val) %*% fit$W + as.matrix(features_val) %*% fit$C, 2, fit$a0, FUN = "+")
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
           metric_train, metric_val, AUC_train, AUC_val, nactive, 
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

