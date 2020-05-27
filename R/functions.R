setupMultiConfigs <- function(configs, genotype_file, phenotype_file, phenotype_names, covariate_names,
                              nlambda, mem,
                              standardize_response, max.iter, rank, prev_iter, batch_size) {
  out.args <- as.list(environment())
  defaults_multi <- list(
    missing.rate = 0.1,
    MAF.thresh = 0.001,
    nCores = 1,
    glmnet.thresh = 1e-07,
    nlams.init = 10,
    nlams.delta = 5,
    vzs=TRUE, # geno.pfile vzs
    increase.size = NULL,
    standardize.variant = FALSE,
    early.stopping = TRUE,
    stopping.lag = 2,
    niter = 10,
    lambda.min.ratio = NULL,
    KKT.verbose = FALSE,
    use.glmnetPlus = NULL,
    save = FALSE,
    save.computeProduct = FALSE,
    prevIter = 0,
    results.dir = NULL,
    meta.dir = 'meta',
    save.dir = 'results',
    verbose = FALSE,
    KKT.check.aggressive.experimental = FALSE,
    gcount.basename.prefix = 'snpnet.train',
    gcount.full.prefix=NULL,
    endian="little",
    metric=NULL,
    plink2.path='plink2',
    zstdcat.path='zstdcat',
    rank = TRUE,
    is.warm.start = TRUE,
    is.A.converge = TRUE,
    thresh = 1e-7
  )
  for (name in setdiff(names(out.args), "configs")) {
    configs[[name]] <- out.args[[name]]
  }
  for (name in names(defaults_multi)) {
    if (!(name %in% names(configs))) {
      configs[[name]] <- defaults_multi[[name]]
    }
  }

  # update settings
  if(is.null(configs[['increase.size']]))  configs[['increase.size']] <- configs[['batch_size']]/2

  # We will write some intermediate files to meta.dir and save.dir.
  # those files will be deleted with snpnet::cleanUpIntermediateFiles() function.
  if (is.null(configs[['results.dir']])) configs[['results.dir']] <- tempdir(check = TRUE)
  dir.create(file.path(configs[['results.dir']], configs[["meta.dir"]]), showWarnings = FALSE, recursive = T)
  dir.create(file.path(configs[['results.dir']], configs[["save.dir"]]), showWarnings = FALSE, recursive = T)
  if(is.null(configs[['gcount.full.prefix']])) configs[['gcount.full.prefix']] <- file.path(
    configs[['results.dir']], configs[["meta.dir"]], configs['gcount.basename.prefix']
  )

  configs
}

#' @importFrom data.table set
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
                               num_covariates, r, thresh = 1E-7, object0, W_init, B_init, A_init, glmnet_thresh = 1e-7,
                               max.iter) {
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
    CC[rownames(B_init), ] <- tcrossprod(as.matrix(B_init), A_init)
    CC[rownames(W_init), ] <- as.matrix(W_init)
  }

  while (!converge || niter > max.iter) {
    niter <- niter + 1
    fit <- glmnetPlus::glmnet(features_matrix, response, family = "mgaussian", lambda = lambda, penalty.factor = penalty_factor,
                              standardize = configs[["standardize.variant"]], standardize.response = FALSE, beta0 = CC, thresh = glmnet_thresh)
    CC <- do.call(cbind, fit$beta)
    pred <- glmnetPlus::predict.mrelnet(fit, newx = features_matrix, type = "response")[, , 1]
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
  W <- CC[seq_len(num_covariates), , drop = F]
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
  residuals <- response - pred
  out <- list(response = response, a0 = a0, W = W, C = C, CC = CC, B = B, A = A, residuals = residuals, obj_values = obj_values)
  out
}


SRRR_iterative_missing_covariates <- function(X, Y, Y_missing, Z, PZ, lambda, r, niter, B0, thresh = 1e-7, object0,
                                              is.warm.start = FALSE, is.A.converge = TRUE, glmnet_thresh = 1e-7) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  B <- matrix(0, ncol(X), r)
  rownames(B) <- colnames(X)
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

    MAXLEN <- 2^31 - 1  # deal with long vector
    ncol.chunk <- floor(MAXLEN / as.double(nrow(X)) / 4)
    numChunks <- ceiling(ncol(X) / ncol.chunk)

    for (jc in 1:numChunks) {
      idx <- ((jc-1)*ncol.chunk+1):min(jc*ncol.chunk, ncol(X))
      if (jc == 1) {
        score <- as.matrix(X[, idx] %*% B[idx, , drop = F])
      } else {
        score <- score + as.matrix(X[, idx] %*% B[idx, , drop = F])
      }
    }

    impute_iter_count <- 0
    projected_score_Z <- Z %*% (PZ %*% score)
    RS <- score - projected_score_Z
    cat("    Start Y-A iteration ...\n")
    B_norm <- sum(row_norm2(B))
    start_Y_A <- Sys.time()
    while (TRUE) {
      A_niter[k] <- A_niter[k] + 1
      impute_iter_count <- impute_iter_count + 1
      projected_Y_Z <- Z %*% (PZ %*% Y)
      RY <- Y - projected_Y_Z
      crossmat <- crossprod(RY, score)
      svd_cross <- svd(crossmat)
      if (impute_iter_count > 1) obj_old <- obj
      A <- tcrossprod(svd_cross$u, svd_cross$v)
      Y_new <- projected_Y_Z + tcrossprod(RS, A)  # implicit W
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

    ZW <- Y_new - tcrossprod(score, A)
    YA <- (Y - ZW) %*% A
    if (k == 1 && !is.warm.start) {
      mfit <- glmnetPlus::glmnet(x = X, y = YA, family = "mgaussian", standardize = F, intercept = F, lambda = lambda, thresh = glmnet_thresh)
    } else {
      mfit <- glmnetPlus::glmnet(x = X, y = YA, family = "mgaussian", standardize = F, intercept = F, lambda = lambda, beta0 = B, thresh = glmnet_thresh)
    }
    if (is.null(dim(mfit$a0))) {
      mfit$a0 <- matrix(mfit$a0, nrow = 1)
      mfit$beta <- list(mfit$beta)
    }
    beta_single <- coef(mfit, s = lambda, x = X, y = YA)
    B <- do.call(cbind, beta_single)[-1, , drop = F]

    end_B <- Sys.time()
    cat("    Finish solving for B. ", "Time elapsed: ",
        time_diff(end_Y_A, end_B), "\n", sep = "")


    if (k > 1) C_old <- C
    C <- tcrossprod(as.matrix(B), A)

    MAXLEN <- 2^31 - 1  # deal with long vector
    ncol.chunk <- floor(MAXLEN / as.double(nrow(X)) / 4)
    numChunks <- ceiling(ncol(X) / ncol.chunk)

    for (jc in 1:numChunks) {
      idx <- ((jc-1)*ncol.chunk+1):min(jc*ncol.chunk, ncol(X))
      if (jc == 1) {
        score <- as.matrix(X[, idx, drop = F] %*% B[idx, , drop = F])
      } else {
        score <- score + as.matrix(X[, idx, drop = F] %*% B[idx, , drop = F])
      }
    }

    Y_new <- ZW + tcrossprod(score, A)
    Y[Y_missing] <- Y_new[Y_missing]
    residuals <- Y - Y_new
    obj_values[k] <- 1/(2*n) * sum((residuals)^2) +
      lambda * sum(row_norm2(B))

    if (k > 1 && abs(obj_values[k] - obj_values[k-1]) < thresh*object0) {
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

  colnames(C) <- colnames(Y)

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

timeDiff <- function(start.time, end.time = NULL) {
  if (is.null(end.time)) end.time <- Sys.time()
  paste(round(end.time-start.time, 4), units(end.time-start.time))
}

# used to compute the new lambda.min.ratio if we want to extend the original lambda sequence
compute_lambda_min_ratio <- function(nlambda.new, nlambda = 100, ratio = 0.01) {
  exp((nlambda.new-1)/(nlambda-1)*log(ratio))
}

check_configs_diff <- function(old_configs, new_configs) {
  msg <- ""
  for (name in intersect(names(old_configs), names(new_configs))) {
    if (!identical(old_configs[[name]], new_configs[[name]])) {
      msg <- paste0(
        msg,
        cat("Changed config for ", name, ": ", old_configs[[name]], " -> ", new_configs[[name]], "\n", sep = "")
        )
    }
  }
  for (name in setdiff(names(old_configs), names(new_configs))) {
    if (!identical(old_configs[[name]], new_configs[[name]])) {
      msg <- paste0(
        msg,
        cat("Deleted config for ", name, ": ", old_configs[[name]], "\n", sep = "")
      )
    }
  }
  for (name in setdiff(names(new_configs), names(old_configs))) {
    if (!identical(old_configs[[name]], new_configs[[name]])) {
      msg <- paste0(
        msg,
        cat("Added config for ", name, ": ", new_configs[[name]], "\n", sep = "")
      )
    }
  }
  if (msg != "") {
    warning(msg)
  }
}

#' Extract Coefficients from the Fitted Object or File
#'
#' @param fit Fit object returned from multisnpnet
#' @param fit_path Path to the file that saves the fit object
#' @param idx Lambda indices where the coefficients are requested
#' @param uv Boolean. Whether U, V are used to represent the decomposed matrices or B, A.
#'
#' @return List of coefficients where each element is a list of one type of coefficients over the
#'   provided lambda indices.
#'
#' @export
coef_multisnpnet <- function(fit = NULL, fit_path = NULL, idx = NULL, uv = TRUE) {
  if (is.null(fit) && is.null(fit_path)) {
    stop("Either fit object or file path to the saved object should be provided.\n")
  }
  if (is.null(fit)) fit <- readRDS(file = fit_path)
  if (is.null(idx)) idx <- seq_along(fit)
  fit <- fit[idx]
  a0 <- lapply(fit, function(x) x$a0)
  W <- lapply(fit, function(x) x$W)
  B <- lapply(fit, function(x) x$B)
  A <- lapply(fit, function(x) x$A)
  if (uv) {
    out <- list(a0 = a0, W = W, U = B, V = A, idx = idx)
  } else {
    out <- list(a0 = a0, W = W, B = B, A = A, idx = idx)
  }
  out
}

#' Predict from the Fitted Object or File
#'
#' @param fit List of fit object returned from multisnpnet.
#' @param saved_path Path to the file that saves the fit object. The full path is constructed as ${saved_path}${idx}.RData.
#' @param new_genotype_file Path to the new suite of genotype files. genotype_file.{pgen, psam, pvar.zst}.
#'   must exist.
#' @param new_phenotype_file Path to the phenotype. The header must include FID, IID.
#' @param idx Lambda indices where the coefficients are requested.
#' @param covariate_names Character vector of the names of the adjustment covariates.
#' @param split_col Name of the split column. If NULL, all samples will be used.
#' @param split_name Vector of split labels where prediction is to be made.
#' @param binary_phenotypes Vector of names of the binary phenotypes. If training split is provided,
#'   logistic regression will be refitted on the covariates and linear prediction score (from
#'   multivariate fit) and the final prediction updated. In addition, AUC will be computed for
#'   binary phenotypes.
#' @param zstdcat_path Path to zstdcat program, needed when loading variants
#'
#' @return A list containing the prediction and the resopnse for which the prediction is made.
#'
#' @export
predict_multisnpnet <- function(fit = NULL, saved_path = NULL, new_genotype_file, new_phenotype_file,
                                idx = NULL, covariate_names = NULL, split_col = NULL, split_name = NULL,
                                binary_phenotypes = NULL, zstdcat_path = "zstdcat") {
  if (is.null(fit) && is.null(saved_path)) {
    stop("Either fit object or file path to the saved object should be provided.\n")
  }
  if (is.null(fit) && is.null(idx)) {
    stop("Lambda indices on which prediction is made must be provided.\n")
  }
  if (is.null(fit)) {
    # last <- max(idx)
    # latest_result <- paste0(saved_path, last, ".RData")
    # e <- new.env()
    # load(latest_result, envir = e)
    # feature_names <- e$active
    #
    fit <- vector("list", length(idx))
    for (i in seq_along(idx)) {
      e <- new.env()
      load(paste0(saved_path, idx[i], ".RData"), envir = e)
      fit[[i]] <- e$fit
    }
  }
  stats <- fit[[length(fit)]][["stats"]]
  std_obj <- fit[[length(fit)]][["std_obj"]]
  weight <- fit[[length(fit)]][["weight"]]
  phenotype_names <- colnames(as.matrix(fit[[length(fit)]][["C"]]))
  is_full_rank <- (ncol(as.matrix(fit[[length(fit)]][["B"]])) == ncol(as.matrix(fit[[length(fit)]][["C"]])))

  covariate_names_fit <- rownames(fit[[length(fit)]][["W"]])
  if (!setequal(covariate_names_fit, covariate_names)) {
    stop("Unequal covariate sets in the fit and the argument.\n",
         "Fit: ", covariate_names_fit, "\n",
         "Argument: ", covariate_names, "\n")
  }

  ids <- list()
  ids[["psam"]] <- snpnet:::readIDsFromPsam(paste0(new_genotype_file, '.psam'))

  ctype <- c("FID" = "character", "IID" = "character")
  if (!is.null(split_col)) ctype[split_col] <- "character"
  phe_master <- data.table::fread(new_phenotype_file, colClasses = ctype, select = c("FID", "IID", split_col, covariate_names, phenotype_names))
  if (length(covariate_names) > 0) {
    cov_master <- as.matrix(phe_master[, covariate_names, with = F])
    cov_no_missing <- apply(cov_master, 1, function(x) all(!is.na(x)))
    phe_master <- phe_master[cov_no_missing, ]
  }
  phe_master[["ID"]] <- paste(phe_master[["FID"]], phe_master[["IID"]], sep = "_")
  phe_master <- phe_master %>%
    dplyr::inner_join(data.frame(ID = ids[["psam"]], sort_order = seq_along(ids[["psam"]]), stringsAsFactors = FALSE), by = "ID") %>%
    dplyr::arrange(sort_order) %>% dplyr::select(-sort_order) %>%
    data.table::as.data.table()

  fill_missing(phe_master, phenotype_names, -9, NA)

  if (is.null(split_col)) {
    split_name <- "train"
    ids[["train"]] <- phe_master$ID
  } else {
    for (split in split_name) {
      ids[[split]] <- phe_master$ID[phe_master[[split_col]] == split]
      if (length(ids[[split]]) == 0) {
        warning(paste("Split", split, "doesn't exist in the phenotype file. Excluded from prediction.\n"))
        split_name <- setdiff(split_name, split)
      }
    }
  }

  if ("train" %in% split_name) {
    split_name <- c("train", setdiff(split_name, "train"))
  }

  for (split in split_name) {
    ids[[split]] <- ids[[split]][!is.na(ids[[split]])]
  }

  phe <- list()
  for (split in split_name) {
    ids_loc <- match(ids[[split]], phe_master[["ID"]])
    phe[[split]] <- phe_master[ids_loc]
  }

  covariates <- list()

  for (split in split_name) {
    if (length(covariate_names) > 0) {
      covariates[[split]] <- phe[[split]][, covariate_names, with = FALSE]
    } else {
      covariates[[split]] <- NULL
    }
  }

  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(zstdcat_path, ' ', paste0(new_genotype_file, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  pvar <- pgenlibr::NewPvar(paste0(new_genotype_file, '.pvar.zst'))
  chr <- list()
  for (split in split_name) {
    chr[[split]] <- pgenlibr::NewPgen(paste0(new_genotype_file, '.pgen'), pvar = pvar, sample_subset = match(ids[[split]], ids[["psam"]]))
  }
  pgenlibr::ClosePvar(pvar)

  feature_names <- c()
  for (i in seq_along(fit)) {
    feature_names <- c(feature_names, which_row_active(fit[[i]]$B))
  }
  feature_names <- unique(feature_names)

  features <- list()
  for (split in split_name) {
    if (!is.null(covariates[[split]]) && is_full_rank) {
      features[[split]] <- data.table::data.table(covariates[[split]])
      features[[split]][, (feature_names) := snpnet:::prepareFeatures(chr[[split]], vars, feature_names, stats)]
    } else {
      features[[split]] <- snpnet:::prepareFeatures(chr[[split]], vars, feature_names, stats)
    }
  }

  pred <- list()
  R2 <- list()
  if (length(binary_phenotypes) > 0) AUC <- list()
  for (split in split_name) {
    pred[[split]] <- array(dim = c(nrow(features[[split]]), length(fit[[length(fit)]][["a0"]]), length(fit)),
                  dimnames = list(ids[[split]], phenotype_names, seq_along(fit)))
    R2[[split]] <- array(dim = c(length(fit[[length(fit)]][["a0"]]), length(fit)),
                         dimnames = list(phenotype_names, seq_along(fit)))
    if (length(binary_phenotypes) > 0) {
      AUC[[split]] <- array(dim = c(length(binary_phenotypes), length(fit)),
                            dimnames = list(binary_phenotypes, seq_along(fit)))
    }
  }

  response <- list()
  variance <- list()
  for (split in split_name) {
    response[[split]] <- as.matrix(phe[[split]][, phenotype_names, with = F])
    if (length(binary_phenotypes) > 0) {
      response[[split]][, binary_phenotypes] <- response[[split]][, binary_phenotypes] - 1
    }
    variance[[split]] <- apply(response[[split]], 2, function(x) mean((x - mean(x, na.rm = T))^2, na.rm = T))
  }

  for (i in seq_along(fit)) {
    if (length(binary_phenotypes) > 0 && "train" %in% split_name) glmfit_bin <- list()
    for (split in split_name) {
      if (!is.null(covariates[[split]]) && !is_full_rank) {
        active_vars <- which_row_active(fit[[i]]$C)
        if (length(active_vars) > 0) {
          features_single <- as.matrix(features[[split]][, active_vars, with = F])
        } else {
          features_single <- matrix(0, nrow = nrow(features[[split]]), ncol = 0)
        }
        pred_var <- safe_product(features_single, fit[[i]]$C[active_vars, , drop = F])
        if (is.null(fit[[i]]$W) || nrow(fit[[i]]$W) < ncol(covariates[[split]])) {
          pred_single <- sweep(pred_var, 2, fit[[i]]$a0, FUN = "+")
        } else {
          pred_single <- as.matrix(covariates[[split]]) %*% fit[[i]]$W + sweep(pred_var, 2, fit[[i]]$a0, FUN = "+")
        }
      } else {
        active_vars <- which_row_active(fit[[i]]$CC)
        if (length(active_vars) > 0) {
          features_single <- as.matrix(features[[split]][, active_vars, with = F])
        } else {
          feature_single <- matrix(0, nrow = nrow(features[[split]]), ncol = 0)
        }
        pred_var <- safe_product(features_single, fit[[i]]$CC[active_vars, , drop = F])
        pred_single <- sweep(pred_var, 2, fit[[i]]$a0, FUN = "+")
      }
      pred_single <- as.matrix(pred_single)
      colnames(pred_single) <- colnames(fit[[i]]$C)
      pred_single <- as.matrix(y_de_standardization(pred_single, std_obj$means, std_obj$sds, weight))
      R2[[split]][, i] <- 1 - apply((pred_single - response[[split]])^2, 2, mean, na.rm = T) / variance[[split]]
      if (length(binary_phenotypes) > 0) {
        for (bphe in binary_phenotypes) {
          if ("train" %in% split_name) {
            if (split == "train") {
              data_logistic_train <- data.frame(response = response[[split]][, bphe], covariates[[split]], score = pred_single[, bphe])
              glmfit_bin[[bphe]] <- glm(response ~ ., data = data_logistic_train, family = binomial())
            }
            data_logistic_split <- data.frame(covariates[[split]], score = pred_single[, bphe])
            pred_prob_split <- predict(glmfit_bin[[bphe]], newdata = data_logistic_split, type = "response")
            pred_single[, bphe] <- pred_prob_split
          }
          not_missing <- !is.na(response[[split]][, bphe])
          pred_obj <- ROCR::prediction(pred_single[not_missing, bphe], response[[split]][not_missing, bphe])
          auc_obj <- ROCR::performance(pred_obj, measure = 'auc')
          AUC[[split]][bphe, i] <- auc_obj@y.values[[1]]
        }
      }
      pred[[split]][, , i] <- pred_single
    }
  }

  out <- list(prediction = pred, response = response, R2 = R2)
  if (length(binary_phenotypes) > 0) out[["AUC"]] <- AUC

  out
}


#' Make plots of the multisnpnet results
#'
#' For 50th lambda, the reduced-rank results are saved at
#' results_dir[i]/${rank_prefix}[i]${rank}[j]/${file_prefix}50${file_suffix}, and snpnet results are
#' saved in ${snpnet_dir}/${phenotype}/${snpnet_subdir}/${snpnet_prefix}50${snpnet_suffix}
#'
#' @param results_dir Character vector each specifies the parent directory of one type (e.g. exact,
#'   lazy) of results
#' @param rank_prefix Character vector each specifies the prefix of the subdirectories holding
#'   per-rank results separately
#' @param type Character vector each specifies the type of the results, e.g., exact, lazy
#' @param rank Numeric vector of ranks for which the results are available
#' @param file_prefix Character vector each specifies the prefix of the result file
#' @param file_suffix Character vector each specifies the suffix of the result file
#' @param snpnet_dir Parent directory the snpnet results
#' @param snpnet_subdir Name of the subdirectory hosting snpnet results. The result files are saved
#'   in
#' @param save_dir Directory to save the plots. If NULL, no plots are generated but the list of plot
#'   objects is returned
#' @param train_name Name of the object storing the training metric
#' @param train_name Name of the object storing the validation metric
#' @param xlim The x limits (x1, x2) of the plots
#' @param ylim The y limits (y1, y2) of the plots
#'
#' @import ggplot2
#'
#' @export
plot_multisnpnet <- function(results_dir, rank_prefix, type, rank,
                             file_prefix, file_suffix,
                             snpnet_dir = NULL, snpnet_subdir = NULL, snpnet_prefix = NULL, snpnet_suffix = NULL,
                             save_dir = NULL, train_name = "metric_train", val_name = "metric_val", test_name = NULL, metric_name = "R2",
                             train_bin_name = "AUC_train", val_bin_name = "AUC_val", test_bin_name = "AUC_test", metric_bin_name = "AUC",
                             xlim = c(NA, NA), ylim = c(NA, NA), mapping_phenotype = NULL) {
  if (!is.null(save_dir)) dir.create(save_dir, recursive = T)
  data_metric_full <- NULL
  bin_names <- c()
  for (dir_idx in seq_along(results_dir)) {
    for (r in rank) {
      dir_rank <- file.path(results_dir[dir_idx], paste0(rank_prefix[dir_idx], r))
      files_in_dir <- list.files(dir_rank)
      result_files <- files_in_dir[startsWith(files_in_dir, file_prefix[dir_idx])]
      max_iter <- max(as.numeric(gsub(file_suffix[dir_idx], "", gsub(pattern = file_prefix[dir_idx], "", result_files))))
      latest_result <- file.path(dir_rank, paste0(file_prefix[dir_idx], max_iter, file_suffix[dir_idx]))

      myenv <- new.env()
      load(latest_result, envir = myenv)
      metric_train <- myenv[[train_name]]
      metric_val <- myenv[[val_name]]
      if (!is.null(test_name)) {
        if (!(test_name %in% names(myenv))) {
          stop("Test result doesn't exist for multisnpnet rank ", r, ".\n")
        }
        metric_test <- myenv[[test_name]]
      }
      if ((train_bin_name %in% names(myenv)) && (val_bin_name %in% names(myenv))) {
        AUC_train <- myenv[[train_bin_name]]
        AUC_val <- myenv[[val_bin_name]]
        if (!is.null(test_name)) AUC_test <- myenv[[test_bin_name]]
        bin_names <- unique(c(bin_names, colnames(AUC_train)))
        metric_train[, colnames(AUC_train)] <- AUC_train
        metric_val[, colnames(AUC_val)] <- AUC_val
        if (!is.null(test_name)) metric_test[, colnames(AUC_test)] <- AUC_test
      }
      imax_train <- max(which(apply(metric_train, 1, function(x) sum(is.na(x))) == 0))
      imax_val <- max(which(apply(metric_val, 1, function(x) sum(is.na(x))) == 0))
      imax <- min(imax_train, imax_val)
      if (!is.null(test_name)) {
        imax_test <- max(which(apply(metric_test, 1, function(x) sum(is.na(x))) == 0))
        imax <- min(imax, imax_test)
        metric_test <- metric_test[1:imax, , drop = F]
        metric_test <- cbind(metric_test, lambda = 1:imax)
        table_test <- reshape2::melt(as.data.frame(metric_test), id.vars = "lambda", variable.name = "phenotype", value.name = "metric_test")
      }
      metric_train <- metric_train[1:imax, , drop = F]
      metric_val <- metric_val[1:imax, , drop = F]
      metric_train <- cbind(metric_train, lambda = 1:imax)
      metric_val <- cbind(metric_val, lambda = 1:imax)

      table_train <- reshape2::melt(as.data.frame(metric_train), id.vars = "lambda", variable.name = "phenotype", value.name = "metric_train")
      table_val <- reshape2::melt(as.data.frame(metric_val), id.vars = "lambda", variable.name = "phenotype", value.name = "metric_val")
      data_metric <- dplyr::inner_join(table_train, table_val, by = c("phenotype", "lambda"))
      if (!is.null(test_name)) {
        data_metric <- dplyr::inner_join(data_metric, table_test, by = c("phenotype", "lambda"))
      }
      data_metric[["type"]] <- type[dir_idx]
      data_metric[["rank"]] <- factor(r, levels = as.character(rank))

      data_metric_full <- rbind(data_metric_full, data_metric)
    }
  }

  if (!is.null(snpnet_dir)) {
    for (phe in as.character(unique(data_metric_full[["phenotype"]]))) {
      print(phe)
      phe_dir <- file.path(snpnet_dir, phe, snpnet_subdir)  # results/results
      files_in_dir <- list.files(phe_dir)
      result_files <- files_in_dir[startsWith(files_in_dir, snpnet_prefix) & endsWith(files_in_dir, snpnet_suffix)]
      max_iter <- max(as.numeric(gsub(snpnet_suffix, "", gsub(pattern = snpnet_prefix, "", result_files))))
      latest_result <- file.path(phe_dir, paste0(snpnet_prefix, max_iter, snpnet_suffix))

      myenv <- new.env()
      load(latest_result, envir = myenv)
      metric_train <- myenv[["metric.train"]]
      metric_val <- myenv[["metric.val"]]
      if (!is.null(test_name)) {
        if (!("metric.test" %in% names(myenv))) {
          stop("Test result doesn't exist for ", phe, ".\n")
        }
        metric_test <- myenv[["metric.test"]]
        imax_test <- max(which(!is.na(metric_test)))
      }

      imax_train <- max(which(!is.na(metric_train)))
      imax_val <- max(which(!is.na(metric_val)))
      imax <- min(imax_train, imax_val)
      if (!is.null(test_name)) imax <- min(imax, imax_test)

      table_snpnet <- data.frame(lambda = 1:imax, phenotype = rep(phe, imax), metric_train = metric_train[1:imax],
                                 metric_val = metric_val[1:imax], type = "exact", rank = "snpnet")
      if (!is.null(test_name)) table_snpnet[["metric_test"]] <- metric_test[1:imax]

      data_metric_full <- rbind(data_metric_full, table_snpnet)
    }
  }

  if (!is.null(mapping_phenotype)) {
    data_metric_full$phenotype <- as.character(data_metric_full$phenotype)
    for (phe in names(mapping_phenotype)) {
      data_metric_full$phenotype[data_metric_full$phenotype == phe] <- mapping_phenotype[phe]
    }
    reverse_mapping <- names(mapping_phenotype)
    names(reverse_mapping) <- mapping_phenotype
  }

  gp <- list(data = data_metric_full)

  if (!is.null(snpnet_dir)) {
    max_metric_reduced_rank <- data_metric_full %>%
      dplyr::filter(rank != "snpnet") %>%
      dplyr::group_by(phenotype) %>%
      dplyr::filter(metric_val == max(metric_val)) %>%
      dplyr::mutate(max_multisnpnet_val = ifelse(!is.null(test_name), metric_test, metric_val)) %>%
      dplyr::select(phenotype, max_multisnpnet_val, rank, lambda)
    max_metric_snpnet <- data_metric_full %>%
      dplyr::filter(rank == "snpnet") %>%
      dplyr::group_by(phenotype) %>%
      dplyr::filter(metric_val == max(metric_val)) %>%
      dplyr::mutate(max_snpnet_val = ifelse(!is.null(test_name), metric_test, metric_val)) %>%
      dplyr::select(phenotype, max_snpnet_val)
    max_metric <- max_metric_reduced_rank %>%
      dplyr::inner_join(max_metric_snpnet, by = "phenotype") %>%
      dplyr::mutate(absolute_change = max_multisnpnet_val - max_snpnet_val,
                    relative_change = max_multisnpnet_val/abs(max_snpnet_val)-1,
                    direction = ifelse(relative_change > 0, "P", "N"))
    val_test_label <- ifelse(!is.null(test_name), "Test ", NULL)
    gp[["max_metric"]] <- max_metric
    gp[["metric_cmp_abs_change"]] <- ggplot(max_metric, aes(x = phenotype, y = absolute_change)) +
      geom_bar(stat = "identity", position = "dodge", aes(fill = direction)) +
      geom_hline(yintercept = 0, colour = "grey90") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      xlab("Phenotype") + ylab(paste0(val_test_label, "Metric Absolute Change"))
    if (!is.null(save_dir)) {
      save_path <- file.path(save_dir, "metric_cmp_abs_change.pdf")
      ggsave(save_path, plot = gp[["metric_cmp_abs_change"]])
    }
    gp[["metric_cmp_rel_change"]] <- ggplot(max_metric, aes(x = phenotype, y = relative_change*100)) +
      geom_bar(stat = "identity", position = "dodge", aes(fill = direction)) +
      geom_hline(yintercept = 0, colour = "grey90") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      xlab("Phenotype") + ylab(paste0(val_test_label, "Metric Relative Change (%)"))
    if (!is.null(save_dir)) {
      save_path <- file.path(save_dir, "metric_cmp_rel_change.pdf")
      ggsave(save_path, plot = gp[["metric_cmp_rel_change"]])
    }
    # relative plot with absoluate value on the second y axis
    abs_range <- range(max_metric$absolute_change)
    rel_range <- range(max_metric$relative_change)
    if (abs_range[1] * abs_range[2] < 0) {
      multiplier <- min(rel_range[2] / abs_range[2], rel_range[1] / abs_range[1]) * 100
    } else {
      multiplier <- max(abs(rel_range)) / max(abs(abs_range)) * 100 / 2
    }
    gp[["metric_cmp_abs_rel_change"]] <- ggplot(max_metric, aes(x = reorder(phenotype, -relative_change), y = relative_change*100)) +
      geom_bar(stat = "identity", position = "dodge", aes(fill = direction)) +
      geom_hline(yintercept = 0, colour = "grey90") +
      geom_point(aes(y = absolute_change * multiplier), size = 1.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      scale_y_continuous(sec.axis = sec_axis(~. * (1.0/multiplier), name = paste0(val_test_label, "Metric Absolute Change"))) +
      xlab("Phenotype") + ylab(paste0(val_test_label, "Metric Relative Change (%)"))
    if (!is.null(save_dir)) {
      save_path <- file.path(save_dir, "metric_cmp_abs_rel_change.pdf")
      ggsave(save_path, plot = gp[["metric_cmp_abs_rel_change"]])
    }
  }

  for (phe in as.character(unique(data_metric_full[["phenotype"]]))) {
    fname_phe <- ifelse(!is.null(mapping_phenotype) && (phe %in% mapping_phenotype), reverse_mapping[phe], phe)
    mname <- ifelse(fname_phe %in% bin_names, metric_bin_name, metric_name)
    gp[[fname_phe]] <- ggplot(dplyr::filter(data_metric_full, phenotype == phe), aes(x = metric_train, y = metric_val, shape = type, colour = rank)) +
      geom_path() + geom_point() +
      xlab(paste(mname, "(train)")) + ylab(paste(mname, "(val)")) +
      xlim(as.numeric(xlim)) + ylim(as.numeric(ylim)) +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
            legend.text=element_text(size=12), legend.title = element_text(size=12),
            legend.position = "bottom",
            strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +
      ggtitle(phe)
    if (!is.null(save_dir)) {
      save_path <- file.path(save_dir, paste0(mname, "_plot_", fname_phe, ".pdf"))
      ggsave(save_path, plot = gp[[fname_phe]])
    }
  }

  gp
}

safe_product <- function(X, Y, MAXLEN = (2^31 - 1) / 2) {
  ncol.chunk <- floor(MAXLEN / as.double(nrow(X)))  # depends on the memory requirements
  numChunks <- ceiling(ncol(X) / as.double(ncol.chunk))
  out <- matrix(0, nrow(X), ncol(Y))
  rownames(out) <- rownames(X)
  colnames(out) <- colnames(Y)
  for (jc in seq_len(numChunks)) {
    idx <- ((jc-1)*ncol.chunk+1):min(jc*ncol.chunk, ncol(X))
    out <- out + X[, idx, drop=FALSE] %*% Y[idx, , drop=FALSE]
  }
  out
}

#' Make biplots of the multisnpnet results
#'
#' Generate biplot visualization based on the decomposed coefficient matrix C.
#' One of the most common use case is: plot_biplot(svd(t(fit$C)), label=list('phenotype'=rownames(A_init), 'variant'=rownames(fit$C)))
#'
#' @param svd_obj A named list containing three matrices with u, d, and v as their names as in the
#'   output from base::svd() function. One can pass the results of base::svd(t(fit$C)).
#'   Please note that this function assumes svd_obj$u and svd_obj$v corresponds to phenotypes and variants, respectively.
#' @param component A named list that specifies the index of the components used in the plot.
#' @param label A named list that specifies the phenotype and variant labels.
#'   The labels needs to be the same order as in svd_obj$u and svd_obj$v.
#' @param n_labels A named list that specifies the number of phenotype and variant labels in the plot.
#' @param color A named list that specifies the color in the plot.
#' @param shape A named list that specifies the color in the plot.
#' @param axis_label A named list that specifies the names used in the axis labels.
#' @param use_ggrepel A binary variable that specifies whether we should use ggrepel to annotate the
#'   labels of the data points.
#'
#' @import ggplot2
#' @importFrom magrittr '%>%'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename select mutate if_else bind_rows
#'
#' @export
plot_biplot <- function(svd_obj, component=list('x'=1, 'y'=2),
                        label=list('phenotype'=NULL, 'variant'=NULL),
                        n_labels=list('phenotype'=5, 'variant'=5),
                        color=list('phenotype'='red', 'variant'='blue'),
                        shape=list('phenotype'=20, 'variant'=4),
                        axis_label=list('main'='variant', 'sub'='phenotype'),
                        use_ggrepel=TRUE) {
    # extract the relevant matrices from the svd object
    u  <- svd_obj$u
    vd <- (svd_obj$v) %*% (diag(svd_obj$d))

    # assign row and col names
    if(is.null(label[['phenotype']])){ label[['phenotype']] <- paste0('phenotype', 1:nrow(u)) }
    if(is.null(label[['variant']])){   label[['variant']]   <- paste0('variant',   1:nrow(vd)) }
    rownames(u)  <- label[['phenotype']]
    rownames(vd) <- label[['variant']]
    colnames(u)  <- 1:length(svd_obj$d)
    colnames(vd) <- 1:length(svd_obj$d)

    # convert the matrices into data frames
    df_u  <- u  %>% as.data.frame() %>% rename('PC_x' := component$x, 'PC_y' := component$y) %>%
    select(PC_x, PC_y) %>% rownames_to_column('label') %>%
    mutate(label = if_else(rank(-(PC_x**2+PC_y**2))<=n_labels[['phenotype']], label, ''))

    df_vd <- vd %>% as.data.frame() %>% rename('PC_x' := component$x, 'PC_y' := component$y) %>%
    select(PC_x, PC_y) %>% rownames_to_column('label') %>%
    mutate(label = if_else(rank(-(PC_x**2+PC_y**2))<=n_labels[['variant']], label, ''))

    # scale u (data on sub-axis) to map to the main-axis
    lim_u_abs   <- 1.1 * max(abs(df_u  %>% select(PC_x, PC_y)))
    lim_vd_abs  <- 1.1 * max(abs(df_vd %>% select(PC_x, PC_y)))

    df_u_scaled <- df_u %>%
    mutate(
        PC_x = PC_x * (lim_vd_abs/lim_u_abs),
        PC_y = PC_y * (lim_vd_abs/lim_u_abs)
    )

    if(! use_ggrepel){
      # generate plot without ggrepel. This is useful when you'd like
      # to conver the plot into plotly object using gplotly.
        p <- ggplot() +
        layer(
          # scatter plot for (VD)
            data=df_vd, mapping=aes(x=PC_x, y=PC_y, shape='variant', label=label),
            geom='point', stat = "identity", position = "identity",
            params=list(size=1, color=color[['variant']])
        )+
        layer(
          # segments (lines) for U
            data=df_u_scaled,
            mapping=aes(x=0, y=0, xend=PC_x, yend=PC_y),
            geom='segment', stat = "identity", position = "identity",
            params=list(size=1, color=color[['phenotype']], alpha=.2)
        )+
        layer(
          # scatter plot for U
            data=df_u_scaled,
            mapping=aes(x=PC_x, y=PC_y, shape='phenotype', label=label),
            geom='point', stat = "identity", position = "identity",
            params=list(size=1, color=color[['phenotype']])
        )

    } else { # use_ggrepel == TRUE
      # generate the plot with ggrepel
        p <- ggplot() +
        layer(
            data=df_vd, mapping=aes(x=PC_x, y=PC_y, shape='variant'),
            geom='point', stat = "identity", position = "identity",
            params=list(size=1, color=color[['variant']])
        )+
        layer(
            data=df_u_scaled,
            mapping=aes(x=0, y=0, xend=PC_x, yend=PC_y),
            geom='segment', stat = "identity", position = "identity",
            params=list(size=1, color=color[['phenotype']], alpha=.2)
        )+
        layer(
            data=df_u_scaled,
            mapping=aes(x=PC_x, y=PC_y, shape='phenotype'),
            geom='point', stat = "identity", position = "identity",
            params=list(size=1, color=color[['phenotype']])
        )+
        ggrepel::geom_text_repel(
            data=bind_rows(
                df_u_scaled %>% mutate(color=color[['phenotype']]),
                df_vd       %>% mutate(color=color[['variant']])
            ),
            mapping=aes(x=PC_x, y=PC_y, label=label, color=color),
            size=3, force=10
        )
    }

    # configure the theme, axis, and axis labels
    p + theme_bw() +
    scale_color_manual(values=setNames(color, color)) +
    scale_shape_manual(values=shape) +
    guides(shape=FALSE,color=FALSE) +
    scale_x_continuous(
        sprintf('Component %s (%s [%s])', component$x, axis_label[['main']], color[['variant']]),
        limits = c(-lim_vd_abs, lim_vd_abs),
        sec.axis = sec_axis(
            ~ . * (lim_u_abs/lim_vd_abs),
            name = sprintf('Component %s (%s [%s])', component$y, axis_label[['sub']], color[['phenotype']])
        )
    ) +
    scale_y_continuous(
        sprintf('Component %s (%s [%s])', component$x, axis_label[['main']], color[['variant']]),
        limits = c(-lim_vd_abs, lim_vd_abs),
        sec.axis = sec_axis(
            ~ . * (lim_u_abs/lim_vd_abs),
            name = sprintf('Component %s (%s [%s])', component$y, axis_label[['sub']], color[['phenotype']])
        )
    )
}
