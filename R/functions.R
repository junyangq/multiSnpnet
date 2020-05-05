setupMultiConfigs <- function(configs, standardize_response, max.iter, rank, prev_iter, batch_size) {
  out.args <- as.list(environment())
  defaults_multi <- list(
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
#' @param split_name Label for the samples where prediction is to be made.
#' @param zstdcat_path Path to zstdcat program, needed when loading variants
#'
#' @return A list containing the prediction and the resopnse for which the prediction is made.
#'
#' @export
predict_multisnpnet <- function(fit = NULL, saved_path = NULL, new_genotype_file, new_phenotype_file,
                                idx = NULL, covariate_names = NULL, split_col = NULL, split_name = "test",
                                zstd_path = "zstdcat") {
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
  phenotype_names <- colnames(fit[[length(fit)]][["C"]])
  is_full_rank <- (ncol(fit[[length(fit)]][["B"]]) == ncol(fit[[length(fit)]][["C"]]))

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

  if (is.null(split_col)) {
    ids[[split_name]] <- phe_master$ID
  } else {
    ids[[split_name]] <- phe_master$ID[phe_master[[split_col]] == split_name]
  }

  phe_test <- phe_master[match(ids[[split_name]], phe_master[["ID"]])]

  if (length(covariate_names) > 0) {
    covariates <- phe_test[, covariate_names, with = FALSE]
  } else {
    covariates <- NULL
  }

  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(zstdcat_path, ' ', paste0(new_genotype_file, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  pvar <- pgenlibr::NewPvar(paste0(new_genotype_file, '.pvar.zst'))
  chr <- pgenlibr::NewPgen(paste0(new_genotype_file, '.pgen'), pvar = pvar, sample_subset = match(ids[[split_name]], ids[["psam"]]))
  pgenlibr::ClosePvar(pvar)

  feature_names <- c()
  for (i in seq_along(fit)) {
    feature_names <- c(feature_names, rownames(fit[[i]]$B))
  }
  feature_names <- unique(feature_names)

  if (!is.null(covariates) && is_full_rank) {
    features <- data.table(covariates)
    features[, (feature_names) := snpnet:::prepareFeatures(chr, vars, feature_names, stats)]
  } else {
    features <- snpnet:::prepareFeatures(chr, vars, feature_names, stats)
  }

  pred <- array(dim = c(nrow(features), length(fit[[length(fit)]][["a0"]]), length(fit)),
                dimnames = list(ids[[split_name]], phenotype_names, seq_along(fit)))

  for (i in seq_along(fit)) {
    if (!is.null(covariates) && !is_full_rank) {
      features_single <- as.matrix(features[, rownames(fit[[i]]$C), with = F])
      pred_single <- as.matrix(covariates) %*% fit[[i]]$W + sweep(features_single %*% fit[[i]]$C, 2, fit[[i]]$a0, FUN = "+")
    } else {
      features_single <- as.matrix(features[, rownames(fit[[i]]$CC), with = F])
      pred_single <- sweep(features_single %*% fit[[i]]$CC, 2, fit[[i]]$a0, FUN = "+")
    }
    pred_single <- multisnpnet:::y_de_standardization(pred_single, std_obj$means, std_obj$sds, weight)
    pred[, , i] <- as.matrix(pred_single)
  }

  list(prediction = pred, response = as.matrix(phe_test[, phenotype_names, with = F]))
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
                             save_dir = NULL, train_name = "metric_train", val_name = "metric_val",
                             xlim = c(NA, NA), ylim = c(NA, NA)) {
  if (!is.null(save_dir)) dir.create(save_dir)
  data_metric_full <- NULL
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
      imax_train <- max(which(apply(metric_train, 1, function(x) sum(is.na(x))) == 0))
      imax_val <- max(which(apply(metric_val, 1, function(x) sum(is.na(x))) == 0))
      imax <- min(imax_train, imax_val)
      metric_train <- metric_train[1:imax, , drop = F]
      metric_val <- metric_val[1:imax, , drop = F]
      metric_train <- cbind(metric_train, lambda = 1:imax)
      metric_val <- cbind(metric_val, lambda = 1:imax)

      table_train <- reshape2::melt(as.data.frame(metric_train), id.vars = "lambda", variable.name = "phenotype", value.name = "metric_train")
      table_val <- reshape2::melt(as.data.frame(metric_val), id.vars = "lambda", variable.name = "phenotype", value.name = "metric_val")
      data_metric <- dplyr::inner_join(table_train, table_val, by = c("phenotype", "lambda"))
      data_metric[["type"]] <- type[dir_idx]
      data_metric[["rank"]] <- factor(r, levels = as.character(rank))

      data_metric_full <- rbind(data_metric_full, data_metric)
    }
  }

  if (!is.null(snpnet_dir)) {
    for (phe in as.character(unique(data_metric[["phenotype"]]))) {
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

      imax_train <- max(which(!is.na(metric_train)))
      imax_val <- max(which(!is.na(metric_val)))
      imax <- min(imax_train, imax_val)

      table_snpnet <- data.frame(lambda = 1:imax, phenotype = rep(phe, imax), metric_train = metric_train[1:imax],
                                 metric_val = metric_val[1:imax], type = "exact", rank = "snpnet")

      data_metric_full <- rbind(data_metric_full, table_snpnet)
    }
  }

  gp <- list()

  if (!is.null(snpnet_dir)) {
    max_metric_reduced_rank <- data_metric_full %>%
      dplyr::filter(rank != "snpnet") %>%
      dplyr::group_by(phenotype) %>%
      dplyr::summarise(max_multisnpnet_val = max(metric_val))
    max_metric_snpnet <- data_metric_full %>%
      dplyr::filter(rank == "snpnet") %>%
      dplyr::group_by(phenotype) %>%
      dplyr::summarise(max_snpnet_val = max(metric_val))
    max_metric <- max_metric_reduced_rank %>%
      dplyr::inner_join(max_metric_snpnet, by = "phenotype") %>%
      dplyr::mutate(absolute_change = max_multisnpnet_val - max_snpnet_val,
                    relative_change = max_multisnpnet_val/abs(max_snpnet_val)-1,
                    direction = ifelse(relative_change > 0, "P", "N"))
    gp[["r2_cmp_rel_change"]] <- ggplot(max_metric, aes(x = phenotype, y = relative_change*100)) +
      geom_bar(stat = "identity", position = "dodge", aes(fill = direction)) +
      geom_hline(yintercept = 0, colour = "grey90") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      xlab("phenotypes") + ylab("R2 Relative Change (%)")
    if (!is.null(save_dir)) {
      save_path <- file.path(save_dir, "R2_cmp_rel_change.pdf")
      ggsave(save_path, plot = gp[["r2_cmp_rel_change"]])
    }
  }

  for (phe in as.character(unique(data_metric[["phenotype"]]))) {
    gp[[phe]] <- ggplot(dplyr::filter(data_metric_full, phenotype == phe), aes(x = metric_train, y = metric_val, shape = type, colour = rank)) +
      geom_path() + geom_point() +
      xlab("metric (train)") + ylab("metric (val)") +
      xlim(as.numeric(xlim)) + ylim(as.numeric(ylim)) +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
            legend.text=element_text(size=12), legend.title = element_text(size=12),
            legend.position = "bottom",
            strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +
      ggtitle(phe)
    if (!is.null(save_dir)) {
      save_path <- file.path(save_dir, paste0("metric_plot_", phe, ".pdf"))
      ggsave(save_path, plot = gp[[phe]])
    }
  }
  gp
}
