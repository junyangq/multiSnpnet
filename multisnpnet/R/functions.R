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

  colnames(C) <- colnames(Y)

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

setup_configs_directories <- function(configs, covariates, save, results.dir, validation, nlambda, lambda.min.ratio, standardize_response) {
  if (configs[["use_plink2"]] && !("mem" %in% names(configs)))
    stop("mem should be provided to guide the memory capacity for PLINK2.")
  configs[["covariates"]] <- covariates
  configs[["nlambda"]] <- nlambda
  configs[["lambda.min.ratio"]] <- lambda.min.ratio
  configs[["standardize_response"]] <- standardize_response
  configs[['gcount.basename.prefix']] <- "snpnet.train"
  default_settings <- list(missing.rate = 0.1, MAF.thresh = 0.001, nCores = 1,
                           nlams.init = 10, nlams.delta = 5,
                           gcount.full.prefix = NULL, vzs = TRUE)
  for (name in names(default_settings)) {
    if (!(name %in% names(configs))) configs[[name]] <- default_settings[[name]]
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
  configs[["gcount.full.prefix"]] <- file.path(configs[['results.dir']], configs[['meta.dir']], configs[['gcount.basename.prefix']])
  configs[["zstdcat.path"]] <- 'zstdcat'
  configs[["save.computeProduct"]] <- FALSE
  configs[["parent.dir"]] <- results.dir
  configs[["endian"]] <- "little"
  configs[["save"]] <- save
  configs
}

###### From Yosuke's implementation using PLINK 2.0 for inner product ######
timeDiff <- function(start.time, end.time = NULL) {
  if (is.null(end.time)) end.time <- Sys.time()
  paste(round(end.time-start.time, 4), units(end.time-start.time))
}

computeStats_P2 <- function(pfile, ids, configs) {
  time.computeStats.start <- Sys.time()
  snpnetLogger('Start computeStats()', indent=2, log.time=time.computeStats.start)
  keep_f       <- paste0(configs[['gcount.full.prefix']], '.keep')
  gcount_tsv_f <- paste0(configs[['gcount.full.prefix']], '.gcount.tsv')

  dir.create(dirname(configs[['gcount.full.prefix']]), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(gcount_tsv_f)) {
    gcount_df <- fread(gcount_tsv_f)
  } else {
    # To run plink2 --geno-counts, we write the list of IDs to a file
    data.frame(ID = ids) %>%
      separate(ID, into=c('FID', 'IID'), sep='_') %>%
      fwrite(keep_f, sep='\t', col.names=F)

    # Run plink2 --geno-counts
    system(paste(
      'plink2',
      '--threads', configs[['nCores']],
      '--memory', configs[['mem']],
      '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
      '--keep', keep_f,
      '--out', configs[['gcount.full.prefix']],
      '--geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs',
      sep=' '
    ), intern=F, wait=T)

    # read the gcount file
    gcount_df <-
      data.table::fread(paste0(configs[['gcount.full.prefix']], '.gcount')) %>%
      rename(original_ID = ID) %>%
      mutate(
        ID = paste0(original_ID, '_', ALT),
        stats_pNAs  = MISSING_CT / (MISSING_CT + OBS_CT),
        stats_means = (HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS ) / OBS_CT,
        stats_msts  = (HAP_ALT_CTS + HET_REF_ALT_CTS + 4 * TWO_ALT_GENO_CTS ) / OBS_CT,
        stats_SDs   = stats_msts - stats_means * stats_means
      )
  }

  out <- list()
  out[["pnas"]]  <- gcount_df %>% select(stats_pNAs) %>% pull()
  out[["means"]] <- gcount_df %>% select(stats_means) %>% pull()
  out[["sds"]]   <- gcount_df %>% select(stats_SDs) %>% pull()

  for(key in names(out)){
    names(out[[key]]) <- gcount_df %>% select(ID) %>% pull()
  }
  out[["excludeSNP"]] <- names(out[["means"]])[(out[["pnas"]] > configs[["missing.rate"]]) | (out[["means"]] < 2 * configs[["MAF.thresh"]])]

  if (configs[['save']]){
    gcount_df %>% fwrite(gcount_tsv_f, sep='\t')
    saveRDS(out[["excludeSNP"]], file = file.path(dirname(configs[['gcount.full.prefix']]), "excludeSNP.rda"))
  }
  snpnetLoggerTimeDiff('End computeStats().', time.computeStats.start, indent=3)
  out
}

computeProduct_P2 <- function(residual, pfile, vars, stats, configs, iter) {
  time.computeProduct.start <- Sys.time()
  snpnetLogger('Start computeProduct()', indent=2, log.time=time.computeProduct.start)

  gc_res <- gc()
  if(configs[['KKT.verbose']]) print(gc_res)

  snpnetLogger('Start plink2 --variant-score', indent=3, log.time=time.computeProduct.start)
  dir.create(file.path(configs[['parent.dir']], configs[["results.dir"]]), showWarnings = FALSE, recursive = T)

  residual_f <- file.path(configs[["parent.dir"]], configs[["results.dir"]], paste0("residuals_iter_", iter, ".tsv"))

  # write residuals to a file
  residual_df <- data.frame(residual)
  # colnames(residual_df) <- paste0('lambda_idx_', colnames(residual))
  residual_df %>%
    rownames_to_column("ID") %>%
    separate(ID, into=c('#FID', 'IID'), sep='_') %>%
    fwrite(residual_f, sep='\t', col.names=T)

  # Run plink2 --geno-counts
  system(paste(
    'plink2',
    '--threads', configs[['nCores']],
    '--memory', as.integer(configs[['mem']]) - ceiling(sum(as.matrix(gc_res)[,2])),
    '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
    '--read-freq', paste0(configs[['gcount.full.prefix']], '.gcount'),
    '--keep', residual_f,
    '--out', str_replace_all(residual_f, '.tsv$', ''),
    '--variant-score', residual_f, 'zs', 'bin',
    sep=' '
  ), intern=F, wait=T)

  prod.full <- readBinMat(str_replace_all(residual_f, '.tsv$', '.vscore'), configs)
  if (! configs[['save']]) system(paste(
    'rm', residual_f, str_replace_all(residual_f, '.tsv$', '.log'), sep=' '
  ), intern=F, wait=T)

  snpnetLoggerTimeDiff('End plink2 --variant-score.', time.computeProduct.start, indent=4)

  rownames(prod.full) <- vars
  if (configs[["standardize.variant"]]) {
    for(residual.col in 1:ncol(residual)){
      prod.full[, residual.col] <- apply(prod.full[, residual.col], 2, "/", stats[["sds"]])
    }
  }
  prod.full[stats[["excludeSNP"]], ] <- NA
  snpnetLoggerTimeDiff('End computeProduct().', time.computeProduct.start, indent=2)
  prod.full
}
###### ------------------------------------------------------------------- ######

# used to compute the new lambda.min.ratio if we want to extend the original lambda sequence
compute_lambda_min_ratio <- function(nlambda.new, nlambda = 100, ratio = 0.01) {
  exp((nlambda.new-1)/(nlambda-1)*log(ratio))
}
