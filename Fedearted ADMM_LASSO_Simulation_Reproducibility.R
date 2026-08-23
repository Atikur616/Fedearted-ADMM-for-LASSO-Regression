# ============================================================================
# Federated ADMM-LASSO
#simulation setting: p = 50, n = 100 per center
# ============================================================================
rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1)

# # Package and numerical settings
if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("The glmnet package is required. Install it with install.packages('glmnet').")
}

suppressPackageStartupMessages(library(glmnet))

# Force glmnet to return the complete user-supplied lambda path.
glmnet::glmnet.control(fdev = 0, devmax = 1)

# Keep numerical libraries single-threaded for more reproducible desktop runs.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)


### Simulation settings
# number of predictors
p <- 50L
##sample per center
n_per_center <- 100L
###number of centers
m <- 2L
## total replications
total_rep <- 100L

# ADMM penalty parameter rho
primary_rho <- 1
run_rho_sensitivity <- FALSE
#### for checking the rho sensitivity analysis
rho_values <- if (run_rho_sensitivity) {
  c(0.5, 0.75, 1:9)
} else {
  primary_rho
}

#### Regularization parameter lambda tuning based on the validation MinMSE
selection_rules <- "MinMSE"
primary_rule <- "MinMSE"

#### True sparse beta vector for ADMM consensus
active_beta_values <- c(5, -2, 3, 1.5, 7)
if (p < length(active_beta_values)) {
  stop("p must be at least 5.")
}
true_beta <- c(active_beta_values, rep(0, p - length(active_beta_values)))

# Homogeneous residual noise
noise_sd <- 1

# Common lambda-grid
n_lambda <- 50L
lambda_min_ratio <- 1e-3
lambda_fraction_grid <- exp(
  seq(log(1), log(lambda_min_ratio), length.out = n_lambda)
)

selection_threshold <- 1e-3

# ADMM convergence criteria
eps_abs <- 1e-4
eps_rel <- 1e-3
max_admm_iter <- 10000L

# seed
base_seed <- 1000L
RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

setting_label <- paste0("p = ", p, ", n = ", n_per_center)
setting_tag <- paste0("p", p, "_n", n_per_center)


cat("\n============================================================\n")
cat("FEDERATED ADMM-LASSO DESKTOP REPRODUCIBILITY RUN\n")
cat("============================================================\n")
cat("Setting            :", setting_label, "\n")
cat("Centers            :", m, "\n")
cat("Replications       :", total_rep, "\n")
cat("Rho values         :", paste(rho_values, collapse = ", "), "\n")
cat("Lambda rule        :", primary_rule, "\n")

# Soft threshold function
soft_threshold <- function(x, threshold) {
  sign(x) * pmax(abs(x) - threshold, 0)
}
#### data generation
simulate_data <- function(n, p, beta, sigma = 1) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- as.vector(X %*% beta + rnorm(n, mean = 0, sd = sigma))
  list(X = X, y = y)
}

stack_clients <- function(clients) {
  list(
    X = do.call(rbind, lapply(clients, function(x) x$X)),
    y = unlist(lapply(clients, function(x) x$y), use.names = FALSE)
  )
}
#### Calculate the selection matrices
calculate_selection_metrics <- function(beta_estimate,
                                        true_beta,
                                        threshold = 1e-3) {
  estimated_active <- abs(beta_estimate) > threshold
  true_active <- abs(true_beta) > threshold

  TP <- sum(estimated_active & true_active)
  FP <- sum(estimated_active & !true_active)
  FN <- sum(!estimated_active & true_active)
  TN <- sum(!estimated_active & !true_active)

  sensitivity <- if ((TP + FN) > 0L) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0L) TN / (TN + FP) else NA_real_

  list(
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Number_Selected = sum(estimated_active)
  )
}

loss_stats_matrix <- function(X, y, beta_estimate) {
  prediction <- as.vector(X %*% beta_estimate)
  error <- y - prediction
  squared_error <- error^2

  list(
    MSE = mean(squared_error),
    SE = if (length(squared_error) > 1L) {
      sd(squared_error) / sqrt(length(squared_error))
    } else {
      NA_real_
    },
    MAE = mean(abs(error))
  )
}

# Compute local loss summaries without transmitting individual-level records.
local_loss_summary <- function(X, y, beta_estimate) {
  prediction <- as.vector(X %*% beta_estimate)
  error <- y - prediction
  squared_error <- error^2

  list(
    N = length(error),
    Sum_Squared_Error = sum(squared_error),
    Sum_Squared_Loss_Squared = sum(squared_error^2),
    Sum_Absolute_Error = sum(abs(error))
  )
}
### aggregate the local losses
aggregate_loss_summaries <- function(local_summaries) {
  total_n <- sum(vapply(local_summaries, function(x) x$N, numeric(1)))
  total_sse <- sum(vapply(
    local_summaries,
    function(x) x$Sum_Squared_Error,
    numeric(1)
  ))
  total_squared_loss_squared <- sum(vapply(
    local_summaries,
    function(x) x$Sum_Squared_Loss_Squared,
    numeric(1)
  ))
  total_sae <- sum(vapply(
    local_summaries,
    function(x) x$Sum_Absolute_Error,
    numeric(1)
  ))

  if (!is.finite(total_n) || total_n <= 0) {
    stop("Federated loss aggregation received no observations.")
  }

  mse <- total_sse / total_n
  mae <- total_sae / total_n

  # Exact pooled standard error of the observation-level squared losses,
    validation_se <- if (total_n > 1L) {
    variance_numerator <- total_squared_loss_squared -
      (total_sse^2 / total_n)

    # Protect against tiny negative values caused by floating-point rounding.
    sample_variance <- max(
      variance_numerator / (total_n - 1L),
      0
    )
    sqrt(sample_variance / total_n)
  } else {
    NA_real_
  }

  list(
    MSE = mse,
    SE = validation_se,
    MAE = mae,
    N = total_n,
    SSE = total_sse,
    SAE = total_sae
  )
}

federated_loss_stats <- function(clients, beta_estimate) {
  local_summaries <- lapply(clients, function(client) {
    local_loss_summary(
      X = client$X,
      y = client$y,
      beta_estimate = beta_estimate
    )
  })

  aggregate_loss_summaries(local_summaries)
}

assert_federated_loss_equivalence <- function(clients,
                                               stacked_data,
                                               beta_estimate,
                                               tolerance = 1e-10) {
  federated_loss <- federated_loss_stats(
    clients = clients,
    beta_estimate = beta_estimate
  )

  stacked_loss <- loss_stats_matrix(
    X = stacked_data$X,
    y = stacked_data$y,
    beta_estimate = beta_estimate
  )

  differences <- c(
    MSE = abs(federated_loss$MSE - stacked_loss$MSE),
    SE = abs(federated_loss$SE - stacked_loss$SE),
    MAE = abs(federated_loss$MAE - stacked_loss$MAE)
  )

  comparison_scale <- max(
    1,
    abs(stacked_loss$MSE),
    abs(stacked_loss$SE),
    abs(stacked_loss$MAE)
  )

  if (any(!is.finite(differences)) ||
      max(differences) > tolerance * comparison_scale) {
    stop(
      "Federated aggregate-loss calculation did not match the stacked ",
      "calculation. Differences: ",
      paste(names(differences), signif(differences, 6), collapse = "; ")
    )
  }

  invisible(TRUE)
}

loss_path_stats <- function(X, y, beta_matrix) {
  prediction_matrix <- X %*% beta_matrix
  error_matrix <- sweep(prediction_matrix, 1L, y, FUN = "-")
  squared_error_matrix <- error_matrix^2

  list(
    MSE = colMeans(squared_error_matrix),
    SE = apply(
      squared_error_matrix,
      2L,
      function(x) if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_
    ),
    MAE = colMeans(abs(error_matrix))
  )
}

make_local_lambda_grid <- function(X, y, lambda_fraction_grid) {
  lambda_max <- max(abs(as.vector(crossprod(X, y))))
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    stop("Could not compute a positive finite local lambda_max.")
  }

  list(
    lambda_max = lambda_max,
    lambda_grid = lambda_max * lambda_fraction_grid
  )
}

make_global_lambda_grid <- function(clients_train, lambda_fraction_grid) {
  # Each center computes X_i^T y_i locally. The coordinator receives the
  # center-level vectors (or their secure aggregate) and uses only their sum.
  local_Xty_messages <- lapply(clients_train, function(client) {
    as.vector(crossprod(client$X, client$y))
  })

  global_Xty <- Reduce(`+`, local_Xty_messages)

  lambda_max <- max(abs(global_Xty))
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    stop("Could not compute a positive finite global lambda_max.")
  }

  list(
    lambda_max = lambda_max,
    lambda_grid = lambda_max * lambda_fraction_grid
  )
}

choose_lambda_indices <- function(tuning_df) {
  eligible <- which(
    tuning_df$Converged &
      is.finite(tuning_df$Validation_MSE) &
      is.finite(tuning_df$Validation_SE)
  )

  if (length(eligible) == 0L) {
    return(list(
      MinMSE = NA_integer_,
      OneSE = NA_integer_,
      OneSE_Cutoff = NA_real_
    ))
  }

  min_index <- eligible[which.min(tuning_df$Validation_MSE[eligible])]
  one_se_cutoff <- tuning_df$Validation_MSE[min_index] +
    tuning_df$Validation_SE[min_index]

  one_se_eligible <- eligible[
    tuning_df$Validation_MSE[eligible] <= one_se_cutoff
  ]

  one_se_index <- one_se_eligible[
    which.max(tuning_df$Lambda[one_se_eligible])
  ]

  list(
    MinMSE = min_index,
    OneSE = one_se_index,
    OneSE_Cutoff = one_se_cutoff
  )
}

find_lambda_index <- function(lambda_vector,
                              target_lambda,
                              tolerance = 1e-8) {
  index <- which.min(abs(lambda_vector - target_lambda))
  difference <- abs(lambda_vector[index] - target_lambda)
  scale <- max(1, abs(target_lambda))

  if (!is.finite(difference) || difference > tolerance * scale) {
    stop(
      "Target lambda was not found on the pooled numerical lambda grid. ",
      "Minimum absolute difference = ", difference
    )
  }

  index
}

bind_rows_safe <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) return(data.frame())
  do.call(rbind, x)
}

#Federated ADMM-LASSO functions
prepare_admm_systems <- function(clients, rho, p) {
  lapply(clients, function(client) {
    XtX <- crossprod(client$X)
    Xty <- as.vector(crossprod(client$X, client$y))

    chol_factor <- tryCatch(
      chol(XtX + rho * diag(p)),
      error = function(e) {
        stop(
          "Cholesky factorization failed for rho = ", rho,
          ". Original message: ", conditionMessage(e)
        )
      }
    )

    list(Xty = Xty, chol_factor = chol_factor)
  })
}

solve_from_cholesky <- function(chol_factor, rhs) {
  backsolve(
    chol_factor,
    forwardsolve(t(chol_factor), rhs)
  )
}
### Local updated
admm_local_update <- function(prepared_system, z, u, rho) {
  rhs <- prepared_system$Xty + rho * (z - u)
  as.vector(solve_from_cholesky(prepared_system$chol_factor, rhs))
}
### run local ADMM
run_admm_lasso <- function(prepared_systems,
                           lambda,
                           rho,
                           p,
                           m,
                           eps_abs,
                           eps_rel,
                           max_iter,
                           z_init = NULL,
                           u_init = NULL,
                           save_history = FALSE) {
  z <- if (is.null(z_init)) rep(0, p) else as.vector(z_init)

  u_list <- if (is.null(u_init)) {
    lapply(seq_len(m), function(k) rep(0, p))
  } else {
    lapply(u_init, as.vector)
  }

  history_list <- if (save_history) vector("list", max_iter) else NULL

  converged <- FALSE
  r_norm <- NA_real_
  s_norm <- NA_real_
  eps_primal <- NA_real_
  eps_dual <- NA_real_

  for (iteration in seq_len(max_iter)) {
    beta_list <- lapply(seq_len(m), function(k) {
      admm_local_update(
        prepared_system = prepared_systems[[k]],
        z = z,
        u = u_list[[k]],
        rho = rho
      )
    })

    z_previous <- z

    # Each center sends only the composite update beta_i + u_i.
    center_update_messages <- Map(`+`, beta_list, u_list)

    # The coordinator averages the received composite vectors.
    average_beta_plus_u <- Reduce(
      `+`,
      center_update_messages
    ) / m

    z_new <- soft_threshold(
      average_beta_plus_u,
      lambda / (m * rho)
    )

    u_list <- lapply(seq_len(m), function(k) {
      u_list[[k]] + beta_list[[k]] - z_new
    })

    r_norm <- sqrt(sum(vapply(
      seq_len(m),
      function(k) sum((beta_list[[k]] - z_new)^2),
      numeric(1)
    )))

    s_norm <- rho * sqrt(m) * sqrt(sum((z_new - z_previous)^2))

    beta_norm <- sqrt(sum(vapply(
      beta_list,
      function(beta_k) sum(beta_k^2),
      numeric(1)
    )))

    z_norm <- sqrt(m) * sqrt(sum(z_new^2))

    u_norm <- sqrt(sum(vapply(
      u_list,
      function(u_k) sum(u_k^2),
      numeric(1)
    )))

    eps_primal <- sqrt(p * m) * eps_abs +
      eps_rel * max(beta_norm, z_norm)

    eps_dual <- sqrt(p * m) * eps_abs +
      eps_rel * rho * u_norm

    if (save_history) {
      history_list[[iteration]] <- data.frame(
        Iteration = iteration,
        Primal_Residual = r_norm,
        Dual_Residual = s_norm,
        Eps_Primal = eps_primal,
        Eps_Dual = eps_dual,
        Primal_Ratio = r_norm / eps_primal,
        Dual_Ratio = s_norm / eps_dual,
        stringsAsFactors = FALSE
      )
    }

    z <- z_new

    if (
      is.finite(r_norm) &&
        is.finite(s_norm) &&
        is.finite(eps_primal) &&
        is.finite(eps_dual) &&
        r_norm <= eps_primal &&
        s_norm <= eps_dual
    ) {
      converged <- TRUE
      break
    }
  }

  history_df <- NULL
  if (save_history) {
    history_df <- bind_rows_safe(history_list[seq_len(iteration)])
  }

  list(
    beta = z,
    u_list = u_list,
    iterations = iteration,
    converged = converged,
    final_primal_residual = r_norm,
    final_dual_residual = s_norm,
    final_eps_primal = eps_primal,
    final_eps_dual = eps_dual,
    final_primal_ratio = r_norm / eps_primal,
    final_dual_ratio = s_norm / eps_dual,
    history = history_df
  )
}

tune_admm_lambda_path <- function(clients_train,
                                  validation_clients,
                                  rho,
                                  lambda_grid,
                                  lambda_max,
                                  p,
                                  m,
                                  eps_abs,
                                  eps_rel,
                                  max_iter) {
  prepared_systems <- prepare_admm_systems(
    clients = clients_train,
    rho = rho,
    p = p
  )

  lambda_order <- sort(unique(lambda_grid), decreasing = TRUE)

  z_warm <- rep(0, p)
  u_warm <- lapply(seq_len(m), function(k) rep(0, p))

  tuning_rows <- vector("list", length(lambda_order))
  path_start <- proc.time()[[3]]

  for (j in seq_along(lambda_order)) {
    lambda_value <- lambda_order[j]
    fit_start <- proc.time()[[3]]

    fit <- run_admm_lasso(
      prepared_systems = prepared_systems,
      lambda = lambda_value,
      rho = rho,
      p = p,
      m = m,
      eps_abs = eps_abs,
      eps_rel = eps_rel,
      max_iter = max_iter,
      z_init = z_warm,
      u_init = u_warm,
      save_history = FALSE
    )

    fit_runtime <- proc.time()[[3]] - fit_start

    if (isTRUE(fit$converged)) {
      z_warm <- fit$beta
      u_warm <- fit$u_list

      # Every center evaluates the consensus model on its own validation
      # records. Only aggregate loss summaries are combined.
      validation_loss <- federated_loss_stats(
        clients = validation_clients,
        beta_estimate = fit$beta
      )
    } else {
      z_warm <- rep(0, p)
      u_warm <- lapply(seq_len(m), function(k) rep(0, p))

      validation_loss <- list(
        MSE = NA_real_,
        SE = NA_real_,
        MAE = NA_real_
      )
    }

    tuning_rows[[j]] <- data.frame(
      Lambda = lambda_value,
      Lambda_Fraction = lambda_value / lambda_max,
      Validation_MSE = validation_loss$MSE,
      Validation_SE = validation_loss$SE,
      Validation_MAE = validation_loss$MAE,
      Converged = fit$converged,
      Iterations = fit$iterations,
      Runtime_Seconds = fit_runtime,
      Final_Primal_Residual = fit$final_primal_residual,
      Final_Dual_Residual = fit$final_dual_residual,
      Final_Eps_Primal = fit$final_eps_primal,
      Final_Eps_Dual = fit$final_eps_dual,
      Final_Primal_Ratio = fit$final_primal_ratio,
      Final_Dual_Ratio = fit$final_dual_ratio,
      stringsAsFactors = FALSE
    )
  }

  tuning_df <- bind_rows_safe(tuning_rows)
  selected <- choose_lambda_indices(tuning_df)

  list(
    prepared_systems = prepared_systems,
    tuning_df = tuning_df,
    selected = selected,
    path_runtime_seconds = proc.time()[[3]] - path_start
  )
}

# Centralized and local glmnet functions
tune_glmnet_lambda_path <- function(X_train,
                                    y_train,
                                    X_validation,
                                    y_validation,
                                    lambda_grid,
                                    lambda_max) {
  lambda_order <- sort(unique(lambda_grid), decreasing = TRUE)
  N <- nrow(X_train)

  path_start <- proc.time()[[3]]

  fit <- glmnet(
    x = X_train,
    y = y_train,
    family = "gaussian",
    alpha = 1,
    lambda = lambda_order / N,
    intercept = FALSE,
    standardize = FALSE,
    thresh = 1e-12,
    maxit = 100000,
    dfmax = ncol(X_train) + 1L,
    pmax = ncol(X_train) + 1L
  )

  total_runtime <- proc.time()[[3]] - path_start

  beta_matrix <- as.matrix(coef(fit))[-1, , drop = FALSE]
  returned_lambda_common <- as.numeric(fit$lambda) * N

  if (ncol(beta_matrix) != length(returned_lambda_common)) {
    stop("Unexpected glmnet coefficient-path dimensions.")
  }

  if (length(returned_lambda_common) != length(lambda_order)) {
    stop(
      "glmnet did not return the complete user-supplied lambda path. ",
      "Check glmnet.control settings."
    )
  }

  maximum_lambda_difference <- max(abs(returned_lambda_common - lambda_order))
  if (maximum_lambda_difference > 1e-8 * max(1, max(lambda_order))) {
    stop("glmnet returned lambda values that do not match the requested grid.")
  }

  validation_loss <- loss_path_stats(
    X = X_validation,
    y = y_validation,
    beta_matrix = beta_matrix
  )

  tuning_df <- data.frame(
    Lambda = returned_lambda_common,
    Lambda_Fraction = returned_lambda_common / lambda_max,
    Validation_MSE = validation_loss$MSE,
    Validation_SE = validation_loss$SE,
    Validation_MAE = validation_loss$MAE,
    Converged = TRUE,
    Iterations = NA_real_,
    Runtime_Seconds = total_runtime,
    Final_Primal_Residual = NA_real_,
    Final_Dual_Residual = NA_real_,
    Final_Eps_Primal = NA_real_,
    Final_Eps_Dual = NA_real_,
    Final_Primal_Ratio = NA_real_,
    Final_Dual_Ratio = NA_real_,
    stringsAsFactors = FALSE
  )

  selected <- choose_lambda_indices(tuning_df)

  list(
    tuning_df = tuning_df,
    beta_matrix = beta_matrix,
    selected = selected,
    path_runtime_seconds = total_runtime
  )
}

get_glmnet_beta_at_index <- function(glmnet_tuning, index) {
  as.vector(glmnet_tuning$beta_matrix[, index])
}

get_glmnet_beta_at_lambda <- function(glmnet_tuning, target_lambda) {
  index <- find_lambda_index(
    glmnet_tuning$tuning_df$Lambda,
    target_lambda
  )
  as.vector(glmnet_tuning$beta_matrix[, index])
}

# Result helper functions
make_success_result_row <- function(replication,
                                    method,
                                    selection_rule,
                                    rho,
                                    beta_estimate,
                                    selected_lambda,
                                    lambda_max,
                                    validation_mse,
                                    validation_se,
                                    true_beta,
                                    selection_threshold,
                                    converged,
                                    test_data = NULL,
                                    test_clients = NULL,
                                    iterations = NA_real_,
                                    tuning_runtime_seconds = NA_real_,
                                    final_runtime_seconds = NA_real_,
                                    final_primal_residual = NA_real_,
                                    final_dual_residual = NA_real_,
                                    final_eps_primal = NA_real_,
                                    final_eps_dual = NA_real_,
                                    final_primal_ratio = NA_real_,
                                    final_dual_ratio = NA_real_,
                                    pooled_beta_difference_same_lambda = NA_real_) {
  selection <- calculate_selection_metrics(
    beta_estimate = beta_estimate,
    true_beta = true_beta,
    threshold = selection_threshold
  )

  if (!is.null(test_clients)) {
    # Data-local test evaluation: centers return only aggregate loss summaries.
    test_loss <- federated_loss_stats(
      clients = test_clients,
      beta_estimate = beta_estimate
    )
  } else if (!is.null(test_data)) {
    # Centralized benchmark evaluation.
    test_loss <- loss_stats_matrix(
      X = test_data$X,
      y = test_data$y,
      beta_estimate = beta_estimate
    )
  } else {
    stop("Either test_data or test_clients must be supplied.")
  }

  data.frame(
    Setting = setting_label,
    p = p,
    n = n_per_center,
    Centers = m,
    Replication = replication,
    Method = method,
    Selection_Rule = selection_rule,
    Rho = rho,
    Lambda = selected_lambda,
    Lambda_Max = lambda_max,
    Lambda_Fraction = selected_lambda / lambda_max,
    Validation_MSE = validation_mse,
    Validation_SE = validation_se,
    TP = selection$TP,
    FP = selection$FP,
    FN = selection$FN,
    TN = selection$TN,
    Sensitivity = selection$Sensitivity,
    Specificity = selection$Specificity,
    Number_Selected = selection$Number_Selected,
    MSE_beta = mean((beta_estimate - true_beta)^2),
    MAE_beta = mean(abs(beta_estimate - true_beta)),
    Test_MSE = test_loss$MSE,
    Test_MAE = test_loss$MAE,
    Converged = as.logical(converged),
    Iterations = iterations,
    Tuning_Runtime_Seconds = tuning_runtime_seconds,
    Final_Runtime_Seconds = final_runtime_seconds,
    Total_Runtime_Seconds = tuning_runtime_seconds + final_runtime_seconds,
    Final_Primal_Residual = final_primal_residual,
    Final_Dual_Residual = final_dual_residual,
    Final_Eps_Primal = final_eps_primal,
    Final_Eps_Dual = final_eps_dual,
    Final_Primal_Ratio = final_primal_ratio,
    Final_Dual_Ratio = final_dual_ratio,
    Max_Abs_Beta_Difference_From_Pooled_At_Same_Lambda =
      pooled_beta_difference_same_lambda,
    Failure_Reason = "",
    stringsAsFactors = FALSE
  )
}

make_failure_result_row <- function(replication,
                                    selection_rule,
                                    rho,
                                    failure_reason,
                                    selected_lambda = NA_real_,
                                    lambda_max = NA_real_,
                                    validation_mse = NA_real_,
                                    validation_se = NA_real_,
                                    iterations = NA_real_,
                                    tuning_runtime_seconds = NA_real_,
                                    final_runtime_seconds = NA_real_,
                                    final_primal_residual = NA_real_,
                                    final_dual_residual = NA_real_,
                                    final_eps_primal = NA_real_,
                                    final_eps_dual = NA_real_,
                                    final_primal_ratio = NA_real_,
                                    final_dual_ratio = NA_real_) {
  data.frame(
    Setting = setting_label,
    p = p,
    n = n_per_center,
    Centers = m,
    Replication = replication,
    Method = "ADMM-LASSO",
    Selection_Rule = selection_rule,
    Rho = rho,
    Lambda = selected_lambda,
    Lambda_Max = lambda_max,
    Lambda_Fraction = if (
      is.finite(selected_lambda) && is.finite(lambda_max) && lambda_max > 0
    ) selected_lambda / lambda_max else NA_real_,
    Validation_MSE = validation_mse,
    Validation_SE = validation_se,
    TP = NA_real_,
    FP = NA_real_,
    FN = NA_real_,
    TN = NA_real_,
    Sensitivity = NA_real_,
    Specificity = NA_real_,
    Number_Selected = NA_real_,
    MSE_beta = NA_real_,
    MAE_beta = NA_real_,
    Test_MSE = NA_real_,
    Test_MAE = NA_real_,
    Converged = FALSE,
    Iterations = iterations,
    Tuning_Runtime_Seconds = tuning_runtime_seconds,
    Final_Runtime_Seconds = final_runtime_seconds,
    Total_Runtime_Seconds = tuning_runtime_seconds + final_runtime_seconds,
    Final_Primal_Residual = final_primal_residual,
    Final_Dual_Residual = final_dual_residual,
    Final_Eps_Primal = final_eps_primal,
    Final_Eps_Dual = final_eps_dual,
    Final_Primal_Ratio = final_primal_ratio,
    Final_Dual_Ratio = final_dual_ratio,
    Max_Abs_Beta_Difference_From_Pooled_At_Same_Lambda = NA_real_,
    Failure_Reason = failure_reason,
    stringsAsFactors = FALSE
  )
}

make_coefficient_rows <- function(replication,
                                  method,
                                  selection_rule,
                                  rho,
                                  beta_estimate,
                                  true_beta) {
  data.frame(
    Setting = setting_label,
    p = p,
    n = n_per_center,
    Replication = replication,
    Method = method,
    Selection_Rule = selection_rule,
    Rho = rho,
    Predictor_Index = seq_along(beta_estimate),
    Beta_Hat = as.vector(beta_estimate),
    Beta_True = as.vector(true_beta),
    Selected = abs(beta_estimate) > selection_threshold,
    Truly_Active = abs(true_beta) > selection_threshold,
    stringsAsFactors = FALSE
  )
}

# Windows version runs the 100 simulation replications.
main_results <- list()
lambda_tuning_results <- list()
coefficient_results <- list()
convergence_history <- list()

main_pointer <- 1L
tuning_pointer <- 1L
coefficient_pointer <- 1L
history_pointer <- 1L

chunk_start <- proc.time()[[3]]

for (replication in seq_len(total_rep)) {
  replication_start <- proc.time()[[3]]
  set.seed(base_seed + replication)

  cat("\n------------------------------------------------------------\n")
  cat("Replication", replication, "of", total_rep, "\n")
  cat("------------------------------------------------------------\n")

  # Identical data are used for all methods and all rho values.
  clients_train <- lapply(seq_len(m), function(center_index) {
    simulate_data(n_per_center, p, true_beta, noise_sd)
  })
  clients_validation <- lapply(seq_len(m), function(center_index) {
    simulate_data(n_per_center, p, true_beta, noise_sd)
  })
  clients_test <- lapply(seq_len(m), function(center_index) {
    simulate_data(n_per_center, p, true_beta, noise_sd)
  })

  # Stacked records are created only for the intentionally centralized
  # LASSO-Pool benchmark. The proposed ADMM method uses the client lists.
  pooled_train <- stack_clients(clients_train)
  pooled_validation <- stack_clients(clients_validation)
  pooled_test <- stack_clients(clients_test)

  # Runtime audit: center-level aggregate losses must exactly reproduce the
  # stacked calculation without transmitting individual-level records.
  audit_beta <- rep(0, p)
  assert_federated_loss_equivalence(
    clients = clients_validation,
    stacked_data = pooled_validation,
    beta_estimate = audit_beta
  )
  assert_federated_loss_equivalence(
    clients = clients_test,
    stacked_data = pooled_test,
    beta_estimate = audit_beta
  )

  global_lambda <- make_global_lambda_grid(
    clients_train = clients_train,
    lambda_fraction_grid = lambda_fraction_grid
  )

  #Pooled LASSO: exact same global lambda grid as ADMM-LASSO
  pooled_tuning <- tune_glmnet_lambda_path(
    X_train = pooled_train$X,
    y_train = pooled_train$y,
    X_validation = pooled_validation$X,
    y_validation = pooled_validation$y,
    lambda_grid = global_lambda$lambda_grid,
    lambda_max = global_lambda$lambda_max
  )

  pooled_tuning_output <- pooled_tuning$tuning_df
  pooled_tuning_output$Setting <- setting_label
  pooled_tuning_output$p <- p
  pooled_tuning_output$n <- n_per_center
  pooled_tuning_output$Replication <- replication
  pooled_tuning_output$Method <- "LASSO-Pool"
  pooled_tuning_output$Rho <- NA_real_
  lambda_tuning_results[[tuning_pointer]] <- pooled_tuning_output
  tuning_pointer <- tuning_pointer + 1L

  for (rule_name in selection_rules) {
    selected_index <- pooled_tuning$selected[[rule_name]]
    if (is.na(selected_index)) {
      stop("Pooled LASSO did not produce a valid ", rule_name, " lambda.")
    }

    beta_hat <- get_glmnet_beta_at_index(pooled_tuning, selected_index)
    tuning_row <- pooled_tuning$tuning_df[selected_index, , drop = FALSE]

    main_results[[main_pointer]] <- make_success_result_row(
      replication = replication,
      method = "LASSO-Pool",
      selection_rule = rule_name,
      rho = NA_real_,
      beta_estimate = beta_hat,
      selected_lambda = tuning_row$Lambda,
      lambda_max = global_lambda$lambda_max,
      validation_mse = tuning_row$Validation_MSE,
      validation_se = tuning_row$Validation_SE,
      test_data = pooled_test,
      true_beta = true_beta,
      selection_threshold = selection_threshold,
      converged = TRUE,
      tuning_runtime_seconds = pooled_tuning$path_runtime_seconds,
      final_runtime_seconds = 0
    )
    main_pointer <- main_pointer + 1L

    coefficient_results[[coefficient_pointer]] <- make_coefficient_rows(
      replication, "LASSO-Pool", rule_name, NA_real_, beta_hat, true_beta
    )
    coefficient_pointer <- coefficient_pointer + 1L
  }

  ######################### Single-center LASSO models ################
  
  for (center_index in seq_len(m)) {
    local_lambda <- make_local_lambda_grid(
      X = clients_train[[center_index]]$X,
      y = clients_train[[center_index]]$y,
      lambda_fraction_grid = lambda_fraction_grid
    )

    local_tuning <- tune_glmnet_lambda_path(
      X_train = clients_train[[center_index]]$X,
      y_train = clients_train[[center_index]]$y,
      X_validation = clients_validation[[center_index]]$X,
      y_validation = clients_validation[[center_index]]$y,
      lambda_grid = local_lambda$lambda_grid,
      lambda_max = local_lambda$lambda_max
    )

    method_name <- paste0("LASSO-", center_index)

    local_tuning_output <- local_tuning$tuning_df
    local_tuning_output$Setting <- setting_label
    local_tuning_output$p <- p
    local_tuning_output$n <- n_per_center
    local_tuning_output$Replication <- replication
    local_tuning_output$Method <- method_name
    local_tuning_output$Rho <- NA_real_
    lambda_tuning_results[[tuning_pointer]] <- local_tuning_output
    tuning_pointer <- tuning_pointer + 1L

    for (rule_name in selection_rules) {
      selected_index <- local_tuning$selected[[rule_name]]
      if (is.na(selected_index)) {
        stop(method_name, " did not produce a valid ", rule_name, " lambda.")
      }

      beta_hat <- get_glmnet_beta_at_index(local_tuning, selected_index)
      tuning_row <- local_tuning$tuning_df[selected_index, , drop = FALSE]

      main_results[[main_pointer]] <- make_success_result_row(
        replication = replication,
        method = method_name,
        selection_rule = rule_name,
        rho = NA_real_,
        beta_estimate = beta_hat,
        selected_lambda = tuning_row$Lambda,
        lambda_max = local_lambda$lambda_max,
        validation_mse = tuning_row$Validation_MSE,
        validation_se = tuning_row$Validation_SE,
        test_clients = clients_test,
        true_beta = true_beta,
        selection_threshold = selection_threshold,
        converged = TRUE,
        tuning_runtime_seconds = local_tuning$path_runtime_seconds,
        final_runtime_seconds = 0
      )
      main_pointer <- main_pointer + 1L

      coefficient_results[[coefficient_pointer]] <- make_coefficient_rows(
        replication, method_name, rule_name, NA_real_, beta_hat, true_beta
      )
      coefficient_pointer <- coefficient_pointer + 1L
    }
  }

  # Federated ADMM-LASSO for every rho for checking the rho sensitivity
  
  for (rho in rho_values) {
    cat("  Running ADMM-LASSO for rho =", rho, "\n")

    admm_tuning <- tune_admm_lambda_path(
      clients_train = clients_train,
      validation_clients = clients_validation,
      rho = rho,
      lambda_grid = global_lambda$lambda_grid,
      lambda_max = global_lambda$lambda_max,
      p = p,
      m = m,
      eps_abs = eps_abs,
      eps_rel = eps_rel,
      max_iter = max_admm_iter
    )

    admm_tuning_output <- admm_tuning$tuning_df
    admm_tuning_output$Setting <- setting_label
    admm_tuning_output$p <- p
    admm_tuning_output$n <- n_per_center
    admm_tuning_output$Replication <- replication
    admm_tuning_output$Method <- "ADMM-LASSO"
    admm_tuning_output$Rho <- rho
    lambda_tuning_results[[tuning_pointer]] <- admm_tuning_output
    tuning_pointer <- tuning_pointer + 1L

    for (rule_name in selection_rules) {
      selected_index <- admm_tuning$selected[[rule_name]]

      if (is.na(selected_index)) {
        main_results[[main_pointer]] <- make_failure_result_row(
          replication = replication,
          selection_rule = rule_name,
          rho = rho,
          failure_reason = "No converged lambda was eligible for selection",
          lambda_max = global_lambda$lambda_max,
          tuning_runtime_seconds = admm_tuning$path_runtime_seconds,
          final_runtime_seconds = 0
        )
        main_pointer <- main_pointer + 1L
        next
      }

      selected_tuning_row <- admm_tuning$tuning_df[
        selected_index, , drop = FALSE
      ]
      selected_lambda <- selected_tuning_row$Lambda

      final_start <- proc.time()[[3]]
      final_fit <- run_admm_lasso(
        prepared_systems = admm_tuning$prepared_systems,
        lambda = selected_lambda,
        rho = rho,
        p = p,
        m = m,
        eps_abs = eps_abs,
        eps_rel = eps_rel,
        max_iter = max_admm_iter,
        z_init = NULL,
        u_init = NULL,
        save_history = TRUE
      )
      final_runtime <- proc.time()[[3]] - final_start

      # converged. This permits transparent diagnosis of slow or failed rho
      history_output <- final_fit$history
      history_output$Setting <- setting_label
      history_output$p <- p
      history_output$n <- n_per_center
      history_output$Replication <- replication
      history_output$Method <- "ADMM-LASSO"
      history_output$Selection_Rule <- rule_name
      history_output$Rho <- rho
      history_output$Lambda <- selected_lambda
      history_output$Lambda_Max <- global_lambda$lambda_max
      history_output$Lambda_Fraction <-
        selected_lambda / global_lambda$lambda_max
      history_output$Converged <- isTRUE(final_fit$converged)
      history_output$Final_Iterations <- final_fit$iterations

      convergence_history[[history_pointer]] <- history_output
      history_pointer <- history_pointer + 1L

      if (!isTRUE(final_fit$converged)) {
        main_results[[main_pointer]] <- make_failure_result_row(
          replication = replication,
          selection_rule = rule_name,
          rho = rho,
          failure_reason = paste0(
            "Final selected-lambda fit did not converge by ",
            max_admm_iter, " iterations"
          ),
          selected_lambda = selected_lambda,
          lambda_max = global_lambda$lambda_max,
          validation_mse = selected_tuning_row$Validation_MSE,
          validation_se = selected_tuning_row$Validation_SE,
          iterations = final_fit$iterations,
          tuning_runtime_seconds = admm_tuning$path_runtime_seconds,
          final_runtime_seconds = final_runtime,
          final_primal_residual = final_fit$final_primal_residual,
          final_dual_residual = final_fit$final_dual_residual,
          final_eps_primal = final_fit$final_eps_primal,
          final_eps_dual = final_fit$final_eps_dual,
          final_primal_ratio = final_fit$final_primal_ratio,
          final_dual_ratio = final_fit$final_dual_ratio
        )
        main_pointer <- main_pointer + 1L
        next
      }

      pooled_beta_same_lambda <- get_glmnet_beta_at_lambda(
        pooled_tuning,
        selected_lambda
      )
      beta_difference <- max(abs(final_fit$beta - pooled_beta_same_lambda))

      main_results[[main_pointer]] <- make_success_result_row(
        replication = replication,
        method = "ADMM-LASSO",
        selection_rule = rule_name,
        rho = rho,
        beta_estimate = final_fit$beta,
        selected_lambda = selected_lambda,
        lambda_max = global_lambda$lambda_max,
        validation_mse = selected_tuning_row$Validation_MSE,
        validation_se = selected_tuning_row$Validation_SE,
        test_clients = clients_test,
        true_beta = true_beta,
        selection_threshold = selection_threshold,
        converged = TRUE,
        iterations = final_fit$iterations,
        tuning_runtime_seconds = admm_tuning$path_runtime_seconds,
        final_runtime_seconds = final_runtime,
        final_primal_residual = final_fit$final_primal_residual,
        final_dual_residual = final_fit$final_dual_residual,
        final_eps_primal = final_fit$final_eps_primal,
        final_eps_dual = final_fit$final_eps_dual,
        final_primal_ratio = final_fit$final_primal_ratio,
        final_dual_ratio = final_fit$final_dual_ratio,
        pooled_beta_difference_same_lambda = beta_difference
      )
      main_pointer <- main_pointer + 1L

      coefficient_results[[coefficient_pointer]] <- make_coefficient_rows(
        replication, "ADMM-LASSO", rule_name, rho, final_fit$beta, true_beta
      )
      coefficient_pointer <- coefficient_pointer + 1L

    }
  }

  cat(
    "Replication", replication, "completed in",
    round(proc.time()[[3]] - replication_start, 2), "seconds.\n"
  )
}

#combined desktop outputs
main_results_df <- bind_rows_safe(main_results)
lambda_tuning_df <- bind_rows_safe(lambda_tuning_results)
coefficient_results_df <- bind_rows_safe(coefficient_results)
convergence_history_df <- bind_rows_safe(convergence_history)

### Primary MinMSE summary at rho = 1
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x)
}

method_order <- c("LASSO-1", "LASSO-2", "LASSO-Pool", "ADMM-LASSO")

primary_data <- rbind(
  subset(
    main_results_df,
    Method != "ADMM-LASSO" & Selection_Rule == primary_rule
  ),
  subset(
    main_results_df,
    Method == "ADMM-LASSO" &
      Selection_Rule == primary_rule &
      Rho == primary_rho
  )
)

summarise_one_method <- function(method_name) {
  x <- primary_data[primary_data$Method == method_name, , drop = FALSE]

  data.frame(
    Method = method_name,
    N_Rep = length(unique(x$Replication)),
    Mean_FP = safe_mean(x$FP),
    SD_FP = safe_sd(x$FP),
    Mean_FN = safe_mean(x$FN),
    SD_FN = safe_sd(x$FN),
    Mean_Sensitivity = safe_mean(x$Sensitivity),
    SD_Sensitivity = safe_sd(x$Sensitivity),
    Mean_Specificity = safe_mean(x$Specificity),
    SD_Specificity = safe_sd(x$Specificity),
    Mean_MSE_beta = safe_mean(x$MSE_beta),
    SD_MSE_beta = safe_sd(x$MSE_beta),
    Mean_MAE_beta = safe_mean(x$MAE_beta),
    SD_MAE_beta = safe_sd(x$MAE_beta),
    Mean_Test_MSE = safe_mean(x$Test_MSE),
    SD_Test_MSE = safe_sd(x$Test_MSE),
    Mean_Test_MAE = safe_mean(x$Test_MAE),
    SD_Test_MAE = safe_sd(x$Test_MAE),
    stringsAsFactors = FALSE
  )
}

primary_summary <- do.call(
  rbind,
  lapply(method_order, summarise_one_method)
)



# A publication-friendly mean (SD) table.
format_mean_sd <- function(mean_value, sd_value, digits = 4L) {
  if (!is.finite(mean_value)) return(NA_character_)
  if (!is.finite(sd_value)) {
    return(formatC(mean_value, format = "f", digits = digits))
  }
  paste0(
    formatC(mean_value, format = "f", digits = digits),
    " (",
    formatC(sd_value, format = "f", digits = digits),
    ")"
  )
}

primary_formatted <- data.frame(
  Method = primary_summary$Method,
  FP = mapply(format_mean_sd,
              primary_summary$Mean_FP,
              primary_summary$SD_FP,
              MoreArgs = list(digits = 2L)),
  FN = mapply(format_mean_sd,
              primary_summary$Mean_FN,
              primary_summary$SD_FN,
              MoreArgs = list(digits = 2L)),
  Sensitivity = mapply(format_mean_sd,
                       primary_summary$Mean_Sensitivity,
                       primary_summary$SD_Sensitivity,
                       MoreArgs = list(digits = 4L)),
  Specificity = mapply(format_mean_sd,
                       primary_summary$Mean_Specificity,
                       primary_summary$SD_Specificity,
                       MoreArgs = list(digits = 4L)),
  MSE_beta = mapply(format_mean_sd,
                    primary_summary$Mean_MSE_beta,
                    primary_summary$SD_MSE_beta,
                    MoreArgs = list(digits = 6L)),
  MAE_beta = mapply(format_mean_sd,
                    primary_summary$Mean_MAE_beta,
                    primary_summary$SD_MAE_beta,
                    MoreArgs = list(digits = 6L)),
  Test_MSE = mapply(format_mean_sd,
                    primary_summary$Mean_Test_MSE,
                    primary_summary$SD_Test_MSE,
                    MoreArgs = list(digits = 6L)),
  stringsAsFactors = FALSE
)

# Paired comparisons at  rho=1
safe_paired_p <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 2L) return(NA_real_)

  d <- x - y
  tol <- sqrt(.Machine$double.eps)

  if (all(abs(d) < tol)) return(NA_real_)
  if (sd(d) < tol) return(0)

  t.test(x, y, paired = TRUE)$p.value
}

admm_primary <- primary_data[
  primary_data$Method == "ADMM-LASSO",
  ,
  drop = FALSE
]

comparators <- c("LASSO-1", "LASSO-2", "LASSO-Pool")
test_metrics <- c("Test_MSE", "MSE_beta", "MAE_beta", "Specificity")
paired_rows <- list()
pointer <- 1L

for (comparison_method in comparators) {
  comparator_data <- primary_data[
    primary_data$Method == comparison_method,
    ,
    drop = FALSE
  ]

  merged <- merge(
    admm_primary,
    comparator_data,
    by = "Replication",
    suffixes = c("_ADMM", "_Comparator")
  )

  for (metric_name in test_metrics) {
    x <- merged[[paste0(metric_name, "_ADMM")]]
    y <- merged[[paste0(metric_name, "_Comparator")]]

    paired_rows[[pointer]] <- data.frame(
      Comparison = paste("ADMM-LASSO vs", comparison_method),
      Metric = metric_name,
      Mean_ADMM = safe_mean(x),
      Mean_Comparator = safe_mean(y),
      P_Value = safe_paired_p(x, y),
      stringsAsFactors = FALSE
    )
    pointer <- pointer + 1L
  }
}

paired_tests <- do.call(rbind, paired_rows)



# Rho-sensitivity summary 
admm_minmse <- subset(
  main_results_df,
  Method == "ADMM-LASSO" & Selection_Rule == primary_rule
)

rho_groups <- split(admm_minmse, admm_minmse$Rho)

rho_summary <- do.call(
  rbind,
  lapply(names(rho_groups), function(rho_name) {
    x <- rho_groups[[rho_name]]
    data.frame(
      Rho = as.numeric(rho_name),
      N_Rep = length(unique(x$Replication)),
      Convergence_Rate = mean(x$Converged %in% TRUE),
      Mean_Iterations = safe_mean(x$Iterations),
      Median_Iterations = median(x$Iterations[is.finite(x$Iterations)]),
      Mean_Test_MSE = safe_mean(x$Test_MSE),
      SD_Test_MSE = safe_sd(x$Test_MSE),
      stringsAsFactors = FALSE
    )
  })
)

rho_summary <- rho_summary[order(rho_summary$Rho), , drop = FALSE]




## rho convergence summary
if (nrow(convergence_history_df) > 0L) {
  primary_history <- convergence_history_df[
    convergence_history_df$Selection_Rule == primary_rule &
      convergence_history_df$Rho == primary_rho &
      convergence_history_df$Converged %in% TRUE,
    ,
    drop = FALSE
  ]

  if (nrow(primary_history) > 0L) {
    primary_history$Primal_Ratio <-
      primary_history$Primal_Residual / primary_history$Eps_Primal
    primary_history$Dual_Ratio <-
      primary_history$Dual_Residual / primary_history$Eps_Dual

    convergence_iterations <- sort(unique(primary_history$Iteration))

    convergence_summary <- do.call(
      rbind,
      lapply(convergence_iterations, function(iteration_value) {
        x <- primary_history[
          primary_history$Iteration == iteration_value,
          ,
          drop = FALSE
        ]

        data.frame(
          Iteration = iteration_value,
          N_Trajectories = length(unique(x$Replication)),
          Mean_Primal_Ratio = safe_mean(x$Primal_Ratio),
          Mean_Dual_Ratio = safe_mean(x$Dual_Ratio),
          stringsAsFactors = FALSE
        )
      })
    )
  }
}

## settings 
settings_df <- data.frame(
  Parameter = c(
    "p", "n_per_center", "m", "total_rep",
    "primary_rho", "rho_values", "selection_rule",
    "active_beta_values", "noise_sd",
    "n_lambda", "lambda_min_ratio",
    "selection_threshold", "eps_abs", "eps_rel",
    "max_admm_iter", "base_seed"
  ),
  Value = c(
    p, n_per_center, m, total_rep,
    primary_rho, paste(rho_values, collapse = ";"), primary_rule,
    paste(active_beta_values, collapse = ";"), noise_sd,
    n_lambda, lambda_min_ratio,
    selection_threshold, eps_abs, eps_rel,
    max_admm_iter, base_seed
  ),
  stringsAsFactors = FALSE
)



cat("\n============================================================\n")
cat("RUN COMPLETED\n")
cat("============================================================\n")
cat("Primary setting:", setting_label, "\n")
cat("Primary rho:", primary_rho, "\n")
cat("Results folder:", normalizePath(output_dir, mustWork = FALSE), "\n")
cat("\nPrimary mean/SD results:\n")
print(primary_summary, row.names = FALSE)
cat("\nFiles written successfully.\n")

#### View Table 1
View(primary_formatted)
### View Table 2
View(paired_tests)

View(rho_summary)


# Figure 1:CONVERGENCE PLOT 
plot(
  convergence_summary_manuscript$Iteration,
  convergence_summary_manuscript$Mean_Dual_Ratio,
  type = "l",
  log = "y",
  col = "#F8766D",
  lwd = 2,
  lty = 1,
  xlab = "ADMM iteration",
  ylab = "Mean residual-to-tolerance ratio",
  main = "p = 50, n = 100",
  ylim = range(
    convergence_summary_manuscript$Mean_Primal_Ratio,
    convergence_summary_manuscript$Mean_Dual_Ratio,
    1
  )
)

# Primal residual
lines(
  convergence_summary_manuscript$Iteration,
  convergence_summary_manuscript$Mean_Primal_Ratio,
  col = "#00BFC4",
  lwd = 2,
  lty = 2
)

# Convergence threshold
abline(
  h = 1,
  col = "gray40",
  lwd = 1.5,
  lty = 3
)

legend(
  "topright",
  legend = c(
    "Dual residual",
    "Primal residual",
    "Convergence threshold"
  ),
  col = c(
    "#F8766D",
    "#00BFC4",
    "gray40"
  ),
  lty = c(1, 2, 3),
  lwd = c(2, 2, 1.5),
  bty = "n"
)
####################################################################################