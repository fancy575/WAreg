#' Compute the Estimating Equation (Internal)
#'
#' @description
#' This internal function computes the estimating equation \eqn{U} based on the expanded covariate
#' matrix from \code{fd_dt} and the coefficient vector \code{beta}. It extracts the relevant
#' spline-expanded covariate columns (using \code{cov_vars}) from \code{fd_dt}, computes the linear
#' predictor, and then calculates \eqn{U} as the column means of the row-wise product of the covariate
#' matrix and a scalar term.
#'
#' @param beta A numeric vector of coefficients.
#' @param fd_dt A data frame containing the required columns (e.g., \code{V_tau}, \code{G_X_tau},
#'   \code{L}, \code{X_tau}) along with the expanded covariate columns.
#' @param cov_vars A character vector of covariate names used to construct the spline-expanded
#'   covariate columns.
#' @param link link function for the while-alive hazard rate, default is "log"
#'
#' @return A numeric vector representing the estimating equation \eqn{U}.
#'
#' @keywords internal
est_eq <- function(beta, fd_dt, cov_vars,link) {
  # Check that required columns exist
  req_cols <- c("V_tau", "G_X_tau", "L", "X_tau")
  missing_cols <- setdiff(req_cols, names(fd_dt))
  if(length(missing_cols) > 0)
    stop("The following required columns are missing in fd_dt: ", paste(missing_cols, collapse = ", "))

  pattern <- paste0("^(sp_int_|(", paste(cov_vars, collapse = "|"), ")_sp_)")
  covariate_columns <- grep(pattern, names(fd_dt), value = TRUE)

  Z_i <- as.matrix(fd_dt[, covariate_columns])

  if(link == "log"){
    exp_Z_beta <- exp(Z_i %*% beta)

  }else if(link == "identity"){
    exp_Z_beta <- Z_i %*% beta

  }

  difference <- fd_dt$L - exp_Z_beta * fd_dt$X_tau

  scalar_part <- (fd_dt$V_tau / fd_dt$G_X_tau) * difference

  U <- colMeans(sweep(Z_i, 1, scalar_part, `*`), na.rm = TRUE)

  return(U)
}

#' Calculate Variance V (Internal)
#'
#' @description
#' This function calculates the variance estimate \eqn{V} for the WA regression estimator.
#' It uses the expanded covariate matrix (constructed based on spline basis functions)
#' from \code{fd_dt}, the coefficient vector \code{beta}, and the covariate names in \code{cov_vars}.
#'
#' @param fd_dt A data frame containing the following required columns:
#'   \itemize{
#'     \item \code{etime} - event time,
#'     \item \code{delta} - censoring indicator (0 indicates censoring),
#'     \item \code{tau} - the evaluated spline time,
#'     \item \code{X_tau} - the minimum event time for each subject,
#'     \item \code{V_tau} - numerator for censoring weights,
#'     \item \code{G_X_tau} - denominator for censoring weights,
#'     \item \code{j} - subject id,
#'     \item \code{i} - cluster id,
#'     \item and the expanded spline covariate columns (with names matching a pattern based on \code{cov_vars}).
#'   }
#' @param beta A numeric vector of estimated coefficients corresponding to the expanded covariate columns.
#' @param cov_vars A character vector of covariate names used to construct the spline-expanded covariate columns.
#'
#' @return A matrix representing the estimated variance \eqn{V}.
#'
#' @keywords internal
#'
cal_V <- function(fd_dt, beta, cov_vars,link) {

  #    (Assumes intercept terms are named with "sp_int_" and covariate expansions with "Z1_sp_" or "Z2_sp_")
  pattern <- paste0("^(sp_int_|(", paste(cov_vars, collapse = "|"), ")_sp_)")
  covariate_columns <- grep(pattern, names(fd_dt), value = TRUE)

  # Form the expanded covariate matrix.
  Z_i <- as.matrix(fd_dt[, covariate_columns])
  n <- nrow(fd_dt)
  p <- ncol(Z_i)

  # 2. Compute linear predictor and related quantities (vectorized).
  Z_beta <- as.vector(Z_i %*% beta)
  if(link == "log"){
    g_inv_Z_beta <- exp(Z_beta)
    g_inv_prime_Z_beta <- g_inv_Z_beta
  }else if(link == "identity"){
    g_inv_Z_beta <- Z_beta
    g_inv_prime_Z_beta <- 1
  }


  fd_dt <- fd_dt %>%
    mutate(
      Z_beta = Z_beta,
      g_inv_Z_beta = g_inv_Z_beta,
      g_inv_prime_Z_beta = g_inv_prime_Z_beta,
      factor = V_tau / G_X_tau,
      adjusted_L = L - g_inv_Z_beta * X_tau
    )

  # 3. Compute A_beta
  W <- fd_dt$factor * fd_dt$g_inv_prime_Z_beta * fd_dt$X_tau
  A_beta <- crossprod(Z_i * W, Z_i) / n

  # 4. Compute Nelson–Aalen quantities
  unique_times <- sort(unique(fd_dt$etime))
  m_time <- length(unique_times)
  Y_counts <- sapply(unique_times, function(t) sum(fd_dt$etime >= t))
  dNC_counts <- sapply(unique_times, function(t) sum(fd_dt$etime == t & fd_dt$delta == 0))
  dLambda_C <- dNC_counts / Y_counts
  avg_at_risk <- Y_counts / n

  # 5. Compute psi for each observation.
  psi <- matrix(0, n, m_time)
  for (i in 1:n) {
    I_t <- as.numeric(fd_dt$etime[i] >= unique_times)
    I_eq <- as.numeric(fd_dt$etime[i] == unique_times & fd_dt$delta[i] == 0)
    dM_i <- I_eq - I_t * dLambda_C
    # Compute psi_i as cumulative sum of (dM / avg_at_risk)
    psi[i,] <- cumsum(dM_i / avg_at_risk)
  }

  # Calculate K for each covariate.
  dK_elems <- list()
  for (i in 1:p) {
    dK_elems[[i]] <- fd_dt$factor * Z_i[, i] * fd_dt$adjusted_L
  }

  KI <- outer(fd_dt$X_tau, unique_times, FUN = "==") * 1

  dK_expands <- lapply(dK_elems, function(dK_elem) {
    dK <- apply(KI, 2, function(k) mean(k * dK_elem))
    matrix(dK, nrow = nrow(fd_dt), ncol = length(dK), byrow = TRUE)
  })

  Q_integrals <- lapply(dK_expands, function(dK_expand) {
    psi * dK_expand
  })

  # For each subject (row), integrate Q_integrals up to their tau.
  Qs <- sapply(1:p, function(i) {
    sapply(1:nrow(fd_dt), function(row_idx) {
      valid_times <- unique_times[unique_times <= fd_dt$tau[row_idx]]
      Q_values <- Q_integrals[[i]][row_idx, unique_times <= fd_dt$tau[row_idx]]
      if (length(valid_times) > 1) {
        integrate_value <- sum(diff(valid_times) * (Q_values[-length(Q_values)] + Q_values[-1]) / 2)
      } else {
        integrate_value <- sum(Q_values)
      }
      dK_elems[[i]][row_idx] - integrate_value
    })
  })

  # Optionally, you could summarize Qs by subject.
  # Here, we create a data frame with Qs and attach subject and cluster identifiers.
  Qs_df <- as.data.frame(Qs) %>%
    mutate(j = fd_dt$j, i = fd_dt$i)

  # Group by cluster (here, 'i') and compute the mean of each Q column.
  summed_by_subject <- Qs_df %>%
    group_by(i) %>%
    summarise(across(starts_with("V"), mean, na.rm = TRUE)) %>%
    select(-i)

  B_beta <- as.matrix(t(summed_by_subject)) %*% as.matrix(summed_by_subject) / length(unique(fd_dt$i))
  V <- solve(A_beta) %*% B_beta %*% solve(A_beta) / length(unique(fd_dt$i))


  return(V)
}

#' Calculate Prediction Error for WA Regression (Internal)
#'
#' @description
#' This internal function calculates a prediction error (PE) measure for the WA regression model.
#' For each subject in the data, it computes the predicted effect at time points specified in
#' \code{t_span} using the spline basis (via \code{gen_cov_prod}) and compares it to the observed
#' outcome (assumed to be in the \code{N_R} column). The overall prediction error is computed as the
#' mean squared error across subjects and time points.
#'
#' @param model A survival model formula (e.g., \code{Surv(time, type) ~ Z1 + Z2}).
#' @param data A data frame containing the data used in the analysis, including the observed outcome
#'   in the \code{N_R} column.
#' @param cluster Optional; a string giving the cluster column name (if applicable).
#' @param id A string giving the subject ID column name.
#' @param wr A numeric vector of weights for recurrent events.
#' @param wd A numeric weight for the terminal event.
#' @param link Link function for while-alive loss rate.
#' @param cens_mod A survival model formula for the censoring distribution.
#' @param sp_knots A numeric vector of spline knot locations.
#' @param degree A non-negative integer specifying the polynomial degree for the spline.
#' @param beta A numeric vector of estimated coefficients corresponding to the spline-expanded basis.
#' @param t_span A numeric vector of time points at which to evaluate prediction error.
#'
#' @return A list with components:
#' \describe{
#'   \item{PE}{The overall prediction error (mean squared error).}
#'   \item{pred_err}{A matrix of squared errors with subjects as rows and times as columns.}
#' }
#'
#' @keywords internal

PE_cal <- function(model, data,cluster=NULL,id,wr,wd,link,
                   cens_mod,
                   sp_knots, degree,beta,t_span) {

  lhs_vars <- all.vars(formula(model)[[2]])
  cov_vars <- all.vars(formula(model)[[3]])

  ## Check censoring model

  lhs_cens <- all.vars(formula(cens_mod)[[2]])
  cov_cens <- all.vars(formula(cens_mod)[[3]])


  fixed_max <- max(data[[ lhs_vars[2] ]], na.rm = TRUE)

  df1 <- data %>%
    filter(.data[[ lhs_vars[2] ]] %in% c(fixed_max, 0)) %>%
    mutate(priority = if_else(.data[[ lhs_vars[2] ]] == fixed_max, 1L, 2L)) %>%
    group_by(across(all_of(id))) %>%
    arrange(priority, .data[[ lhs_vars[1] ]]) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(-priority) %>%
    mutate(delta=if_else(.data[[ lhs_vars[2] ]] == fixed_max, 1,0 ))

  df2 <- data %>%
    group_by(across(all_of(id))) %>%
    slice(1) %>%
    ungroup() %>%
    left_join(
      data %>%
        filter(!(!!sym(lhs_vars[2]) %in% c(0, fixed_max))) %>%
        select(all_of(id),
               time_new = !!sym(lhs_vars[1]),
               type_new = !!sym(lhs_vars[2])),
      by = id
    ) %>%
    mutate(
      !!lhs_vars[1] := time_new,
      !!lhs_vars[2] := type_new
    ) %>%
    select(-time_new, -type_new)

  ## manipulate data based  on tau

  diff_pre <- 0
  PE_sum <- 0

  for(t in 1:length(t_span)){
    tau = t_span[t]


    df_long <- df1 %>%
      mutate(
        tau=tau,
        X_tau = pmin(!!sym(lhs_vars[1]), tau),
        V_tau = if_else(!!sym(lhs_vars[1]) <= tau, delta, 1)
      )

    cens_mod_n <- reformulate(
      attr(terms(cens_mod), "term.labels"),
      response = call("Surv",
                      formula(cens_mod)[[2]][[2]],
                      call("==", as.name("delta"), 0))
    )

    cens_mod_n <- Surv(time, delta==0) ~ 1
    cens_mod_n <- Surv(time, delta==0) ~ Z2


    if (length(attr(terms(cens_mod_n), "term.labels")) > 0) {
      fit_model <- coxph(cens_mod_n, data = df1)
    } else {
      fit_model <- survfit(cens_mod_n, data = df1)
    }


    if (inherits(fit_model, "coxph")) {
      # Compute the baseline survival curve from the Cox model fit
      base_fit <- survfit(fit_model)
      # Build a step function for the baseline survival probability
      S0 <- approxfun(base_fit$time, base_fit$surv, method = "constant", f = 0, rule = 2)
      # Predict the linear predictors for all rows in df2_longer at once
      lp <- predict(fit_model, newdata = df_long, type = "lp")
      # Compute individual survival probabilities at X_tau
      df_long$G_X_tau <- S0(df_long$X_tau)^(exp(lp))
    } else {
      # For a Kaplan-Meier model (no covariates), just build the baseline step function
      S0 <- approxfun(fit_model$time, fit_model$surv, method = "constant", f = 0, rule = 2)
      df_long$G_X_tau <- S0(df_long$X_tau)
    }

    df_long <- df_long  %>%
      mutate(N_D = ifelse(time <= tau,1,0)*wd )

    df_long$N_R <- mapply(function(subj, tau_val) {
      events <- df2[df2[[ id ]] == subj & df2[[ lhs_vars[1] ]] <= tau_val, ]

      if(nrow(events) == 0) {
        return(0)
      }
      counts <- table(events[[ lhs_vars[2] ]])
      types <- as.numeric(names(counts))

      weighted_counts <- counts * wr[types]

      sum(weighted_counts)
    }, df_long[[ id ]], df_long$X_tau)


    df_long <- df_long %>% mutate(L = N_D + N_R)

    ## ----------------- Create spline covaraites ----------------- ##


    offsets <- c(0, sp_knots[-length(sp_knots)])

    # In df_long, compute helper columns for the current knot index and previous knot.
    df_long <- df_long %>%
      mutate(
        idx = findInterval(tau, sp_knots),
        idx = if_else(tau %in% sp_knots, pmax(idx - 1, 0), idx),
        prev_knot = ifelse(idx == 0, 0, sp_knots[idx]),
        spline_indicator = as.numeric(tau > prev_knot)
      )

    # Determine if the original model includes an intercept.
    has_int <- attr(terms(model), "intercept") == 1


    if (has_int) {
      # Create one intercept spline term per knot.
      for (i in seq_along(sp_knots)) {
        knot_val <- sp_knots[i]
        off <- offsets[i]
        new_name <- paste0("sp_int_", knot_val)
        # For each row: (tau - off) * I(tau > off)
        df_long[[new_name]] <- (df_long$tau - off) * as.numeric(df_long$tau >= off)
      }

      # For each covariate, create spline basis for powers 1 to degree (i.e., degree terms per knot).
      for (v in cov_vars) {
        for (i in seq_along(sp_knots)) {
          knot_val <- sp_knots[i]
          off <- offsets[i]
          for (p in 1:degree) {
            new_name <- paste0(v, "_sp_", knot_val, "_p", p)
            df_long[[new_name]] <- df_long[[v]] * ((df_long$tau - off)^p) * as.numeric(df_long$tau >= off)
          }
        }
      }

    } else {
      # Without an intercept, for each covariate, create spline basis for powers 0 to degree (i.e. degree+1 terms per knot).
      for (v in cov_vars) {
        for (i in seq_along(sp_knots)) {
          knot_val <- sp_knots[i]
          off <- offsets[i]
          for (p in 0:degree) {
            new_name <- paste0(v, "_sp_", knot_val, "_p", p)
            df_long[[new_name]] <- df_long[[v]] * ((df_long$tau - off)^p) * as.numeric(df_long$tau > off)
          }
        }
      }
    }

    df_long <- df_long %>% select(-prev_knot, -spline_indicator)
    pattern <- paste0("^(sp_int_|(", paste(cov_vars, collapse = "|"), ")_sp_)")
    covariate_columns <- grep(pattern, names(df_long), value = TRUE)

    Z_i <- as.matrix(df_long[, covariate_columns])

    if(link == "log"){
      diff <- sum(as.vector(df_long$V_tau * (df_long$L - exp(Z_i %*% beta) )^2))

    }else if(link=="identity"){
      diff <- sum(as.vector(df_long$V_tau * (df_long$L - Z_i %*% beta )^2))

    }

    if(t==1){
      time_diff <- t_span[t]
      int_part <- (diff_pre + 0) * time_diff /2

    }else{
      time_diff <- t_span[t] - t_span[t-1]
      int_part <- (diff_pre + diff) * time_diff /2

    }

    diff_pre <- diff
    PE_sum <- PE_sum + ifelse(is.na(int_part),0,int_part)

  }

  return(PE_sum)
}

#' Cross-Validation Fit for WA Regression (Internal)
#'
#' @description
#' This internal function performs K-fold cross validation to assess the prediction error
#' of the WA regression model. For each fold, the model is fitted on the training data using
#' the specified spline parameters and then the prediction error is computed on the test data.
#' The final cross-validation error is the average prediction error across all folds.
#'
#' @param model A survival model formula (e.g. \code{Surv(time, type) ~ Z1 + Z2}).
#' @param data A data frame containing the data used for fitting the model.
#' @param cluster Optional. A string specifying the column name for cluster ID (if applicable).
#' @param id A string specifying the column name for subject ID.
#' @param wr A numeric vector of weights for recurrent events.
#' @param wd A numeric value representing the weight for the terminal event.
#' @param link Link function for while alive loss rate.
#' @param cens_mod A survival model formula for the censoring distribution (e.g. \code{Surv(time, delta) ~ Z2}).
#' @param degree A non-negative integer specifying the polynomial degree for the spline.
#' @param time_range A numeric vector specifying the spline knot locations.
#' @param L A number of knots for cross-validation.
#' @param K An integer indicating the number of folds for cross-validation (default is 10).
#' @param t_span A numeric vector of time points at which the prediction error is evaluated.
#'
#' @return A numeric value representing the average prediction error across folds.
#'
#' @keywords internal


cv_fit <- function(model, data,cluster=NULL,id,wr,wd,link,
                   cens_mod, degree,time_range,L,
                   K=10,t_span){
  sp_knots <- seq(time_range[1],time_range[2],length=L)

  # Step 1: Identify unique subjects
  unique_cluster <- unique(data[[cluster]])

  # Step 2: Shuffle the unique subjects
  unique_cluster <- sample(unique_cluster)

  # Step 3: Calculate the sizes of each split
  size <- floor(length(unique_cluster) / K)
  remainder <- length(unique_cluster) %% K

  # Initialize an empty list to store each split
  splits <- vector("list", K)

  # Initialize the starting index
  start <- 1

  # Step 4: Loop through each part to create splits based on subjects
  for (i in 1:K) {
    # Calculate the ending index for each split
    end <- start + size - 1 + ifelse(i <= remainder, 1, 0)

    # Get the subjects for the current split
    current_subjects <- unique_cluster[start:end]

    # Subset the original data based on these subjects
    splits[[i]] <- data[data[[cluster]] %in% unique_cluster, ]

    # Update the start index for the next split
    start <- end + 1
  }

  PE_sum <- 0

  for(i in 1:K){
    temp_data <- do.call(rbind, splits[-i])

    temp_fit <- WA_fun(model,temp_data,cluster,id,wr,wd,link,cens_mod,
                       sp_knots,degree,cal_var = F)

    temp_PE <- PE_cal(model,splits[[i]],cluster,id,wr,wd,link,cens_mod,sp_knots,
                      degree,temp_fit$est, t_span)

    PE_sum <- PE_sum + temp_PE
    print(i)

  }

  return(PE_sum)


}








