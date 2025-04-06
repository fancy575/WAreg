#' Process WA Regression Results (Internal)
#'
#' @description
#' This internal function wraps the WA regression fitting procedure. It calls
#' \code{WA_fit} to fit the model and, optionally, \code{cal_V} to compute a variance estimate.
#' The function then assembles a list with all relevant inputs and outputs, and returns an object
#' of class \code{"WA_obj"}.
#'
#' @param model A survival model formula (e.g., \code{Surv(time, type) ~ Z1 + Z2}).
#' @param data A data frame containing the data.
#' @param cluster Optional. A string specifying the cluster variable name (if applicable).
#' @param id A string specifying the subject ID column name.
#' @param wr A numeric vector of weights for recurrent events.
#' @param wd A numeric value for the terminal event weight.
#' @param link Link function for while alive loss rate.
#' @param cens_mod A survival model formula for the censoring distribution (e.g., \code{Surv(time, delta) ~ Z2}).
#' @param sp_knots A numeric vector of spline knot locations.
#' @param degree A non-negative integer specifying the polynomial degree for the spline (default is 3).
#' @param cal_var Logical; if \code{TRUE}, the variance estimate is computed using \code{cal_V} (default is \code{FALSE}).
#'
#' @return A list of class \code{"WA_obj"} containing:
#' \describe{
#'   \item{model}{The original model formula.}
#'   \item{data}{The input data.}
#'   \item{cluster}{The cluster variable name (if provided).}
#'   \item{id}{The subject ID column name.}
#'   \item{wr}{The recurrent event weights.}
#'   \item{wd}{The terminal event weight.}
#'   \item{cens_model}{The censoring model formula.}
#'   \item{sp_knots}{The spline knot locations used.}
#'   \item{degree}{The polynomial degree for the spline.}
#'   \item{WA_out}{A list containing the estimated WA outcomes (e.g., estimated coefficients and covariance matrix).}
#'   \item{variance}{The variance estimate (if computed, otherwise \code{NULL}).}
#'   \item{fit}{The fitted object returned by \code{WA_fit}.}
#' }
#'
#' @keywords internal
#'
WA_fun <- function(model, data,cluster=NULL,id,wr,wd,link,
                   cens_mod,
                   sp_knots, degree = 3, cal_var=F){


  lhs_vars <- all.vars(formula(model)[[2]])
  cov_vars <- all.vars(formula(model)[[3]])


  lhs_cens <- all.vars(formula(cens_mod)[[2]])
  cov_cens <- all.vars(formula(cens_mod)[[3]])


  ## -------------------- Prep for EE -------------------- ##

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

  ## stack data based on spline basis


  df_long <- df1 %>%
    tidyr::crossing(tau = sp_knots) %>%
    mutate(
      X_tau = pmin(!!sym(lhs_vars[1]), tau),
      V_tau = if_else(!!sym(lhs_vars[1]) <= tau, delta, 1)
    )

  cens_mod_n <- reformulate(
    attr(terms(cens_mod), "term.labels"),
    response = call("Surv",
                    formula(cens_mod)[[2]][[2]],
                    call("==", as.name("delta"), 0))
  )

  # cens_mod_n <- Surv(time, delta==0) ~ 1
  # cens_mod_n <- Surv(time, delta==0) ~ Z2


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
      knot_index = match(tau, sp_knots),
      prev_knot = sapply(knot_index, function(k) if (k == 1) 0 else sp_knots[k - 1]),
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

  df_long <- df_long %>% select(-knot_index, -prev_knot, -spline_indicator)

  k <- length(sp_knots)
  m <- length(cov_vars)

  if (has_int) {
    n_beta <- m * k * degree + k
  } else {
    n_beta <- m * k * (degree + 1)
  }

  if(!is.null(cluster)){
    feed_dt <- df_long %>%
      select(-all_of(lhs_vars[2])) %>%
      rename(j = !!sym(id),
             i = !!sym(cluster),
             etime = !!sym(lhs_vars[1]))
  }else{
    feed_dt <- df_long %>%
      select(-all_of(lhs_vars[2])) %>%
      rename(j = !!sym(id),
             i = Cluster,
             etime = !!sym(lhs_vars[1]))
  }


  int_beta <- rep(0,n_beta)

  solution <- nleqslv::nleqslv(int_beta, est_eq, fd_dt = feed_dt, cov_vars=cov_vars,link=link)

  beta <- solution$x

  pattern <- paste0("^(sp_int_|(", paste(cov_vars, collapse = "|"), ")_sp_)")
  cov_names <- grep(pattern, names(feed_dt), value = TRUE)
  names(beta) <- cov_names


  ## Calculate naive variance
  if(cal_var==T){
    variance <- tryCatch({
      variance <- cal_V(feed_dt, beta,cov_vars,link)
    }, error = function(e) {
      # Return a 2x2 matrix of zeros if an error occurs
      matrix(0, nrow = length(int_beta), ncol = length(int_beta))
    })

    return(list(est = beta ,se = sqrt(diag(variance)),cov=variance) )

  }else{
    return(list(est = beta) )

  }


}
