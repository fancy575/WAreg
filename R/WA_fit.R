#' Fit Weighted Composite Endpoint Model with While-Alive Frequency Measures
#'
#' @description
#' This function fits a regression model for weighted composite endpoints using
#' while-alive cumulative frequency measures. It supports both individual- and
#' cluster-randomized trial settings and accommodates time-varying effects via splines.
#' Optionally performs cross-validation to select optimal spline complexity.
#' @import dplyr
#' @import survival
#' @import nleqslv
#' @import ggplot2
#' @param model A survival formula, e.g., \code{Surv(time, type) ~ Z1 + Z2}.
#' @param data A data frame containing the necessary columns for analysis.
#' @param cluster Optional string. Name of the column indicating cluster ID (for cluster-randomized trials), if NULL, it indicates independent data.
#' @param id A string. Name of the column indicating subject ID.
#' @param wr A numeric vector of weights for recurrent event types (excluding terminal event).
#' @param wd A single numeric value. Weight for terminal event.
#' @param link A link function for the while-alive loss rate ,could be "log","identity" (default is \code{"log"})
#' @param cens_mod A survival model formula for the censoring distribution (e.g., \code{Surv(time, delta) ~ Z1}).
#' @param sp_knots A numeric vector specifying the internal knots for spline basis (e.g., \code{c(5, 10, 15)}).
#' @param degree A non-negative integer specifying the spline polynomial degree (default is 1).
#' @param cv Logical. Whether to perform cross-validation to select number of knots (default is \code{FALSE}).
#' @param K Integer. Number of folds for cross-validation (default is 10).
#' @param tr A numeric vector of length 2 giving the time range for cross-validation (\code{c(start, end)}).
#' @param nk Integer vector. Candidate numbers of knots to consider during cross-validation (default is \code{2:6}).
#' @param ts A numeric vector of time points at which to evaluate treatment/covariate effects (default is \code{seq(0.001, 40, length = 1000)}).
#'
#' @return A list with class \code{"WA_obj"} containing:
#' \describe{
#'   \item{model}{Original model formula.}
#'   \item{cens_model}{Censoring model formula.}
#'   \item{degree}{Spline degree used.}
#'   \item{sp_knots}{Spline knot locations used.}
#'   \item{cov_vars}{Covariate names extracted from model.}
#'   \item{WA_out}{A list containing estimated coefficients, covariance matrix, and influence function.}
#'   \item{results_df}{A data frame of estimated time-varying effects and confidence intervals.}
#'   \item{call}{The original function call.}
#' }
#'
#' @examples
#' \dontrun{
#' data(ex_dt)
#' fit <- WA_fit(Surv(time, type) ~ Z1 + Z2, data = ex_dt,
#'               id = "Subject", cluster = "Cluster",
#'               wr = c(1, 1), wd = 1,
#'               cens_mod = Surv(time, delta) ~ Z2,
#'               sp_knots = c(5, 10, 15, 20, 25, 30),
#'               degree = 1)
#' summary(fit)
#' plot(fit, ts = seq(0, 40, by = 1))
#' }
#'
#' @export

WA_fit <- function(model,data,cluster=NULL,id,wr,wd,link="log",
                   cens_mod,
                   sp_knots, degree = 1,cv=F,K=10,tr = c(5,35),
                   nk = 2:6,ts = seq(0.001,40,length=1000)){


  # data <- ex_dt
  # model <- Surv(time,type)~Z1+Z2
  # cluster <- "Cluster"
  # wr <- c(1,1)
  # wd <- 1
  # sp_knots=c(5,10,15,20,25,30,35,40)
  # degree = 1
  # id = "Subject"
  # cens_mod <- Surv(time,type)~Z2
  ## ------------------------- Diagnosis ------------------------- ##

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  na_rows <- sum(!complete.cases(data))
  if(na_rows > 0){
    message("Found ", na_rows, " rows with missing values. These rows will be removed.")
    data <- data[complete.cases(data), ]
  }

  lhs_vars <- all.vars(formula(model)[[2]])
  na_lhs <- setdiff(lhs_vars, names(data))
  if (length(na_lhs) > 0) {
    stop("The following outcomes in the model are missing from data: ",
         paste(na_lhs, collapse = ", "))
  }

  # Extract covariate names from the right-hand side (e.g. "Z1+Z2")
  cov_vars <- all.vars(formula(model)[[3]])
  na_cov <- setdiff(cov_vars, names(data))
  if (length(na_cov) > 0) {
    stop("The following covariates are missing from data: ",
         paste(na_cov, collapse = ", "))
  }


  ## Check censoring model

  lhs_cens <- all.vars(formula(cens_mod)[[2]])
  na_lhs_c <- setdiff(lhs_cens, names(data))
  if (length(na_lhs) > 0) {
    stop("The following outcomes in the censoring model are missing from data: ",
         paste(na_lhs_c, collapse = ", "))
  }

  cov_cens <- all.vars(formula(cens_mod)[[3]])
  if(!is.null(cov_cens)){
    na_cov_c <- setdiff(cov_cens, names(data))
    if (length(na_cov_c) > 0) {
      stop("The following covariates in censoring model are missing from data: ",
           paste(na_cov_c, collapse = ", "))
    }

  }


  if (!is.null(cov_cens)) {
    if (!(cov_cens %in% names(data))) {
      stop("Covariates in censoring model '", cens_cov, "' not found in data")
    }
  }

  # Check wr length: using lhs_vars[2] (i.e. "type")
  if (length(wr) != (max(data[[ lhs_vars[2] ]], na.rm = TRUE) - 1)) {
    stop("Length of 'wr' must equal to (max(type) - 1)")
  }

  if (!link %in% c("log", "identity")) {
    stop("Error: Unsupported link function. Only 'log' and 'identity' are allowed.")
  }



  # Check if id exists in data
  if (!(id %in% names(data))) {
    stop("ID column '", id, "' not found in data")
  }

  # Check if cluster exists, if provided
  if (!is.null(cluster)) {
    if (!(cluster %in% names(data))) {
      stop("Cluster column '", cluster, "' not found in data")
    }
  }else{
    data[["Cluster"]] <- data[[id]]
  }


  # Check that degree is numeric and >= 0
  if (!is.numeric(degree) || degree < 0) {
    stop("'degree' must be numeric and greater than or equal to 0")
  }

  # Check that the maximum sp_knot is smaller than the maximum event time
  if (max(sp_knots) >= max(data[[ lhs_vars[1] ]], na.rm = TRUE)) {
    stop("The maximum value in knots must be smaller than the maximum event time")
  }


  # Check that the maximum sp_knot is smaller than the maximum event time
  if (max(tr) >= max(data[[lhs_vars[1] ]], na.rm = TRUE)) {
    stop("The maximum value in knots seletion must be smaller than the maximum event time")
  }

  # Clean up temporary diagnosis variables
  invisible({
    rm(na_cov, na_lhs, na_cov_c, na_lhs_c)
    gc()
  })
  ## -------------------- Prep for EE -------------------- ##

  if(cv==T){
    cv_val <- c()
    for(k in length(nk)){
      cv <- cv_fit(model,data,cluster,id,wr,wd,link,cens_mod,degree,
                   tr,nk[k],K,ts)

      cv_val <- c(cv_val,cv)
      message("Number of knots: ", nk[k], ", prediction error: ", cv)

    }

    cv_sum <- data.frame(knots_num = nk, pe = cv_val)

    op <- nk[which.min(cv_val)]
    sp_knots = seq(tr[1],tr[2],length=op)

  }else{

    WA_reg <- WA_fun(model, data,cluster,id,wr,wd,link,
                     cens_mod,
                     sp_knots, degree,cal_var=T)

  }

  result <- list(WA_out=WA_reg, data = data,cens_mod=cens_mod,model=model,
                 sp_knots=sp_knots,degree=degree,link=link)

  class(result) <- "WA_obj"

  return(result)

}

#' Summary Method for WA_obj Objects
#'
#' @description
#' This function provides a summary of a WA_obj object including the spline knots,
#' the primary model, the censoring model, the degree of the spline, and a summary
#' table of the estimated parameters (along with their standard errors, z-values,
#' and two-sided p-values).
#'
#' @param object A WA_obj object returned by WA_reg.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns the summary table.
#' @export
summary.WA_obj <- function(object, ...) {
  cat("--------------------------------------------------\n")
  cat("Summary of WA_obj\n")
  cat("--------------------------------------------------\n")

  cat("\nSpline Knots:\n")
  print(object$sp_knots)

  cat("\nModel:\n")
  print(object$model)

  cat("\nCensoring Model:\n")
  print(object$cens_model)

  cat("\nSpline Degree:\n")
  print(object$degree)

  if (!is.null(object$WA_out)) {
    # Extract estimates and compute standard errors.
    est <- object$WA_out$est
    se <- sqrt(diag(object$WA_out$cov))
    z_val <- est / se
    p_val <- 2 * (1 - pnorm(abs(z_val)))

    summary_table <- data.frame(
      Estimate = est,
      SE = se,
      z = z_val,
      p_value = p_val
    )
    rownames(summary_table) <- names(est)

    cat("\nSummary Table for WA_out:\n")
    print(summary_table)
  } else {
    cat("\nNo WA_out information available.\n")
    summary_table <- NULL
  }

  invisible(summary_table)
}


#' Summary Method for WA_obj Objects
#'
#' @description
#' This function provides a summary of a WA_obj object including the spline knots,
#' the primary model, the censoring model, the degree of the spline, and a summary
#' table of the estimated parameters (along with their standard errors, z-values,
#' and two-sided p-values).
#'
#' @param object A WA_obj object returned by WA_reg.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns the summary table.
#' @export

summary.WA_obj <- function(object){
  cat("--------------------------------------------------\n")
  cat("Summary of While-alive regression\n")
  cat("--------------------------------------------------\n")

  cat("\nSpline Knots:\n")
  print(object$sp_knots)

  cat("\nModel:\n")
  print(object$model)

  cat("\nCensoring Model:\n")
  print(object$cens_mod)

  cat("\nDegree:\n")
  print(object$degree)

  # Create summary table from WA_out, if present.
  if (!is.null(object$WA_out)) {
    est <- object$WA_out$est
    se <- object$WA_out$se
    z_val <- est / se
    p_val <- 2 * (1 - pnorm(abs(z_val)))

    summary_table <- data.frame(
      Estimate = est,
      SE = se,
      z = z_val,
      p_value = p_val
    )

    # Original row names might be like "sp_int_5" or "Z1_sp_5_p1".
    original_names <- rownames(summary_table)
    sp_knots <- object$sp_knots

    new_names <- sapply(original_names, function(x) {
      if (grepl("^sp_int_", x)) {
        # Extract the numeric part after "sp_int_"
        num <- as.numeric(sub("sp_int_([0-9]+).*", "\\1", x))
        # Find the index of this knot in sp_knots.
        knot_index <- which(sp_knots == num)
        paste("Intercept", knot_index)
      } else if (grepl("^[A-Za-z0-9]+_sp_", x)) {
        # For example, "Z1_sp_5_p1"
        cov <- sub("^([A-Za-z0-9]+)_sp_.*", "\\1", x)
        knot_val <- as.numeric(sub("^[A-Za-z0-9]+_sp_([0-9]+)_p.*", "\\1", x))
        power <- sub(".*_p([0-9]+)$", "\\1", x)
        knot_index <- which(sp_knots == knot_val)
        paste(cov, "knot", knot_index, "degree", power)
      } else {
        x
      }
    })

    rownames(summary_table) <- new_names


    cat("\nSummary Table for While-alive regression results:\n")
    print(summary_table)
  } else {
    cat("\nNo WA_out information available.\n")
  }

  invisible(object)


}

#' Build Spline Basis Vector
#'
#' @description
#' Constructs a spline basis vector evaluated at time \code{t} based on the provided spline knots,
#' polynomial degree, and intercept flag. When \code{intercept} is \code{TRUE}, the function returns
#' only the treatment-effect portion (i.e. powers 1:\code{degree} for each interval). When \code{intercept} is \code{FALSE},
#' it returns the full basis (i.e. powers 0:\code{degree} for each interval).
#'
#' @param t A numeric value specifying the time at which to evaluate the basis.
#' @param sp_knots A numeric vector of knot locations.
#' @param degree A non-negative integer specifying the polynomial degree.
#' @param intercept Logical; if \code{TRUE}, assumes an intercept was included in the model.
#'
#' @return A numeric vector representing the spline basis evaluated at time \code{t}.
#'
#' @keywords internal

build_L <- function(t, sp_knots, degree, intercept) {
  k <- length(sp_knots)

  time_intervals <- c(0, sp_knots[-k])

  if (intercept) {
    L_cov <- unlist(lapply(time_intervals, function(lower) {
      if (t > lower) sapply(1:degree, function(p) (t - lower)^p) else rep(0, degree)
    }))
    return(c(L_cov))
  } else {
    L_all <- unlist(lapply(time_intervals, function(lower) {
      if (t > lower) sapply(0:degree, function(p) (t - lower)^p) else rep(0, degree)
    }))
    return(L_all)
  }
}

#' Plot WA Object
#'
#' @description
#' This function produces plots of the estimated time-varying covariate effects contained
#' in a WA object. It first plots the estimates (with confidence bands). The plots are faceted by the variables.
#'
#' @param object A WA object (of class \code{WA_obj}) containing a \code{results_df} element.
#' @param ts A numeric vector of time points at which the effect is evaluated.
#' @param alpha A numeric value indicating the significance level for confidence bands (default is 0.05).
#'
#' @return A list with two elements: \code{raw}, the raw results data frame, and \code{smoothed},
#' the loess-smoothed results data frame.
#'
#' @examples
#' \dontrun{
#' # Assuming wa_obj is an object returned by WA_reg with a results_df component:
#' ts <- seq(0, 40, by = 1)
#' plot_results <- plot.WA_obj(wa_obj, ts)
#' }
#'
#' @export


plot.WA_obj <- function(object, ts, alpha = 0.05) {
  library(ggplot2)
  library(dplyr)

  # In this version, we assume that object is already of class "WA_obj"
  # and that object$WA_out contains the estimated parameters and covariance,
  # and that the plot data frame is constructed to include:
  #   t, beta_arm, var_arm, lb, ub, and arm (a grouping variable, e.g. "0", "1").
  # For demonstration, assume that object$results_df has been created.

  # Here we simulate a results data frame for illustration.
  # In your actual code, results_df is created by your WA object.
  results_df <- object$results_df

  # First, plot the raw (unsmoothed) fitted curve.
  p_raw <- ggplot(results_df, aes(x = t, y = beta_arm, color = factor(arm), fill = factor(arm))) +
    geom_line(size = 1) +                    # Plot the fitted curve
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +  # Add confidence band
    labs(x = "Time (t)", y = "Beta Arm", color = "Arm", fill = "Arm") +
    theme_minimal() +
    theme(legend.position = "top")

  # Next, create a smoothed version of the curve using loess.
  smoothed_df <- results_df %>%
    group_by(arm) %>%
    mutate(
      beta_smooth = predict(loess(beta_arm ~ t, span = 0.3)),
      lb_smooth = predict(loess(lb ~ t, span = 0.3)),
      ub_smooth = predict(loess(ub ~ t, span = 0.3))
    ) %>%
    ungroup()

  p_smooth <- ggplot(smoothed_df, aes(x = t, y = beta_smooth, color = factor(arm), fill = factor(arm))) +
    geom_line(size = 0.5) +  # Smoothed fitted curve
    geom_ribbon(aes(ymin = lb_smooth, ymax = ub_smooth), alpha = 0.2, color = NA) +  # Smoothed confidence band
    labs(x = "Time (month)", y = "Treatment effect", color = "Treatment group", fill = "Treatment group") +
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )

  print(p_smooth)

  # Return a list with both the raw and smoothed results.
  return(list(results_df))
}


plot.WA_obj <- function(object, ts, alpha = 0.05, span = 0.3) {

  sp_knots <- object$sp_knots
  degree   <- object$degree
  cov_vars <- all.vars(formula(object$model)[[3]])
  k <- length(sp_knots)
  m <- length(cov_vars)

  # Determine whether the original model included an intercept.
  has_int <- attr(terms(object$model), "intercept") == 1

  # Initialize a list to store the results for each covariate.
  results_list <- vector("list", m)

  for (j in 1:m) {
    if (has_int) {
      start_idx <- k + (j - 1) * k * degree + 1
      end_idx <- k + j * k * degree
    } else {
      # When no intercept, each covariate's block is of length k*(degree+1).
      start_idx <- (j - 1) * k * (degree + 1) + 1
      end_idx <- j * k * (degree + 1)
    }
    beta_cov <- object$WA_out$est[start_idx:end_idx]
    cov_cov <- object$WA_out$cov[start_idx:end_idx, start_idx:end_idx]

    fitted_params <- numeric(length(ts))
    variances <- numeric(length(ts))

    for (i in seq_along(ts)) {
      t_val <- ts[i]
      # build_L returns the covariate spline basis (ignoring intercept) when has_int == TRUE.
      L_effect <- build_L(t_val, sp_knots, degree, intercept = has_int)
      fitted_params[i] <- sum(L_effect * beta_cov)
      variances[i] <- as.numeric(t(L_effect) %*% cov_cov %*% L_effect)
    }

    df_temp <- data.frame(t = ts,
                          beta_effect = fitted_params,
                          var_effect = variances)
    df_temp$lb <- df_temp$beta_effect - qnorm(1 - alpha/2) * sqrt(df_temp$var_effect)
    df_temp$ub <- df_temp$beta_effect + qnorm(1 - alpha/2) * sqrt(df_temp$var_effect)
    df_temp$covariate <- cov_vars[j]
    results_list[[j]] <- df_temp
  }

  # Combine the results for all covariates.
  results_df <- do.call(rbind, results_list)


  # Smooth the results by covariate using loess.
  smoothed_df <- results_df %>%
    group_by(covariate) %>%
    mutate(
      beta_smooth = predict(loess(beta_effect ~ t, span = span)),
      lb_smooth = predict(loess(lb ~ t, span = span)),
      ub_smooth = predict(loess(ub ~ t, span = span))
    ) %>%
    ungroup()

  # Plot the smoothed results.
  p_smooth <- ggplot(smoothed_df, aes(x = t, y = beta_smooth, color = covariate, fill = covariate)) +
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin = lb_smooth, ymax = ub_smooth), alpha = 0.2, color = NA) +
    facet_wrap(~ covariate, scales = "free_y") +
    labs(x = "Time (month)", y = "Covariate Effect") +
    theme_minimal() +
    theme(legend.position = "none")

  print(p_smooth)

  return(list(results_df,plot=p_smooth))
}
