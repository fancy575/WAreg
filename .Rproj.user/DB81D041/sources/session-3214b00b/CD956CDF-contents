# R/zzz-internal.R -----------------------------------------------------------
# Internal helpers (not exported)

#' Pipe operator
#'
#' Re-exports the `%>%` operator from magrittr so that it is available
#' when WAreg is loaded.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs %>% rhs
NULL


#' @keywords internal
#' @noRd
.WA_link <- function(link = c("log","identity")) {
  link <- match.arg(link)
  if (link == "log") {
    list(link = link,
         linkinv = function(eta) exp(eta),
         dmu_deta = function(mu) mu)
  } else {
    list(link = link,
         linkinv = function(eta) pmax(eta, 0),
         dmu_deta = function(mu) as.numeric(mu >= 0))
  }
}

#' @keywords internal
#' @noRd
.WA_parse_formula <- function(formula, data) {
  tf <- terms(formula, data = data)
  lhs_call  <- attr(tf, "variables")[[2]]
  lhs_vars  <- all.vars(lhs_call)
  if (length(lhs_vars) < 2L)
    stop("Left side must be Surv(time, status).")
  time_name   <- lhs_vars[1]
  status_name <- lhs_vars[2]
  if (!time_name %in% names(data) || !status_name %in% names(data))
    stop("Columns specified in Surv() not found in data.")
  time_vec   <- data[[time_name]]
  status_vec <- data[[status_name]]
  if (!is.integer(status_vec)) {
    if (is.numeric(status_vec)) {
      status_vec <- as.integer(status_vec)
    } else {
      if (any(suppressWarnings(!is.na(as.numeric(levels(as.factor(status_vec))))))) {
        status_vec <- as.integer(as.numeric(status_vec))
      } else {
        status_vec <- as.integer(as.factor(status_vec)) - 1L
      }
    }
  }
  Xmm <- model.matrix(delete.response(tf), data = data)
  if (ncol(Xmm) > 0 && colnames(Xmm)[1] == "(Intercept)") {
    Xmm <- Xmm[, -1, drop = FALSE]
  }
  Z_cols <- colnames(Xmm)
  list(
    time_vec    = time_vec,
    status_vec  = status_vec,
    Xmm         = Xmm,
    Z_cols      = Z_cols,
    time_name   = time_name,
    status_name = status_name,
    tf          = tf
  )
}

#' @keywords internal
#' @noRd
.WA_time_basis <- function(t, knots,
                           basis = c("il","pl","bz","ns","ms","st","tl","tf"),
                           degree = 3,
                           include_intercept = FALSE) {
  basis <- match.arg(basis)
  t  <- as.numeric(t)
  ti <- sort(unique(as.numeric(knots)))
  stopifnot(length(ti) >= 2L)
  tmin <- min(ti); tmax <- max(ti)
  te   <- pmin(pmax(t, tmin), tmax)
  inter <- if (length(ti) > 2L) ti[-c(1L, length(ti))] else NULL
  newM <- function(nc) matrix(0, nrow = length(te), ncol = nc)

  if (basis == "bz") {
    B <- splines::bs(te, degree = degree, knots = inter,
                     Boundary.knots = c(tmin, tmax), intercept = include_intercept)
    tun <- (te - tmin) / pmax(tmax - tmin, .Machine$double.eps)
    out <- sweep(B, 1, tun, `*`)
    colnames(out) <- sprintf("t0_seg%d", seq_len(ncol(out))); return(out)
  }
  if (basis == "ns") {
    B <- splines::ns(te, knots = inter, Boundary.knots = c(tmin, tmax))
    tun <- (te - tmin) / pmax(tmax - tmin, .Machine$double.eps)
    out <- sweep(B, 1, tun, `*`)
    colnames(out) <- sprintf("t0_seg%d", seq_len(ncol(out))); return(out)
  }
  if (basis == "ms") {
    if (!requireNamespace("splines2", quietly = TRUE)) {
      warning("splines2 not installed; fallback to 'bz'.")
      return(.WA_time_basis(t, knots, basis = "bz", degree = degree,
                            include_intercept = include_intercept))
    }
    B <- splines2::mSpline(te, degree = degree, knots = inter,
                           Boundary.knots = c(tmin, tmax), intercept = include_intercept)
    colnames(B) <- sprintf("t0_seg%d", seq_len(ncol(B))); return(B)
  }
  if (basis == "st") {
    R <- length(ti) - 1L
    out <- newM(R); for (r in seq_len(R)) out[, r] <- as.numeric(te >= ti[r])
    colnames(out) <- sprintf("t0_seg%d", seq_len(R)); return(out)
  }
  if (basis == "tl") {
    R <- length(ti) - 1L
    out <- newM(R); for (r in seq_len(R)) out[, r] <- pmax(te - ti[r], 0)
    colnames(out) <- sprintf("t0_seg%d", seq_len(R)); return(out)
  }
  if (basis == "il") {
    R <- length(ti) - 1L
    out <- newM(R)
    for (r in seq_len(R)) {
      kL <- ti[r]; kR <- ti[r+1L]
      sel <- (te > kL) & (te < kR)
      out[sel, r] <- te[sel] - kL
    }
    colnames(out) <- sprintf("t0_seg%d", seq_len(R)); return(out)
  }
  if (basis == "pl") {
    R <- length(ti) - 1L; P <- degree + 1L
    out <- newM(R*P); cn <- character(R*P); col <- 1L
    for (r in seq_len(R)) {
      kL <- ti[r]; kR <- ti[r+1L]
      ind <- (te > kL) & (te <= kR); base <- pmax(te - kL, 0)
      for (p in 0:degree) {
        out[ind, col] <- base[ind]^p
        cn[col] <- sprintf("t%d_seg%d", p, r); col <- col + 1L
      }
    }
    colnames(out) <- cn; return(out)
  }
  if (basis == "tf") {
    out <- matrix(1, nrow = length(te), ncol = 1)
    colnames(out) <- "t0_seg1"; return(out)
  }
  stop("Unknown basis.")
}

#' @keywords internal
#' @noRd
.WA_design_Z <- function(df, Z_cols, knots, basis, degree, include_intercept = FALSE) {
  stopifnot("tau" %in% names(df))
  B <- .WA_time_basis(df$tau, knots = knots, basis = basis,
                      degree = degree, include_intercept = include_intercept)
  for (z in Z_cols) {
    zv <- df[[z]]
    for (j in seq_len(ncol(B))) {
      nm <- sprintf("%s_%s", z, colnames(B)[j])
      df[[nm]] <- zv * B[, j]
    }
  }
  list(
    data = df,
    cov_pattern = paste0("^(", paste(Z_cols, collapse="|"), ")_t[0-9]+_seg[0-9]+$")
  )
}

#' @keywords internal
#' @noRd
.WA_ee <- function(beta, data, cov_pattern,
                   L_col="L", Dmin_col="X_min_tau", V_col="V_i_tau", G_col="G_X_min_tau",
                   link = c("log","identity")) {
  lf <- .WA_link(link)
  Zi_cols <- grep(cov_pattern, names(data), value = TRUE)
  if (length(Zi_cols) != length(beta)) stop("Length(beta) != design columns.")
  Zi  <- as.matrix(data[, Zi_cols, drop = FALSE])
  eta <- as.vector(Zi %*% beta)
  mu  <- lf$linkinv(eta)
  fac  <- data[[V_col]] / pmax(data[[G_col]], .Machine$double.eps)
  diff <- data[[L_col]] - mu * data[[Dmin_col]]
  colMeans(sweep(Zi, 1, fac * diff, `*`))
}

#' @keywords internal
#' @noRd
.WA_var <- function(data, beta, cov_pattern,
                    L_col="L", Dmin_col="X_min_tau", V_col="V_i_tau", G_col="G_X_min_tau",
                    id_col, cluster_col = NULL, link = c("log","identity")) {
  lf <- .WA_link(link)
  Zi_cols <- grep(cov_pattern, names(data), value = TRUE)
  Zi <- as.matrix(data[, Zi_cols, drop = FALSE])
  p  <- ncol(Zi); nR <- nrow(Zi)

  Zb   <- as.vector(Zi %*% beta)
  mu   <- lf$linkinv(Zb)
  fac  <- data[[V_col]] / pmax(data[[G_col]], .Machine$double.eps)
  Xmin <- data[[Dmin_col]]
  Lval <- data[[L_col]]
  adjL <- Lval - mu * Xmin

  wA     <- fac * mu * Xmin
  A      <- (t(Zi) %*% (Zi * wA)) / nR

  tgrid <- sort(unique(data$obs_T))
  m <- length(tgrid); if (m == 0L) stop("No obs_T")
  dt <- c(0, diff(tgrid))

  idx_T   <- findInterval(data$obs_T, tgrid, all.inside = TRUE)
  idx_tau <- findInterval(data$tau,   tgrid, all.inside = TRUE)

  counts_at <- tabulate(idx_T, nbins = m)
  Y_at_risk <- rev(cumsum(rev(counts_at)))

  idx_cens <- idx_T[data$Delta == 0L]
  dN_cens  <- tabulate(idx_cens, nbins = m)

  dLambda_C <- dN_cens / pmax(Y_at_risk, 1L)
  avg_risk  <- Y_at_risk / nR
  inv_avg   <- 1 / pmax(avg_risk, .Machine$double.eps)
  S_cum     <- cumsum(dLambda_C * inv_avg)
  A_i       <- ifelse(data$Delta == 0L, inv_avg[idx_T], 0)

  dK_mat <- (fac * adjL) * Zi
  j_idx  <- match(Xmin, tgrid, nomatch = 0L)

  sums_by_time <- matrix(0, nrow = m, ncol = p)
  if (any(j_idx > 0L)) {
    sums <- rowsum(dK_mat[j_idx > 0L, , drop = FALSE],
                   group = j_idx[j_idx > 0L], reorder = FALSE)
    sums_by_time[as.integer(rownames(sums)), ] <- sums
  }
  dK_by_time <- sums_by_time / nR

  cumtrap <- function(V) {
    Vprev <- rbind(V[1, , drop = FALSE], V[-m, , drop = FALSE])
    apply(dt * (Vprev + V) / 2, 2, cumsum)
  }
  T_dK  <- cumtrap(dK_by_time)
  T_SdK <- cumtrap(S_cum * dK_by_time)

  r    <- idx_tau
  kT   <- idx_T
  rmin <- pmin(r, kT)
  psi1 <- (as.integer(data$Delta == 0L) * inv_avg[kT]) - S_cum[pmin(1L, kT)]

  T_dK0  <- rbind(0, T_dK)
  T_SdK0 <- rbind(0, T_SdK)
  S0     <- c(0, S_cum)
  dt0    <- c(0, dt)

  cl_raw <- if (!is.null(cluster_col) && cluster_col %in% names(data)) data[[cluster_col]] else data[[id_col]]
  cl_id   <- as.integer(factor(cl_raw))
  M       <- length(unique(cl_id))

  cluster_sums  <- matrix(0, nrow = M, ncol = p)
  cluster_sizes <- tabulate(cl_id, nbins = M)

  step <- 200000L
  for (a in seq(1L, nR, by = step)) {
    b <- min(a + step - 1L, nR)
    rr <- r[a:b]; kk <- kT[a:b]; rrmin <- rmin[a:b]; Ai <- A_i[a:b]

    Td_r  <- T_dK0[rr + 1L, , drop = FALSE]
    Td_kT <- T_dK0[kk + 1L, , drop = FALSE]
    TS_r  <- T_SdK0[rrmin + 1L, , drop = FALSE]

    half     <- 0.5 * dt0[kk + 1L]
    dK_kT    <- dK_by_time[kk, , drop = FALSE]
    halfTerm <- sweep(dK_kT, 1L, half, `*`)

    TermA <- sweep(Td_r - Td_kT + halfTerm, 1L, Ai * (rr >= kk), `*`)
    TermB <- TS_r + sweep(Td_r - Td_kT, 1L, S0[kk + 1L] * (rr > kk), `*`)
    integral <- TermA - TermB

    if (any(rr == 1L)) {
      mask1 <- (rr == 1L)
      dK1  <- matrix(dK_by_time[1L, ], nrow = sum(mask1), ncol = p, byrow = TRUE)
      integral[mask1, ] <- psi1[a:b][mask1] * dK1
    }

    Qblock <- dK_mat[a:b, , drop = FALSE] - integral
    sums_block <- rowsum(Qblock, group = cl_id[a:b], reorder = FALSE)
    rows <- as.integer(rownames(sums_block))
    cluster_sums[rows, ] <- cluster_sums[rows, ] + as.matrix(sums_block)
  }

  cluster_means <- sweep(cluster_sums, 1L, pmax(cluster_sizes, 1L), `/`)
  B <- crossprod(cluster_means) / M

  Ainv <- tryCatch(solve(A), error = function(e) {
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Singular A and MASS not installed.")
    MASS::ginv(A)
  })
  Ainv %*% B %*% Ainv / M
}

#' @keywords internal
#' @noRd
.WA_design_at_t <- function(newdata, t, Z_cols, knots, basis, degree, include_intercept=FALSE) {
  B <- .WA_time_basis(t, knots = knots, basis = basis, degree = degree,
                      include_intercept = include_intercept)
  nb <- ncol(B)
  out_cols <- unlist(lapply(Z_cols, function(z) sprintf("%s_t0_seg%d", z, seq_len(nb))), use.names = FALSE)
  X <- matrix(0, nrow = nrow(newdata), ncol = length(out_cols),
              dimnames = list(NULL, out_cols))
  ctr <- 1L
  for (z in Z_cols) {
    blk <- outer(newdata[[z]], B)  # n x nb
    X[, ctr:(ctr+nb-1L)] <- blk
    ctr <- ctr + nb
  }
  X
}

#' @keywords internal
#' @noRd
.WA_ipcw_fit <- function(subj, method = c("km","cox"), ipcw_formula = ~ 1) {
  method <- match.arg(method)
  if (method == "km") {
    sf <- survival::survfit(survival::Surv(time = subj$obs_T, event = 1 - subj$Delta) ~ 1)
    list(method = "km", times = sf$time, surv = sf$surv)
  } else {
    rhs <- paste(deparse(ipcw_formula[[2]]), collapse = "")
    cfit <- survival::coxph(stats::as.formula(paste0("Surv(obs_T, I(1-Delta)) ~ ", rhs)),
                            data = subj, ties = "breslow", x = TRUE, y = TRUE, model = FALSE)
    bh <- survival::basehaz(cfit, centered = FALSE)
    list(method = "cox", cfit = cfit, basehaz = bh, rhs = rhs)
  }
}

#' @keywords internal
#' @noRd
.WA_ipcw_predict_G <- function(ipcw_fit, X_min_tau, newdata = NULL) {
  if (identical(ipcw_fit$method, "km")) {
    stats::approx(ipcw_fit$times, ipcw_fit$surv, xout = X_min_tau,
                  method = "constant", f = 0,
                  yleft = 1, yright = tail(ipcw_fit$surv, 1))$y
  } else {
    bh <- ipcw_fit$basehaz
    Lambda0 <- stats::approxfun(bh$time, bh$hazard, method = "linear", rule = 2)
    mm <- stats::model.matrix(stats::as.formula(paste0("~ -1 + ", ipcw_fit$rhs)), data = newdata)
    beta <- stats::coef(ipcw_fit$cfit)
    if (length(beta)) mm <- mm[, names(beta), drop = FALSE]
    lp <- if (length(beta)) drop(mm %*% beta) else 0
    pmax(exp(-Lambda0(X_min_tau) * exp(lp)), 0)
  }
}
