#' Example Dataset for WA Regression
#'
#' @description
#' This dataset \code{ex_dt} is a small example demonstrating a structure typically used by
#' \code{WA_reg} functions. It contains columns for cluster ID, subject ID, time, event type,
#' and two covariates (Z1 and Z2).
#'
#' @format A data frame with 7 rows and 6 columns:
#' \describe{
#'   \item{Cluster}{Integer. The cluster ID (if applicable).}
#'   \item{Subject}{Integer. The subject/patient ID.}
#'   \item{time}{Numeric. The observed event time.}
#'   \item{type}{Integer. The event type (e.g., 1 = recurrent event, 2 = terminal event).}
#'   \item{Z1}{Integer. A binary covariate (for illustration).}
#'   \item{Z2}{Numeric. A continuous covariate.}
#' }
#'
#' @usage data(ex_dt)
#'
#' @examples
#' \dontrun{
#' data(ex_dt)
#' head(ex_dt)
#' }
#'
"ex_dt"
