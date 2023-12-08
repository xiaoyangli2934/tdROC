#' Example data: Mayo Data
#'
#' This example dataset is included for illustration. The Mayo PBC data is publicly
#'    available and has been used in many statistical researches (e.g., Zheng and Heagerty 2005).
#'    This example dataset is a subset of the full PBC data.
#'
#' @docType data
#' @usage data(mayo)
#' @keywords datasets
#' @format A data frame with 312 observations and 4 variables:
#' @format \code{time}: event time or censoring time
#' @format \code{censor} : censoring indicator.
#' @format \code{mayoscore4} and \code{mayoscore5}: derived from 4 and 5 covariates respectively.
#' @references Heagerty, P. J., & Zheng, Y. (2005). Survival model predictive accuracy and ROC curves. Biometrics, 61(1), 92-105.
#'
#' @examples
#' data(mayo)
#' head(mayo)
"mayo"
