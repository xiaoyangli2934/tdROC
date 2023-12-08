#' Calculate the area under a ROC curve (AUC) by trapezoidal integration
#'
#' This function reads in a vector of sensitivity and a vector of specificity to calculates
#' the area under the curve (AUC) by trapezoidal integration.
#' @param sens a numerical vector of sensitivity values within the range of (0, 1).
#' @param spec a numerical vector of specificity values within the range of (0, 1).
#' @return It returns AUC as a numerical scalar.
#' @note This function sorts \code{sens} and \code{1-spec} in an increasing order.
#'      A 0 and 1 will be added to the two ends of the sorted vectors. The Area Under the Curve (AUC) is obtained by trapezoidal
#'       integration of the area under the piecewise linear curve obtained by connecting
#'        points in \code{sens} and \code{1-spec}.
#'
#' @keywords internal

AUC_calc_integral <- function(sens, spec) {
  o <- order(sens, 1 - spec)
  y <- sens[o]
  x <- 1 - spec[o]
  x <- c(0, x, 1)
  y <- c(0, y, 1)
  m <- length(x)
  x1 <- x[-m]
  x2 <- x[-1]
  y1 <- y[-m]
  y2 <- y[-1]
  AUC <- sum((y1 + y2) * (x2 - x1) / 2)
  AUC
}
