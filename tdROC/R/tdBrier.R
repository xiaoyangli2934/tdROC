#' Calculate the Brier Score
#'
#' This function reads in a vector of estimated weight and same length biomarker to calculate
#'    Brier Score by formula (Wu and Li, 2018). This function is used internally by other functions
#'    in this package.
#' @param W a numerical vector of weight estimated by nonparametric weight adjustments (Li \emph{et al.}, 2015). Same length with \code{X}.
#' @param X a numerical vector of risk score values. Same length with \code{W}.
#' @return Brier Score as a numerical scalar.
#' @note This function estimates brier score by using the formula from Wu and Li, 2018.
#'
#' @keywords internal
#'
tdBrier <- function(W, X) {
  squared_diff <- ((1 - X)^2 * W + X^2 * (1 - W))
  mean_squared_diff <- mean(squared_diff)
  return(mean_squared_diff)
}
