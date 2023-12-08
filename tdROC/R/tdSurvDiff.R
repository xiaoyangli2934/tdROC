#' Calculate the Survival Difference
#'
#' This function reads in a vector of estimated weight and same length risk score to calculate
#' survival difference by formula (Wu and Li, 2018).This function is used internally by other functions
#' in this package.
#' @param W a numerical vector of weight estimated by nonparametric weight adjustments (Li \emph{et al.}, 2015). Same length with \code{X}.
#' @param X a numerical vector of risk score values. Same length with \code{W}.
#' @return Survival difference as a numerical scalar.
#'
#' @keywords internal

tdSurvDiff <- function(W, X) {
  diff <- ((1 - X) * W + (0 - X) * (1 - W))
  mean_diff <- mean(diff)
  return(mean_diff)
}
