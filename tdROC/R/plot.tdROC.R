#' Plot the time-dependent ROC curve
#'
#' This function reads in object returned by \code{tdROC()} and plot ROC curve for it.
#'
#' @param x the object returned by \code{tdROC()}.
#' @param lwd user-specified line width. Default is \code{2}.
#' @param xlab user-specified label for x-axis. Default is "\code{1-specificity}".
#' @param ylab user-specified label for y-axis. Default is "\code{sensitivity}".
#' @param xlim user-specified limit for x axis. Default is \code{c(0,1)}.
#' @param ylim user-specified limit for y axis. Default is \code{c(0,1)}.
#' @param main user-specified title for the plot. Default is "\code{ROC curve}"
#' @param col user-specified color for ROC curve. Defualt is "\code{black}".
#' @param abline user-specified reference diagnol line. Default is \code{True}.
#' @param \dots for future methods
#'
#' @return Returns a plot of ROC curve. If the tdROC object comes with bootstrap result,
#'  then the ROC curve will be plotted with confidence interval.
#'
#' @importFrom graphics plot legend lines
#'
#' @examples
#' library(survival)
#' data(mayo)
#' dat <- mayo[, c("time", "censor", "mayoscore5")]
#' fm <- tdROC(
#'   X = dat$mayoscore5, Y = dat$time, delta = dat$censor,
#'   tau = 365 * 6, span = 0.1, nboot = 0, alpha = 0.05, n.grid = 1000, cut.off = 5:9
#' )
#' # plot the object "fm" from tdROC()
#' plot_tdROC(fm)
#'
#' @export

plot_tdROC <- function(x, lwd = 2, xlab = "1-specificity", ylab = "sensitivity", xlim = c(0, 1), ylim = c(0, 1), main = "ROC curve", col = "black", abline = T, ...) {
  # Plot the ROC curve estimated by tdROC()
  # Arguments:
  #  -- x: the object returned by tdROC()
  # Return:
  #  -- a plot of ROC curve
  spec <- 1 - x$main_res$ROC$spec
  sens <- x$main_res$ROC$sens
  tmp <- order(spec, sens)
  x1 <- spec[tmp]
  y1 <- sens[tmp]
  plot(
    x = x1, y = y1, xlab = xlab, ylab = ylab,
    type = "s", lwd = lwd, xlim = xlim, ylim = ylim,
    main = main,
    col = col
  )
  if (inherits(x[[3]], "list") == FALSE) {
    legend("bottomright", legend = paste0("AUC: ", round(x$main_res$AUC.integral, 4)), col = col, lty = 1, cex = 1)
  } else {
    spec <- 1 - x$boot_res$bROC$spec.mean

    sens.upper <- x$boot_res$bROC$sens.upper
    sens.lower <- x$boot_res$bROC$sens.lower

    tmp2.u <- order(spec, sens.upper)
    x2.upper <- spec[tmp2.u]
    y2.upper <- sens.upper[tmp2.u]
    tmp2.l <- order(spec, sens.lower)
    x2.lower <- spec[tmp2.l]
    y2.lower <- sens.lower[tmp2.l]

    lines(x2.upper, y2.upper, col = col[1], lty = 2)
    lines(x2.lower, y2.lower, col = col[1], lty = 2)
    legend("bottomright", legend = c(paste0("AUC: ", round(x$main_res$AUC.integral, 4), "(", round(x$boot_res$bAUC$CIlow, 4), ",", round(x$boot_res$bAUC$CIhigh, 4), ")")), col = col[1], lty = 1, cex = 1)
  }
  if (abline == T) {
    abline(0, 1, col = "gray", lwd = 2, lty = 2)
  }
}
