#' Plot the time-dependent ROC curve with competing risk
#'
#' This function reads in object returned by \code{tdROC.cr()} and plot ROC curve for it.
#'
#' @param x the object returned by \code{tdROC.cr()}.
#' @param lwd user-specified line width. Default is \code{2}.
#' @param xlab user-specified label for x-axis. Default is "\code{1-specificity}".
#' @param ylab user-specified label for y-axis. Default is "\code{sensitivity}".
#' @param xlim user-specified limit for x axis. Default is \code{c(0,1)}.
#' @param ylim user-specified limit for y axis. Default is \code{c(0,1)}.
#' @param main user-specified title for the plot. Default is "\code{ROC curve}"
#' @param col user-specified color for ROC curve. Defualt is "\code{c("red", "blue")}" for the primary event and competing event.
#' @param abline user-specified reference diagnol line. Default is \code{True}.
#' @param \dots for future methods
#'
#' @return Returns several plots of ROC curve. For competing risk data, there are two definitions of controls introduced by Zheng et al, which was listed below
#'
#'  \deqn{
#'  \text{Definition A:} \text{Case} k:T \le \tau, \delta = k; \text{Control}_A: (T>\tau)\cup (T \le \tau \cap \delta \ne k)
#'  }
#'
#'  \deqn{
#'  \text{Definition B:} \text{Case} k:T \le \tau, \delta = k; \text{Control}_B: (T>\tau)
#'  }
#'
#'  For more details about above two definitions, please read details of function \code{tdROC.cr}.
#'  If the \code{tdROC.cr} object comes without bootstrap result, the ROC curve for above two definitions will be plotted together and indicated by the specified \code{col}.
#'  If the \code{tdROC.cr} object with bootstrap result, one more ROC curve with confidence interval will be plotted for each definition.
#'
#' @importFrom graphics plot legend lines
#'
#' @references Zheng Y, Cai T, Jin Y, Feng Z. Evaluating prognostic accuracy of biomarkers under competing risk. Biometrics. 2012;68(2):388-396. doi:10.1111/j.1541-0420.2011.01671.x
#'
#' @examples
#' library(survival)
#' data(Melano)
#' tdROC.cr_res <- tdROC.cr(
#'   X = Melano$thick, Y = Melano$time,
#'   delta = Melano$status, tau = 1800, nboot = 100
#' )
#' plot_tdROC_cr(tdROC.cr_res)
#'
#' @export

plot_tdROC_cr <- function(x, lwd = 2, xlab = "1-specificity", ylab = "sensitivity", xlim = c(0, 1), ylim = c(0, 1), col = c("red", "blue"), main = "ROC curve", abline = T, ...) {
  # Plot the ROC curve estimated by tdROC.cr()
  # Arguments:
  #  -- x: the object returned by tdROC()
  # Return:
  #  -- a plot of ROC curve
  specA <- 1 - x$main_res$ROC$specA
  specB <- 1 - x$main_res$ROC$specB
  sens <- x$main_res$ROC$sens
  tmp1 <- order(specA, sens)
  x1 <- specA[tmp1]
  y1 <- sens[tmp1]
  tmp2 <- order(specB, sens)
  x2 <- specB[tmp2]
  y2 <- sens[tmp2]
  # windows() ;
  plot(
    x = x1, y = y1, xlab = xlab, ylab = ylab,
    type = "s", lwd = lwd, xlim = xlim, ylim = ylim,
    main = main,
    col = col[1]
  )
  lines(
    x = x2, y = y2, type = "s",
    col = col[2]
  )
  if (abline == T) {
    abline(0, 1, col = "gray", lwd = 2, lty = 2)
  }

  legend("bottomright", legend = c(paste0("AUC.A: ", round(x$main_res$AUC.A.integral, 4)), paste0("AUC.B: ", round(x$main_res$AUC.B.integral, 4))), col = col, lty = 1, cex = 1)

  if (inherits(x[[3]], "list") == T) {
    plot(
      x = x1, y = y1, xlab = xlab, ylab = ylab,
      type = "s", lwd = lwd, xlim = xlim, ylim = ylim,
      main = main,
      col = col[1]
    )
    specA <- 1 - x$boot_res$bROC$specA.mean
    specB <- 1 - x$boot_res$bROC$specB.mean

    sens.upper <- x$boot_res$bROC$sens.upper
    sens.lower <- x$boot_res$bROC$sens.lower


    tmp3.u <- order(specA, sens.upper)
    x3.upper <- specA[tmp3.u]
    y3.upper <- sens.upper[tmp3.u]
    tmp3.l <- order(specA, sens.lower)
    x3.lower <- specA[tmp3.l]
    y3.lower <- sens.lower[tmp3.l]

    lines(x3.upper, y3.upper, col = col[1], lty = 2)
    lines(x3.lower, y3.lower, col = col[1], lty = 2)
    legend("bottomright", legend = c(paste0("AUC.A: ", round(x$main_res$AUC.A.integral, 4), "(", round(x$boot_res$bAUC.A.integral$CIlow, 4), ",", round(x$boot_res$bAUC.A.integral$CIhigh, 4), ")")), col = col[1], lty = 1, cex = 1)
    if (abline == T) {
      abline(0, 1, col = "gray", lwd = 2, lty = 2)
    }

    plot(
      x = x2, y = y2, xlab = xlab, ylab = ylab,
      type = "s", lwd = lwd, xlim = xlim, ylim = ylim,
      main = main,
      col = col[2]
    )

    tmp4.u <- order(specB, sens.upper)
    x4.upper <- specB[tmp4.u]
    y4.upper <- sens.upper[tmp4.u]
    tmp4.l <- order(specB, sens.lower)
    x4.lower <- specB[tmp4.l]
    y4.lower <- sens.lower[tmp4.l]

    lines(x4.upper, y4.upper, col = col[2], lty = 2)
    lines(x4.lower, y4.lower, col = col[2], lty = 2)
    legend("bottomright", legend = c(paste0("AUC.B: ", round(x$main_res$AUC.B.integral, 4), "(", round(x$boot_res$bAUC.B.integral$CIlow, 4), ",", round(x$boot_res$bAUC.B.integral$CIhigh, 4), ")")), col = col[2], lty = 1, cex = 1)
    if (abline == T) {
      abline(0, 1, col = "gray", lwd = 2, lty = 2)
    }
  }

  invisible(0)
}
