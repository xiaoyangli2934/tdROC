% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.tdROC.R
\name{plot_tdROC}
\alias{plot_tdROC}
\title{Plot the time-dependent ROC curve}
\usage{
plot_tdROC(
  x,
  lwd = 2,
  xlab = "1-specificity",
  ylab = "sensitivity",
  xlim = c(0, 1),
  ylim = c(0, 1),
  main = "ROC curve",
  col = "black",
  abline = T,
  ...
)
}
\arguments{
\item{x}{the object returned by \code{tdROC()}.}

\item{lwd}{user-specified line width. Default is \code{2}.}

\item{xlab}{user-specified label for x-axis. Default is "\code{1-specificity}".}

\item{ylab}{user-specified label for y-axis. Default is "\code{sensitivity}".}

\item{xlim}{user-specified limit for x axis. Default is \code{c(0,1)}.}

\item{ylim}{user-specified limit for y axis. Default is \code{c(0,1)}.}

\item{main}{user-specified title for the plot. Default is "\code{ROC curve}"}

\item{col}{user-specified color for ROC curve. Defualt is "\code{black}".}

\item{abline}{user-specified reference diagnol line. Default is \code{True}.}

\item{\dots}{for future methods}
}
\value{
Returns a plot of ROC curve. If the tdROC object comes with bootstrap result,
 then the ROC curve will be plotted with confidence interval.
}
\description{
This function reads in object returned by \code{tdROC()} and plot ROC curve for it.
}
\examples{
library(survival)
data(mayo)
dat <- mayo[, c("time", "censor", "mayoscore5")]
fm <- tdROC(
  X = dat$mayoscore5, Y = dat$time, delta = dat$censor,
  tau = 365 * 6, span = 0.1, nboot = 0, alpha = 0.05, n.grid = 1000, cut.off = 5:9
)
# plot the object "fm" from tdROC()
plot_tdROC(fm)

}
