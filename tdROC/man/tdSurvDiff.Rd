% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdSurvDiff.R
\name{tdSurvDiff}
\alias{tdSurvDiff}
\title{Calculate the Survival Difference}
\usage{
tdSurvDiff(W, X)
}
\arguments{
\item{W}{a numerical vector of weight estimated by nonparametric weight adjustments (Li \emph{et al.}, 2015). Same length with \code{X}.}

\item{X}{a numerical vector of risk score values. Same length with \code{W}.}
}
\value{
Survival difference as a numerical scalar.
}
\description{
This function reads in a vector of estimated weight and same length risk score to calculate
survival difference by formula (Wu and Li, 2018).This function is used internally by other functions
in this package.
}
\keyword{internal}
