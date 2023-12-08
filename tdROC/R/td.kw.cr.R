#' Calculate conditional probability of being a case at time tau
#'
#' This is key function to estimate the weight, the conditional probability of being a case at time tau
#'  given the observed time to event, event status, and prognostic risk score, as described in Wu and Li, 2018.
#'
#' @param X a numeric vector of risk score for each subject. Higher value of \code{X} indicates higher risk of the event.
#'        It can be biomarker, a function of multiple biomarker, or predicted cumulative incidence function by other methods.
#'        Same length with \code{Y} and \code{delta}.
#' @param Y a numeric vector of time to event. Same length with \code{X} and \code{delta}.
#' @param delta a vector of numeric indicator of event type. The primary event you want to study should be coded as 1,
#'        the competing event should be coded as 2, and censoring should be coded as 0. Same length with \code{X} and \code{Y}.
#' @param event.code numeric indicator of event (1), or competing event (2), it specifies you are going to calculate the
#'        conditional probability for which event.
#' @param tau a scalar, the prediction horizon at which the prediction is evaluated.
#' @param span a numeric value, the proportion of neighbour observations used in nearest neighbor method, default is 0.1.
#' @param h a numeric value, the bandwidth of kernel weights, defualt is \code{NULL}. If not specified, the function will use the value of
#'        \code{span} to calculate kernel weights. In case both \code{span} and \code{h} are specified, the function will use \code{h}.
#' @param type a character value, indicating the type of kernel function used to calculate kernel weights. Default is "\code{uniform}" kernel.
#'        Other options are "\code{Epanechnikov}" and "\code{normal}".
#' @param epsilon the precision parameter for weight calculation using neighborhood approximation. If not specified, default will be
#'        calculating weights for all right censored points individually.
#'
#' @details This function read in the risk score value \code{X}, the time-to-event data \code{Y} and censoring indicator \code{delta}
#'        to estimate the weight, the conditional probability of being a case at time tau when there is competing event.
#'        The weight estimation serves for the further prediction accuracy estimation, including AUC, Brier score and so on.
#' @seealso \code{\link[survival]{survfit}}
#' @return Returns the estimated conditional probability of being a case at time tau for the specified event code.
#' @importFrom survival survfit
#' @keywords internal

td.kw.cr <- function(X, Y, delta, event.code, tau,
                     span = 0.1, h = NULL, type = "uniform", epsilon = 0.01) {
  n <- length(X)
  W <- rep(NA, n)
  delta.cr <- as.factor(delta)

  ## outcome for overall survival
  event.all <- ifelse(delta == 0, 0, 1)

  W[which(Y > tau |
    (Y <= tau & (!delta %in% c(0, event.code))))] <- 0
  W[which(Y <= tau & delta == event.code)] <- 1

  idx <- which(Y <= tau & delta == 0)

  # if epsilon did not specified, calculate all
  if (is.null(epsilon)) {
    for (i in idx) {
      kw <- calc.kw(X = X, x0 = X[i], span = span, h = h, type = type)
      wgt.fit.os <- survfit(Surv(Y, event.all) ~ 1, weights = as.vector(kw))
      wgt.fit.cs <- survfit(Surv(Y, delta.cr) ~ 1, weights = as.vector(kw))
      os.pred <- wgt.fit.os$surv[wgt.fit.os$time == Y[i]]
      # calculate the cumulative incidence estimator from the probability of states
      cif1 <- wgt.fit.cs$pstate[wgt.fit.cs$time == Y[i], event.code + 1]
      time_idx <- which(abs(Y - tau) == sort(abs(Y - tau))[1])[1]
      cif2 <- wgt.fit.cs$pstate[wgt.fit.cs$time == Y[time_idx], event.code + 1]

      W[i] <- ifelse(os.pred == 0, 1, (cif2 - cif1) / os.pred)
    }
  } else { # if epsilon specified
    m <- 0
    grid_points <- vector()
    fit.os_store <- list()
    fit.cs_store <- list()

    for (i in idx) {
      check <- FALSE

      if (length(grid_points) > 0) {
        check <- sum(abs(grid_points - X[i]) / epsilon <= 1) > 0
      }

      if (check == TRUE) {
        # print(check)
        grid_idx <- which.min(abs(grid_points - X[i]) / epsilon)

        os.pred <- fit.os_store[[grid_idx]]$surv[round(fit.os_store[[grid_idx]]$time, 6) == round(Y[i], 6)][1] # round avoid error it did not match
        # calculate the cumulative incidence estimator from the probability of states
        cif1 <- fit.cs_store[[grid_idx]]$pstate[round(fit.os_store[[grid_idx]]$time, 6) == round(Y[i], 6), event.code + 1][1]
        time_idx <- which(abs(Y - tau) == sort(abs(Y - tau))[1])[1]
        cif2 <- fit.cs_store[[grid_idx]]$pstate[round(fit.cs_store[[grid_idx]]$time, 6) == round(Y[time_idx], 6), event.code + 1]
        W[i] <- ifelse(os.pred == 0, 1, (cif2 - cif1) / os.pred)
        # print(i)
      } else { # calculate Wi if it has't been calculated
        # print(check)
        m <- m + 1
        grid_points[m] <- X[i]

        kw <- calc.kw(X = X, x0 = X[i], span = span, h = h, type = type)
        fit.os_store[[m]] <- survfit(Surv(Y, event.all) ~ 1, weights = as.vector(kw))
        fit.cs_store[[m]] <- survfit(Surv(Y, delta.cr) ~ 1, weights = as.vector(kw))

        os.pred <- fit.os_store[[m]]$surv[round(fit.os_store[[m]]$time, 6) == round(Y[i], 6)][1]
        # calculate the cumulative incidence estimator from the probability of states
        cif1 <- fit.cs_store[[m]]$pstate[round(fit.cs_store[[m]]$time, 6) == round(Y[i], 6), event.code + 1][1]
        time_idx <- which(abs(Y - tau) == sort(abs(Y - tau))[1])[1]
        cif2 <- fit.cs_store[[m]]$pstate[round(fit.cs_store[[m]]$time, 6) == round(Y[time_idx], 6), event.code + 1]
        W[i] <- ifelse(os.pred == 0, 1, (cif2 - cif1) / os.pred)
        # print(i)
      }
    }
  }
  return(W)
}
