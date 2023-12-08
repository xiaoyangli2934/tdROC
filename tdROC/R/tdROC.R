#' Estimate time-dependent prediction accuracy measures, including the ROC, AUC, Brier score,
#'    and survival difference, with right-censored survival data.
#'
#' This is a core function of the ‘tdROC‘ package. It uses the nonparametric weights proposed by Li (Li et al., 2015)
#' to estimate a number of time-dependent prediction accuracy measures for right-censored survival outcomes,
#' including ROC curve, AUC, Brier score, and survival difference. For each measure, the variance can be estimated
#' through bootstrap resampling.
#'
#' @param X  a numeric vector of risk score in the same length as \code{Y} and \code{delta}, such as biomarker or predicted probability. A higher value indicates higher risk of the event.
#'  The calibration results (Brier score, survival difference) are applicable only when the risk score has the predicted probability interpretation.
#' @param Y a numeric vector of time to event in the same length as \code{X} and \code{delta}.
#' @param delta a vector of binary indicator of event (1) or censoring (0) in the same length as \code{X} and \code{Y}.
#' @param tau a scalar, the prediction horizon at which the prediction is evaluated.
#' @param span a numeric value, the proportion of neighbour observations used in nearest neighbor method, default to 0.1.
#' @param h a numeric value, the bandwidth of kernel weights, the defualt is \code{NULL}. If not specified, the function will use the value of
#'  \code{span} to calculate kernel weights. In case both \code{span} and \code{h} are specified, the function will use \code{h}.
#' @param type a character value, indicating the type of kernel function used to calculate kernel weights. The default is \code{"uniform"} kernel. Other options are \code{"Epanechnikov"} and \code{"normal"}.
#'  It will only be used when the bandwidth \code{h} is specified.
#' @param n.grid a positive integer, the number of grid points used when calculating the ROC curve. The default is \code{1000}.
#' @param X.min the lower boundary of grid cut-off points for biomarker \code{X}. The default is the minimum of \code{X}.
#' @param X.max the upper boundary of grid cut-off points for biomarker \code{X}. The default is the maximum of \code{X}.
#' @param cut.off a vector of \code{X} cut-off values at which sensitivity and specificity will be calculated.
#' @param nboot the number of bootstrap replications to be used for variance estimation. The default is \code{nboot = 0}, corresponding to no variance estimation.
#' @param alpha It is (1 - level of confidence interval)/2, default is \code{0.05}. It is used only when \code{nboot > 0}.
#' @param epsilon The precision parameter used in an approximation to the weight calculation when the sample size is large. If a weight corresponding to a specific risk score is already calculated,
#'   then the weights corresponding to adjacent risk scores, within the distance specified by epsilon, will be the same under the approximation. This approximation avoids repeated
#'   calculation of weights that are almost the same, and hence increases the speed of computation in this situation. The default is NULL, which means no approximation is used. A higher
#'   value indicates less precision.
#'
#' @param method It is used to specify which method you would like to use to estimate AUC, default to \code{"both"}. Other options are \code{"integral"} and \code{"empirical"}.
#' @param output It is used to specify which kind of output you want, default to \code{"both"}. Other options are \code{"AUC"}, including AUC, sensitivity, and specificity are included,
#'  and \code{"calibration"} including Brier Score and survival difference.
#'
#' @details This function takes the risk score value \code{X}, the time-to-event data \code{Y} and censoring indicator \code{delta} as input to estimate
#'  a number of time-dependent prediction accuracy measures for right-censored survival outcomes, including ROC curve, AUC, Brier score, and survival difference.
#'  The confidence intervals of above quantities will be estimated by bootstrap.
#'
#'  This function offer two options to estimate AUC. The first one make use of estimated sensitivity and specificity to calculate the AUC via trapezoidal integration
#'  by setting a series of cutoff point. The output will also include corresponding sensitivity and specificity for our plot function. The other one estimate AUC by the empirical estimator
#'  of the proportion of concordance pairs with proposed weight estimator (Li et al, 2015). These two methods will generate quite similar estimates. The option can be set by argument \code{method}.
#'
#'  We also include Brier Score and survival difference to evaluate the calibration metrics. Their definitions are included below.
#'  They can be estimated with the proposed conditional probability weight (Wu and Li, 2018).
#'  Both of them are measures to assess the accuracy of probabilistic predictions \code{X}. The calibration result makes sense only
#'  when the risk score \code{X} is a predicted probability, and should be ignored otherwise.
#'
#'  \deqn{
#'  \text{Brier Score} = E{[1(T \le \tau, \delta = 1) - X]^2}
#'  }
#'
#'  \deqn{
#'  \text{Survival difference} = E[1(T \le \tau, \delta = 1) - X]
#'  }
#'
#'  As mentioned in arguments, we introduced a small precision parameter \code{epsilon} to speed up the computation when the sample size is large.
#'  For each subject with a risk score, \eqn{X_i}, we assess whether there exists a previously processed grid point, \eqn{X_{grid,m}} where \eqn{1\le m \le j},
#'  within the proximity of \eqn{X_i} such that \eqn{|X_i - X_{grid,m}| < \epsilon}. In the absence of such a point, we designate \eqn{X_i} as a new grid point,
#'  \eqn{X_{grid,j+1}}, and store the corresponding \code{survfit} object for subsequent weight estimation and mark it as a processed grid point. Conversely,
#'  if a previously processed grid point is found, we directly utilize the stored \code{survfit} object associated with it for weight calculation.
#'  Given that the most time-consuming step in our estimation process is the \code{survfit} computation, this method significantly reduces computing time
#'  without incurring notable bias especially when dealing with large sample sizes.
#'
#' @return Returns a list of the following items:
#
#' @return \code{main_res:} a list of \code{AUC.integral} estimated by trapezoidal integration, \code{AUC.empirical} estimated by empirical estimator of the proportion of concordance pairs.
#'  and a data frame \code{ROC} with dimension \code{(2+n.grid) x 3} with columns \code{cut.off}, \code{sens}, and \code{spec}.
#' @return \code{calibration_res:} brier score and survival difference estimated based on the formula similar to Wu and Li (2018). When the risk score \code{X}
#' is a biomarker value instead of a predicted cumulative incidence probability, the brier score and survival difference cannot be calculated. In this case, please disregard the calibration results.
#' @return \code{boot_res:} a list of bootstrap results, including \code{bAUC}, \code{bAUC2}, \code{bBS}, \code{bSurvDiff}, \code{bROC}.
#'  For \code{bAUC}, \code{bAUC2}, \code{bBS}, \code{bSurvDiff}, each one is a list including corresponding mean, standard deviation, and confidence interval.
#'  \code{bROC} is a data frame with colomns \code{sens.mean}, \code{sens.sd}, \code{sens.lower}, \code{sens.upper}, \code{spec.mean}, \code{spec.sd}, \code{spec.lower}, \code{spec.upper}

#' @importFrom survival survfit
#' @importFrom stats sd quantile
#'
#' @examples
#' library(survival)
#' data(mayo)
#' dat <- mayo[, c("time", "censor", "mayoscore5")]
#' fm <- tdROC(
#'   X = dat$mayoscore5, Y = dat$time, delta = dat$censor,
#'   tau = 365 * 6, span = 0.1, nboot = 0, alpha = 0.05,
#'   n.grid = 1000, cut.off = 5:9
#' )
#' # In the following example, We use biomarker mayoscore5 to estimate predicted probability
#' # tipycally a monotone transformation function such as expit() is used to transform biomarker
#' # with range out of range into estimated probability between 0 and 1
#' expit <- function(x){ 1/(1+exp(-x)) }
#'
#' tdROC(
#'   X = expit(dat$mayoscore5), Y = dat$time, delta = dat$censor,
#'   tau = 365 * 6, span = 0.1, nboot = 0, alpha = 0.05,
#'   n.grid = 1000, cut.off = 5:9
#' )
#'
#' tdROC(
#'   X = expit(dat$mayoscore5), Y = dat$time, delta = dat$censor,
#'   tau = 365 * 6, span = 0.1, nboot = 0, alpha = 0.05,
#'   n.grid = 1000, cut.off = 5:9, epsilon = 0.05
#' )
#'
#'
#' @export

tdROC <- function(X, Y, delta, tau, span = 0.1, h = NULL, type = "uniform",
                  n.grid = 1000, X.min = NULL, X.max = NULL, cut.off = NULL,
                  nboot = 0, alpha = 0.05, epsilon = NULL,
                  method = "both", output = "both") {

  n <- length(X)
  positive <- rep(NA, n)
  if (!is.null(epsilon)) { # If epsilon specified, start to use approximation method
    # Shuffle all points
    # set.seed(123)
    idx_random <- sample(n)
    X <- X[idx_random]
    Y <- Y[idx_random]
    delta <- delta[idx_random]
    # Initialize list object to save fitted survfit() results
    res_times <- res_survs <- list()
    grid_points <- c()
    for (i in 1:n) {
      if (Y[i] > tau) {
        positive[i] <- 0
      } else {
        if (delta[i] == 1) {
          positive[i] <- 1
        } else { # In the case of right censoring
          # checked == True only if X_i is close enough to one of the grid points
          checked <- FALSE
          J <- length(grid_points)
          if (J > 0) {
            for (j in 1:J) {
              # If X_i is close enough to one of grid points, use saved result
              if (abs(X[i] - grid_points[j]) < epsilon) {
                idxs <- sapply(c(Y[i], tau), function(t) sum(t > res_times[[j]]))
                tmp <- sapply(
                  idxs,
                  function(idx) {
                    ifelse(idx > 0, res_survs[[j]][idx], 1)
                  }
                )
                checked <- TRUE
                break
              }
            }
          }
          # If far apart from grid points, add a new grid point
          if (checked == FALSE) {
            grid_points <- c(grid_points, X[i])
            fm <- survfit(Surv(Y, delta) ~ 1,
              weights = calc.kw(X = X, x0 = X[i], span = span, h = h, type = type),
              se.fit = F
            )
            fm_times <- fm$time[which(fm$n.event == 1)]
            fm_survs <- fm$surv[which(fm$n.event == 1)]

            res_times[[J + 1]] <- fm_times
            res_survs[[J + 1]] <- fm_survs

            idxs <- sapply(c(Y[i], tau), function(t) sum(t > fm_times))
            tmp <- sapply(
              idxs,
              function(idx) {
                ifelse(idx > 0, fm_survs[idx], 1)
              }
            )
          }

          if (tmp[1] == 0) {
            positive[i] <- 1
          } else {
            positive[i] <- 1 - tmp[2] / tmp[1]
          }
        }
      }
    }
    # Remove saved K-M models
    rm(res_times, res_survs)
  } else {
    # By default (null epsilon), do the calculation for each right censored point
    for (i in 1:n) {
      if (Y[i] > tau) {
        positive[i] <- 0
      } else {
        if (delta[i] == 1) {
          positive[i] <- 1
        } else {
          kw <- calc.kw(X = X, x0 = X[i], span = span, h = h, type = type)
          fm <- survfit(Surv(Y, delta) ~ 1, weights = kw)
          tmp <- summary(fm, times = c(Y[i], tau))$surv
          if (tmp[1] == 0) {
            positive[i] <- 1
          } else {
            positive[i] <- 1 - tmp[2] / tmp[1]
          }
        }
      }
    }
  }

  negative <- 1 - positive

  if (is.null(X.min)) {
    X.min <- min(X)
  }
  if (is.null(X.max)) {
    X.max <- max(X)
  }

  if (is.null(cut.off)) {
    cut.off <- c(-Inf, seq(X.min, X.max, length = n.grid), Inf)
  }

  sens <- spec <- NULL
  for (this.c in cut.off) {
    sens <- c(sens, sum(positive * as.numeric(X > this.c)) / sum(positive))
    # sensitivity that incorporates fractional "positive"
    spec <- c(spec, sum(negative * as.numeric(X <= this.c)) / sum(negative))
    # specificity that incorporates fractional "negative"
  }
  ROC <- data.frame(
    cut.off = cut.off,
    sens = sens,
    spec = spec
  )
  W <- positive
  main_res <- calibration_res <-AUC <- AUC2 <-  NA

  if (output != "calibration") {
    if(method != "integral"){
      AUC2 <- AUC_calc(X, W)
    }
    if(method != "empirical"){
      AUC <- AUC_calc_integral(sens, spec)
    }
    main_res <- list(AUC.integral = AUC, AUC.empirical = AUC2, ROC = ROC)
  }
  if(output != "AUC"){
    calibration_res <- c(
      BrierScore = tdBrier(W, X),
      SurvDiff = tdSurvDiff(W, X)
    )
  }


  # if bootstrap number is specified
  if (nboot > 0) {
    # start bootstrap for AUC
    boot.AUC.integral <- boot.AUC.empirical <- rep(NA, nboot)
    boot.BS <- boot.SurvDiff <- rep(NA, nboot)

    n.grid <- min(n.grid, length(X))
    # cut.off <- NULL
    if (is.null(cut.off)) {
      cut.off <- seq(min(X), max(X), length.out = n.grid)
    }

    boot.sens <- matrix(NA, nrow = nboot, ncol = length(cut.off))
    boot.spec <- matrix(NA, nrow = nboot, ncol = length(cut.off))

    set.seed(123)
    # the random number seed is hardcoded
    for (b in 1:nboot) {
      loc <- sample(x = 1:n, size = n, replace = T)
      X2 <- X[loc]
      Y2 <- Y[loc]
      delta2 <- delta[loc]
      out <- tdROC(X2, Y2, delta2, tau, span,
        nboot = 0, alpha, n.grid,
        cut.off = cut.off, X.min = X.min, X.max = X.max, epsilon = epsilon, method
      )
      boot.AUC.integral[b] <- out$main_res$AUC.integral
      boot.AUC.empirical[b] <- out$main_res$AUC.empirical
      boot.BS[b] <- out$calibration_res[1]
      boot.SurvDiff[b] <- out$calibration_res[2]

      boot.sens[b, ] <- out$main_res$ROC$sens
      boot.spec[b, ] <- out$main_res$ROC$spec
    }

    boot.AUC.integral <- unlist(boot.AUC.integral)
    boot.AUC.empirical <- unlist(boot.AUC.empirical)

    boot.BS <- unlist(boot.BS)
    boot.SurvDiff <- unlist(boot.SurvDiff)


    bAUC.integral <- bAUC.empirical <-bBS <-bSurvDiff <-bROC <-   NA
    if(output != "calibration"){

      if (method != "integral") {
        bAUC.empirical <- list(
          mean = mean(boot.AUC.empirical),
          sd = sd(boot.AUC.empirical),
          CIlow = quantile(boot.AUC.empirical, prob = alpha / 2),
          CIhigh = quantile(boot.AUC.empirical, prob = 1 - alpha / 2)
        )
      }
      if(method != "empirical"){
        bAUC.integral <- list(
          mean = mean(boot.AUC.integral),
          sd = sd(boot.AUC.integral),
          CIlow = quantile(boot.AUC.integral, prob = alpha / 2),
          CIhigh = quantile(boot.AUC.integral, prob = 1 - alpha / 2)
        )
      }
      sens.mean <- apply(boot.sens, 2, mean)
      sens.sd <- apply(boot.sens, 2, sd)
      sens.lower <- apply(boot.sens, 2, quantile, prob = alpha / 2)
      sens.upper <- apply(boot.sens, 2, quantile, prob = 1 - alpha / 2)

      spec.mean <- apply(boot.spec, 2, mean)
      spec.sd <- apply(boot.spec, 2, sd)
      spec.lower <- apply(boot.spec, 2, quantile, prob = alpha / 2)
      spec.upper <- apply(boot.spec, 2, quantile, prob = 1 - alpha / 2)
      bROC <- data.frame(
        sens.mean, sens.sd, sens.lower, sens.upper,
        spec.mean, spec.sd, spec.lower, spec.upper
      )
    }
    if (output != "AUC") {
      bBS <- list(
        mean = mean(boot.BS),
        sd = sd(boot.BS),
        CIlow = quantile(boot.BS, prob = alpha / 2),
        CIhigh = quantile(boot.BS, prob = 1 - alpha / 2)
      )
      bSurvDiff <- list(
        mean = mean(boot.SurvDiff),
        sd = sd(boot.SurvDiff),
        CIlow = quantile(boot.SurvDiff, prob = alpha / 2),
        CIhigh = quantile(boot.SurvDiff, prob = 1 - alpha / 2)
      )
    }


    boot_res <- list(
      bAUC.integral = bAUC.integral, bAUC.empirical = bAUC.empirical,
      bBS = bBS, bSurvDiff = bSurvDiff, bROC = bROC
    )
  } else {
    boot_res <- "Please set `nboot` to obtain bootstrap result."
  }
  message("The calibration results (Brier score, survival difference) are applicable only when the risk score has the predicted probability interpretation.")
  output <- list(main_res = main_res, calibration_res = calibration_res, boot_res = boot_res)
  return(output)
}
