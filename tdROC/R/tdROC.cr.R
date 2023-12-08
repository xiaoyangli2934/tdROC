#' Estimate time-dependent prediction accuracy measures, including the ROC, AUC, Brier score, and survival probability difference, with competing risk data.
#'
#' This is a core function of the ‘tdROC‘ package. It uses the nonparametric weights proposed by Wu (Wu and Li, 2018)
#' to estimate a number of time-dependent prediction accuracy measures for right-censored survival outcomes,
#' including ROC curve, AUC, Brier score, and survival difference, with competing risk data. For each measure, the variance can be estimated
#' through bootstrap resampling.
#'
#' @param X  a numeric vector of risk score in the same length as \code{Y} and \code{delta}, such as biomarker or predicted probability. A higher value indicates higher risk of the event.
#'  The calibration results (Brier score, survival difference) are applicable only when the risk score has the predicted probability interpretation.
#' @param Y a numeric vector of time to event in the same length as \code{X} and \code{delta}.
#' @param delta a vector of numeric indicator of event type in the same length as \code{X} and \code{Y}. The primary event should be coded as 1,
#'        the competing event should be coded as 2, and censoring should be coded as 0.
#' @param tau a scalar, the prediction horizon at which the prediction is evaluated.
#' @param span a numeric value, the proportion of neighbour observations used in nearest neighbor method. The default is 0.1.
#' @param h a numeric value, the bandwidth of kernel weights, the defualt is \code{NULL}. If not specified, the function will use the value of
#'  \code{span} to calculate kernel weights. In case both \code{span} and \code{h} are specified, the function will use \code{h}.
#' @param type a character value, indicating the type of kernel function used to calculate kernel weights. The default is \code{"uniform"} kernel. Other options are \code{"Epanechnikov"} and \code{"normal"}.
#'  It will only be used when the bandwidth \code{h} is specified.
#' @param n.grid a positive integer, the number of grid points used when calculating the ROC curve. The default is \code{1000}.
#' @param cut.off a vector of \code{X} cut-off values at which sensitivity and specificity will be calculated.
#' @param nboot the number of bootstrap replications to be used for variance estimation. The default is \code{nboot = 0}, corresponding to no variance estimation.
#' @param alpha It is (1 - level of confidence interval)/2, default is \code{0.05}. It is used only when \code{nboot > 0}.
#' @param epsilon The precision parameter used in an approximation to the weight calculation when the sample size is large. If a weight corresponding to a specific risk score is already calculated,
#'   then the weights corresponding to adjacent risk scores, within the distance specified by epsilon, will be the same under the approximation. This approximation avoids repeated
#'   calculation of weights that are almost the same, and hence increases the speed of computation in this situation. The default is NULL, which means no approximation is used. A higher
#'   value indicates less precision.
#' @param method It is used to specify which method you would like to use to estimate AUC, default to \code{"both"}. Other options are \code{"integral"} and \code{"empirical"}.
#' @param output It is used to specify which kind of output you want, default to \code{"both"}. Other options are \code{"AUC"}, including AUC, sensitivity, and specificity are included,
#'  and \code{"calibration"} including Brier Score and survival difference.
#'
#'
#' @details This function takes the risk score value \code{X}, the time-to-event data \code{Y} and censoring indicator \code{delta} as input to estimate
#'  a number of time-dependent prediction accuracy measures for survival outcomes, including ROC curve, AUC, Brier score, and survival difference, with competing risk.
#'  The confidence intervals of above quantities are estimated by bootstrap.
#'
#'  For competing risk data, there are two definition of controls introduced by Zheng et al, which are listed below
#'
#'  \deqn{
#'  \text{Definition A:} \text{Case} k:T \le \tau, \delta = k; \text{Control}_A: (T>\tau)\cup (T \le \tau \cap \delta \ne k)
#'  }
#'
#'  \deqn{
#'  \text{Definition B:} \text{Case} k:T \le \tau, \delta = k; \text{Control}_B: (T>\tau)
#'  }
#'
#'  Based on the definition A, both the event-free subjects and subjects who experience other competing events were included as controls. While definition B include only event-free subjects.
#'  This function offers two options to estimate AUC. The first one make use of estimated sensitivity and specificity to calculate the AUC via trapezoidal integration
#'  by setting a series of cutoff point. For the two different definition, we separately calculate the sensitivity, specificity and AUC. The output will also include the sensitivities
#'  and specificities for our plot function. The other one estimates AUC by the empirical estimator of the proportion of concordance pairs with proposed weight estimator (Wu and Li, 2018).
#'  These two methods generate quite similar estimates. The option can be set by the argument \code{method}.
#'
#'  In addition to the above prediction measures, we include Brier Score and survival difference to evaluate the calibration metrics. Their definitions are included below.
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
#'  This function uses the same approximation as \code{tdROC} with the argument \code{epsilon}
#'
#' @return Returns a list of the following items:
#
#' @return \code{main_res:} a list of \code{AUC.A.integral} estimated by trapezoidal integration for definition A,
#'        \code{AUC.A.empirical} estimated by empirical estimator for definition A,
#'         \code{AUC.B.integral} estimated by trapezoidal integration for definition B,
#'        \code{AUC.B.empirical} estimated by empirical estimator for definition B,
#'        and a data frame \code{ROC} with dimension \code{(2+n.grid) x 4} with columns \code{cut.off}, \code{sens}, \code{specA} and \code{specB}.
#' @return \code{calibration_res:} brier score and survival difference estimated based on the formula similar to Wu and Li (2018). When the risk score \code{X}
#' is a biomarker value instead of a predicted cumulative incidence probability, the brier score and survival difference cannot be calculated. In this case, please disregard the calibration results.
#' @return \code{boot_res:} a list of bootstrap results, including \code{bAUC.A.integral}, \code{bAUC.A.empirical},  \code{bAUC.B.integral}, \code{bAUC.B.empirical}, \code{bBS}, \code{bSurvDiff}, \code{bROC}.
#'  For \code{bAUC.A.integral}, \code{bAUC.A.empirical},  \code{bAUC.B.integral}, \code{bAUC.B.empirical}, \code{bBS}, \code{bSurvDiff}, each one is a list including corresponding mean, standard deviation, confidence interval.
#'  \code{bROC} is a data frame with colomns \code{sens.mean}, \code{sens.sd}, \code{sens.lower}, \code{sens.upper},
#'  \code{specA.mean}, \code{specA.sd}, \code{specA.lower}, \code{specA.upper}, \code{specB.mean}, \code{specB.sd}, \code{specB.lower}, \code{specB.upper}

#' @importFrom survival survfit
#' @importFrom stats sd quantile
#'
#' @references Zheng Y, Cai T, Jin Y, Feng Z. Evaluating prognostic accuracy of biomarkers under competing risk. Biometrics. 2012;68(2):388-396. doi:10.1111/j.1541-0420.2011.01671.x
#'
#' @examples
#' library(survival)
#' data(Melano)
#' expit <- function(x){ 1/(1+exp(-x)) }
#' tdROC.cr(X = expit(Melano$thick) , Y = Melano$time, delta = Melano$status, tau = 1800, nboot = 100)
#'
#' @export


tdROC.cr <- function(X, Y, delta, tau,
                     span = 0.1, h = NULL, type = "uniform",
                     epsilon = 0.01, cut.off = NULL, n.grid = 1000, nboot = 0, alpha = 0.05,
                     method = "both",
                     output = "both") {
  n <- length(X)
  W1 <- td.kw.cr(X, Y, delta, event.code = 1, tau, span, h, type, epsilon)
  W2 <- td.kw.cr(X, Y, delta, event.code = 2, tau, span, h, type, epsilon)
  main_res <- calibration_res <- NA
  if (output != "calibration") {
    main_res <- AUC.cr(X, W.prim = W1, W.cmp = W2, cut.off = cut.off, n.grid = n.grid, method)
  }
  if(output != "AUC"){
    calibration_res <- c(
      BrierScore = tdBrier(W1, X),
      SurvDiff = tdSurvDiff(W1, X)
    )
  }

  if (nboot > 0) {
    boot.AUC.A.integral <- boot.AUC.A.empirical <- boot.AUC.B.integral <- boot.AUC.B.empirical <- rep(NA, nboot)
    boot.BS <- boot.SurvDiff <- rep(NA, nboot)

    n.grid <- min(n.grid, length(X))
    # cut.off <- NULL
    if (is.null(cut.off)) {
      cut.off <- seq(min(X), max(X), length.out = n.grid)
    }

    boot.sens <- matrix(NA, nrow = nboot, ncol = length(cut.off))
    boot.specA <- matrix(NA, nrow = nboot, ncol = length(cut.off))
    boot.specB <- matrix(NA, nrow = nboot, ncol = length(cut.off))

    set.seed(123)
    for (i in 1:nboot) {
      loc <- sample(x = 1:n, size = n, replace = T)
      Xb <- X[loc]
      Yb <- Y[loc]
      deltab <- delta[loc]
      W1b <- td.kw.cr(Xb, Yb, deltab, event.code = 1, tau, span, h, type, epsilon)
      W2b <- td.kw.cr(Xb, Yb, deltab, event.code = 2, tau, span, h, type, epsilon)

      if (output != "calibration") {
        tdROC_res <- AUC.cr(X = Xb, W.prim = W1b, W.cmp = W2b, cut.off, n.grid, method)

        boot.AUC.A.integral[i] <- tdROC_res$AUC.A.integral
        boot.AUC.A.empirical[i] <- tdROC_res$AUC.A.empirical
        boot.AUC.B.integral[i] <- tdROC_res$AUC.B.integral
        boot.AUC.B.empirical[i] <- tdROC_res$AUC.B.empirical
        boot.sens[i, ] <- tdROC_res$ROC$sens
        boot.specA[i, ] <- tdROC_res$ROC$specA
        boot.specB[i, ] <- tdROC_res$ROC$specB
      }
      if(output != "AUC"){
        boot.BS[i] <- tdBrier(W1b, Xb)
        boot.SurvDiff[i] <- tdSurvDiff(W1b, Xb)
      }


    }

    boot.AUC.A.integral <- unlist(boot.AUC.A.integral)
    boot.AUC.A.empirical <- unlist(boot.AUC.A.empirical)
    boot.AUC.B.integral <- unlist(boot.AUC.B.integral)
    boot.AUC.B.empirical <- unlist(boot.AUC.B.empirical)
    boot.BS <- unlist(boot.BS)
    boot.SurvDiff <- unlist(boot.SurvDiff)


    bAUC.A.integral <- bAUC.A.empirical <- bAUC.B.integral <- bAUC.B.empirical <-bBS <-bSurvDiff <-bROC <-   NA
    if(output != "calibration"){

      if (method != "integral") {
        bAUC.A.empirical <- list(
          mean = mean(boot.AUC.A.empirical),
          sd = sd(boot.AUC.A.empirical),
          CIlow = quantile(boot.AUC.A.empirical, prob = alpha / 2),
          CIhigh = quantile(boot.AUC.A.empirical, prob = 1 - alpha / 2)
        )
        bAUC.B.empirical <- list(
          mean = mean(boot.AUC.B.empirical),
          sd = sd(boot.AUC.B.empirical),
          CIlow = quantile(boot.AUC.B.empirical, prob = alpha / 2),
          CIhigh = quantile(boot.AUC.B.empirical, prob = 1 - alpha / 2)
        )
      }
      if(method != "empirical"){
        bAUC.A.integral <- list(
          mean = mean(boot.AUC.A.integral),
          sd = sd(boot.AUC.A.integral),
          CIlow = quantile(boot.AUC.A.integral, prob = alpha / 2),
          CIhigh = quantile(boot.AUC.A.integral, prob = 1 - alpha / 2)
        )
        bAUC.B.integral <- list(
          mean = mean(boot.AUC.B.integral),
          sd = sd(boot.AUC.B.integral),
          CIlow = quantile(boot.AUC.B.integral, prob = alpha / 2),
          CIhigh = quantile(boot.AUC.B.integral, prob = 1 - alpha / 2)
        )
      }

      sens.mean <- apply(boot.sens, 2, mean)
      sens.sd <- apply(boot.sens, 2, sd)
      sens.lower <- apply(boot.sens, 2, quantile, prob = alpha / 2)
      sens.upper <- apply(boot.sens, 2, quantile, prob = 1 - alpha / 2)

      specA.mean <- apply(boot.specA, 2, mean)
      specA.sd <- apply(boot.specA, 2, sd)
      specA.lower <- apply(boot.specA, 2, quantile, prob = alpha / 2)
      specA.upper <- apply(boot.specA, 2, quantile, prob = 1 - alpha / 2)

      specB.mean <- apply(boot.specB, 2, mean)
      specB.sd <- apply(boot.specB, 2, sd)
      specB.lower <- apply(boot.specB, 2, quantile, prob = alpha / 2)
      specB.upper <- apply(boot.specB, 2, quantile, prob = 1 - alpha / 2)
      bROC <- data.frame(
        sens.mean, sens.sd, sens.lower, sens.upper,
        specA.mean, specA.sd, specA.lower, specA.upper,
        specB.mean, specB.sd, specB.lower, specB.upper
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
      bAUC.A.integral = bAUC.A.integral, bAUC.A.empirical = bAUC.A.empirical, bAUC.B.integral = bAUC.B.integral, bAUC.B.empirical = bAUC.B.empirical,
      bBS = bBS, bSurvDiff = bSurvDiff, bROC = bROC
    )
  } else {
    boot_res <- "Please set `nboot` to obtain bootstrap result."
  }
  message("The calibration results (Brier score, survival difference) are applicable only when the risk score has the predicted probability interpretation.")
  res <- list(main_res = main_res, calibration_res = calibration_res, boot_res = boot_res)
  return(res)
}
