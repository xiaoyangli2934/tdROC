#' Example data: Malignant Melanoma Data
#'
#' This example dataset is included for illustration. The Melano data is publicly available.
#'    In 1962-1977, 205 patients with malignant melanoma (skin cancer) had a radical
#'    operation performed at an academic medical center. At the end of the follow-up,
#'    57 died from cancer, 14 died from other causes, and the other 134 patients were
#'    alive (censored). This example dataset illustrates the prediction accuracy of
#'    competing risk outcomes with baseline age and tumor thickness.

#'
#' @docType data
#' @usage data(Melano)
#' @keywords datasets
#' @format A data frame with 312 observations and 4 variables: time (event time/censoring), time), censor (censoring
#'         mayoscore4, mayoscore5. The two scores are derived from 4 and 5 covariates
#'        respectively.
#' @references Andersen, P. K. , & Skovgaard, L. T. (2010). Regression with linear predictors. New York, NY: Springer.
#'
#' @examples
#' data(Melano)
#' head(Melano)
"Melano"
