#' PKPDindex_data
#'
#' Example dataset for Emax model fitting.
#'
#' This dataset contains information about drug concentrations (AUC/MIC, Cmax/MIC, and Time above MIC)
#' and their corresponding response values, used for modelling the drug's effect based on the Emax model.
#'
#' @name PKPDindex_data
#' @docType data
#' @format A data frame with 20 rows and 4 columns:
#' \describe{
#'   \item{auc_mic}{Area under the concentration-time curve (numeric)}
#'   \item{cmax_mic}{Maximum concentration of the drug (numeric)}
#'   \item{t_mic}{Time above minimum inhibitory concentration (numeric)}
#'   \item{response}{Observed drug response (numeric; Delta log10 CFU/ml)}
#' }
#' @source Generated for package example purposes.
#' @keywords datasets
#' @examples
#' data(PKPDindex_data)
"PKPDindex_data"
