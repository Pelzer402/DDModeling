#' \code{DDFitPar} class definition
#' @name DDFitPar-class
#' @rdname DDFitPar-class
#' @slot INP_P      \code{data.frame} containing the to be fitted parameter
#' @slot FIT_P      \code{data.frame} containing the fitted parameter
#' @slot FIT_V      \code{numeric} representing the value of the Fit
#' @slot FIT_N      \code{numeric} representing the number of evaluation points used in the fit
setClass("DDFitPar",
         slots      = list(
           INP_P      = "data.frame",
           FIT_P      = "data.frame",
           FIT_V      = "numeric",
           FIT_N      = "numeric"
         )
)
