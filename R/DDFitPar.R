#' An S4 class to represent a drift diffusion model
#'
#' @name DDFitPar-class
#' @rdname DDFitPar-class
#' @slot Parameter \code{data.frame} containing the fitted parameter
#' @slot Fit       \code{numeric} representing the value of the Fit
#' @slot nFit      \code{numeric} representing the number of evaluation points used in the fit
setClass("DDFitPar",
         slots      = list(
           Parameter     = "data.frame",
           Fit           = "numeric",
           nFit          = "numeric"
         )
)
