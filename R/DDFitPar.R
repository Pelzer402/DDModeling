#' \code{DDFitPar} class definition
#' @name DDFitPar-class
#' @rdname DDFitPar-class
#' @slot FIT_V      \code{numeric} representing the value of the Fit
#' @slot FIT_N      \code{numeric} representing the number of evaluation points used in the fit
setClass("DDFitPar",
         slots      = list(
           FIT_V      = "numeric",
           FIT_N      = "numeric"
         )
)
