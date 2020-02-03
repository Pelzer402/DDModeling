#' \code{DDFit} class definition
#' @include DDRep.R DDModel.R DDFitPar.R
#' @name DDFit-class
#' @rdname DDFit-class
#' @slot INP_REP    \code{DDRep} object containing the to be fitted data
#' @slot FIT_REP    \code{DDRep} object containing the fitted data
#' @slot MODEL      \code{DDModel} object containing the model that was used in the fit
#' @slot FIT        \code{DDFitPar} object containing information regarding the fit
setClass("DDFit",
         slots      = list(
           INP_REP      = "DDRep",
           FIT_REP      = "DDRep",
           MODEL        = "DDModel",
           FIT          = "DDFitPar"
         )
)
