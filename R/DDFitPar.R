#' \code{DDFitPar} class definition
#' @name DDFitPar-class
#' @rdname DDFitPar-class
#' @slot FIT_V      \code{numeric} representing the value of the Fit
#' @slot FIT_N      \code{numeric} representing the number of evaluation points used in the fit
#' @slot START_VALUE \code{character} specifying the method used for the start value problem
#' @slot METHOD     \code{character} specifying the method used for fitting
#' @slot S_SAMPLING \code{logical} if super sampling was used
#' @slot TRIALS_TOAL \code{numeric} specifying the total numbers of trials that were simulated
#' @slot TRIALS_SAMPLE \code{numeric} specifying the number of trials inside a given sample
setClass("DDFitPar",
         slots      = list(
           FIT_V         = "numeric",
           FIT_N         = "numeric",
           START_VALUE   = "character",
           METHOD        = "character",
           S_SAMPLING    = "logical",
           TRIALS_TOTAL  = "numeric",
           TRIALS_SAMPLE = "numeric"
         )
)

