#' An S4 class to represent a Fit of a given DDRep to a model
#'
#' @name DDFit-class
#' @rdname DDFit-class
#' @slot Input_Rep  \code{DDRep} object containing the to be fitted data
#' @slot Fit_Rep    \code{DDRep} object containing the fitted data
#' @slot Model      \code{DDModel} object containing the model that was used in the fit
#' @slot Fit        \code{DDFitPar} object containing information regarding the fit
setClass("DDFit",
         slots      = list(
           Input_Rep    = "DDRep",
           Fit_Rep      = "DDRep",
           Model        = "DDModel",
           Fit          = "DDFitPar"
         )
)
