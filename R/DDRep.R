#' An S4 class to represent a representation of a drift diffusion simulation
#'
#' @name DDRep-class
#' @rdname DDRep-class
#' @slot RAW  \code{list} of data.frames that contain the RAW data (i.e. 3 coloumns: $cond $Resp $time)
#' @slot REP  \code{list} of data.frames that contain data representations (CDF and CAF)
#' @slot RF \code{list} of numeric vectors that contain the percentiles of the representation.
setClass("DDRep",
         slots      = list(
           RAW      = "list",
           REP      = "list",
           RF     = "list"
         )
)

DDRep <- function(model = NULL,raw=NULL){
  return(.DDRep_cpp(model,raw))
}
