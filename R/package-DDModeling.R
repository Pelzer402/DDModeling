# Clean up .dll
.onUnload <- function (libpath) {library.dynam.unload("DDModeling", libpath)}
#' DDModeling pakage
#'
#' A pakage for the easy integration of dirft diffusion models in cognitive psychology
#'
#' @docType package
#' @name DDModeling
#' @title DDModeling
#' @author Thomas Pelzer
#' @useDynLib DDModeling
#' @import Rcpp data.table
#' @importFrom stats predict time
#' @importFrom utils read.table
NULL


#' Custom implementations of the show method for DDModeling
#'
#' @name show-methods
#' @param object object to be shown (see usage for data types)
NULL

#' Custom implementations of the plot method for DDModeling
#'
#' @name plot-methods
#' @param x to be plotted (see usage for data types)
NULL

#' Custom implementations of the summary method for DDModeling
#'
#' @name summary-methods
#' @param object object to be summarized (see usage for data types)
NULL

#' Custom implementations of the Compare method for DDModeling
#'
#' @name Compare-methods
#' @param e1 a DDRep object
#' @param e2 either another DDRep object or a character specifying a grid_path
NULL
