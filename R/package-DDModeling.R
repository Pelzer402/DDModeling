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
#' @import Rcpp
NULL

