# Clean up .dll
.onUnload <- function (libpath) {library.dynam.unload("DDModeling", libpath)}

#' DDModeling pakage
#'
#' A pakage for the easy integration of dirft diffusion models in cognitive psychology
#'
#' @docType package
#' @name DDModeling
#' @title DDModeling Pakage
#' @author Thomas Pelzer
#' @useDynLib DDModeling
NULL

