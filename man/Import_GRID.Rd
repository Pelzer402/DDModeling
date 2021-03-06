% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRID.R
\name{Import_GRID}
\alias{Import_GRID}
\title{Imports a GRID into the R environment}
\usage{
Import_GRID(grid_path = NULL, to = "frame")
}
\arguments{
\item{grid_path}{\code{path} to a directory containing a .GRID fileset.}

\item{to}{\code{character} that specifies the format to import to (choose from "frame", "keras_data" and "DDRep")}
}
\value{
\code{list} of \code{DDRep}, \code{list} of \code{matrix} or \code{data.frame} dependend on the "to" parameter
}
\description{
Imports a GRID into the R environment
}
\details{
"frame" will return a \code{data.frame} like structure to the total Grid stored in the grid_path. "keras_data" will return a \code{list}
containing an INPUT and OUTPUT element. The INPUT constitutes all relevant values from the CDF/CAF (namely reaction times and accuracy) in a \code{matrix} with each row
corresponding to one parameter set in the OUTPUT. "DDRep" will convert the complete Grid to DDReps in a \code{list}.
}
