% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fit.R
\name{Fit_DDModel}
\alias{Fit_DDModel}
\title{Function to fit a given DDRep to a given DDModel}
\usage{
Fit_DDModel(
  model = NULL,
  data = NULL,
  DL_model = NULL,
  DL_scale = list(INPUT = list(Center = FALSE, Stddev = FALSE), OUTPUT = list(Center =
    FALSE, Stddev = FALSE)),
  grid_path = NULL,
  trials = 10000L,
  s_sampling = FALSE,
  simplex_struc = c(0L),
  simplex_coef = list(alpha = 1, beta = 0.5, gamma = 2, sigma = 0.5, tol = 1e-07, nfunc
    = 120L, nshrink = 3L)
)
}
\arguments{
\item{model}{\code{DDModel} object}

\item{data}{\code{DDRep} object or list of \code{DDRep} objects}

\item{DL_model}{(optional) \code{Model} in the form of a keras neural network model}

\item{DL_scale}{(optional) \code{list} containing mean and sd of the transformation used in the deep learning model while training.}

\item{grid_path}{(optional) \code{path} to a directory containing a .GRID fileset. If NULL the model will be fitted using 20 randomly drawn startparametersets from the model-DOMAIN.}

\item{trials}{\code{integer} indicating the number of trials used while fitting (s_sampling = FALSE) or the maximum number of trials used while super sampling (s_sampling = TRUE)}

\item{s_sampling}{\code{bool} indicating super sampling while fitting}

\item{simplex_struc}{\code{numeric vector} containing the number of simplex iterations per sorting cycle.}

\item{simplex_coef}{\code{list} containing named parameters for the SIMPLEX}
}
\value{
\code{DDFit} object
}
\description{
Function to fit a given DDRep to a given DDModel
}
\details{
This is a rather complex function that is very flexible depending on your input!
\describe{
\item{Random Fitting:}{Exclude all optional parameters}
\item{Grid Fitting:}{Specify a grid_path}
\item{Deep Learning Fitting}{Specify DL_model and DL_scale}
}
}
