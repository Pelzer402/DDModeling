#' Function to fit  a given DDRep to a given DDModel
#'
#' @name Fit_DDModel
#'
#' @param model \code{DDModel} object
#' @param rep   \code{DDRep} object
#' @param grid  optional \code{path} to a direcotey containing a .GRID fileset. If NULL the model will be fitted using 20 randomly drawn startparametersets from the model-DOMAIN.
#'
#' @return \code{DDFit} object
Fit_DDModel <- function(model = NULL, rep = NULL, grid = NULL){
  Flag <- 1
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model))
  {
    ArgumentCheck::addError(msg = "'model' is missing!",argcheck = Check)
    Flag <- 99
  }
  if (!isS4(model))
  {
    ArgumentCheck::addError(msg = "'model' must be a DDModel object!",argcheck = Check)
    Flag <- 99
  }
  if (is.null(rep))
  {
    ArgumentCheck::addError(msg = "'rep' is missing!",argcheck = Check)
    Flag <- 99
  }
  if(!isS4(rep))
  {
    ArgumentCheck::addError(msg = "'rep' must be a DDRep object!",argcheck = Check)
    Flag <- 99
  }
  ArgumentCheck::finishArgCheck(Check)
  if (Flag == 99)
  {
    return(cat("Fit_DDModel failed"))
  }
  else
  {
    if (is.null(grid))
    {
      model@FORM <- rep@FORM
      return(.Fit_DDModel(model,rep))
    }
    else
    {
      model@FORM <- rep@FORM
      return(.Fit_DDModel(model,rep))
    }

  }
}

