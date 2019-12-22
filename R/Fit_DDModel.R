#' Function to fit  a given DDRep to a given DDModel
#'
#' @name Fit_DDModel
#'
#' @param model \code{DDModel} object
#' @param rep   \code{DDRep} object
#' @param grid_path  \code{path} to a direcotey containing a .GRID fileset. If NULL the model will be fitted using 20 randomly drawn startparametersets from the model-DOMAIN.
#'
#' @return \code{DDFit} object
Fit_DDModel <- function(model = NULL, rep = NULL, grid_path = NULL){
  Flag <- 1
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model) || !is(model,"DDModel"))
  {
    ArgumentCheck::addError(msg = "'model' is missing or in the wrong format!",argcheck = Check)
    Flag <- 99
  }
  if (is.null(rep) || !is(rep,"DDRep"))
  {
    ArgumentCheck::addError(msg = "'rep' is missing or in the wrong format!",argcheck = Check)
    Flag <- 99
  }
  if (is.null(grid_path))
  {
    ArgumentCheck::addError(msg = "'grid_path' is missing!",argcheck = Check)
    Flag <- 99
  }
  else
  {
    Grid_model <- readRDS(list.files(grid_path,full.names = TRUE,pattern = ".Gcfg"))
    if (!identical(Grid_model,model))
    {
      ArgumentCheck::addError(msg = "The grid under 'grid_path' does not comform to the model under 'model'!",argcheck = Check)
      Flag <- 99
    }
  }
  ArgumentCheck::finishArgCheck(Check)
  if (Flag == 99)
  {
    return(cat("Fit_DDModel failed"))
  }
  else
  {
      Grid_parts <- list.files(grid_path,full.names = TRUE,pattern = ".GRID")
      return(.Fit_DDModel_grid(model,rep,Grid_parts))
  }
}

