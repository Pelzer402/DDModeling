#' Function to fit  a given DDRep to a given DDModel
#'
#' @name Fit_DDModel
#'
#' @param model \code{DDModel} object
#' @param rep   \code{DDRep} object or list of \code{DDRep} objects
#' @param grid_path  \code{path} to a direcotey containing a .GRID fileset. If NULL the model will be fitted using 20 randomly drawn startparametersets from the model-DOMAIN.
#' @param s_sampling \code{bool} indicating super sampling while fitting
#' @param trials \code{integer} indicating the number of trials used while fitting (s_sampling = FALSE) or the maximum number of trials used while super sampling (s_sampling = TRUE)
#' @return \code{DDFit} object
Fit_DDModel <- function(model = NULL, rep = NULL, grid_path = NULL, s_sampling = FALSE, trials = 10000){
  Flag <- 1
  Check <- ArgumentCheck::newArgCheck()
  if (!is(model,"DDModel"))
  {
    ArgumentCheck::addError(msg = "'model' is missing or in the wrong format!",argcheck = Check)
    Flag <- 99
  }
  if (!is(rep,"DDRep"))
  {
    if (is.list(rep))
    {
      if (all(unlist(lapply(rep,function(x){is(x,"DDRep")}))))
      {
        Flag <- 10
      }
      else
      {
        ArgumentCheck::addError(msg = "'rep' is missing or in the wrong format!",argcheck = Check)
        Flag <- 99
      }
    }
    else
    {
      ArgumentCheck::addError(msg = "'rep' is missing or in the wrong format!",argcheck = Check)
      Flag <- 99
    }
  }
  if (is.null(grid_path))
  {
    ArgumentCheck::addError(msg = "'grid_path' is missing!",argcheck = Check)
    Flag <- 99
  }
  if (!is.logical(s_sampling))
  {
    ArgumentCheck::addError(msg = "'s_sampling' needs to be of format boolean!",argcheck = Check)
    Flag <- 99
  }
  if (!trials%%1==0)
  {
    ArgumentCheck::addError(msg = "'trials' needs to be of format integer!",argcheck = Check)
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
    if (Flag == 10)
    {
      ncores <- parallel::detectCores() - 1
      clust <- parallel::makeCluster(ncores)
      COMP_List <- list()
      for (i in 1:length(rep))
      {
        COMP_List[[i]] <-list(model,rep[[i]],Grid_parts,s_sampling,trials)
      }
      RES <- parallel::clusterApply(clust,COMP_List,.Fit_DDModel_grid)
      parallel::stopCluster(clust)
      return(RES)
    }
    else
    {
      COMP_List <- list(model,rep,Grid_parts,s_sampling,trials)
      return(.Fit_DDModel_grid(COMP_List))
    }


  }
}

