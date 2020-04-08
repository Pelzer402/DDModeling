#' Function to fit a given DDRep to a given DDModel
#' @name Fit_DDModel
#' @param model \code{DDModel} object
#' @param data   \code{DDRep} object or list of \code{DDRep} objects
#' @param grid_path  \code{path} to a directory containing a .GRID fileset. If NULL the model will be fitted using 20 randomly drawn startparametersets from the model-DOMAIN.
#' @param s_sampling \code{bool} indicating super sampling while fitting
#' @param DL_model \code{Model} in the form of a keras neural network model
#' @param trials \code{integer} indicating the number of trials used while fitting (s_sampling = FALSE) or the maximum number of trials used while super sampling (s_sampling = TRUE)
#' @param simplex_struc \code{numeric vector} containing the number of simplex iterations per sorting cycle.
#' @return \code{DDFit} object
#' @export
Fit_DDModel <- function(model = NULL, data = NULL, DL_model = NULL, grid_path = NULL, s_sampling = FALSE, trials = 10000, simplex_struc= c(20,2)){
  Method <- 0 # 1 == random, 2 == grid, 3 == DL, 99 == Error
  Output <- 0 # 1 == single, 2 == list
  Check <- ArgumentCheck::newArgCheck()
  if (!methods::is(model,"DDModel"))
  {
    ArgumentCheck::addError(msg = "'model' is missing or in the wrong format!",argcheck = Check)
    Method <- 99
  }
  if (!methods::is(data,"DDRep"))
  {
    if (is.list(data))
    {
      if (all(unlist(lapply(data,function(x){methods::is(x,"DDRep")}))))
      {
        Output <- 2
      }
      else
      {
        ArgumentCheck::addError(msg = "'data' is missing or in the wrong format!",argcheck = Check)
        Method <- 99
      }
    }
    else
    {
      ArgumentCheck::addError(msg = "'data' is missing or in the wrong format!",argcheck = Check)
      Method <- 99
    }
  }
  else
  {
    Output <- 1
  }
  if (is.null(grid_path) & is.null(DL_model))
  {
    ArgumentCheck::addWarning(msg = "'grid_path' or 'DL_model' is missing!\n Startvalues will determined using a uniform random distribution!",argcheck = Check)
    Method <- 1
  }
  else if (!is.null(grid_path) & !is.null(DL_model))
  {
    ArgumentCheck::addError(msg = "'grid_path' and 'DL_model' is specified!\n You may only use one!",argcheck = Check)
    Method <- 99
  }
  else if (!is.null(grid_path))
  {
    Grid_model <- readRDS(list.files(grid_path,full.names = TRUE,pattern = "\\.Gcfg$"))
    #if (!identical(Grid_model,model))
    #{
    #  ArgumentCheck::addError(msg = "The grid under 'grid_path' does not comform to the model under 'model'!",argcheck = Check)
    #  Method <- 99
    #}
    #else
    #{
    #  Method <- 2
    #}
    Method <- 2
  }
  else if (!is.null(DL_model))
  {
    Method <- 3
  }
  if (!is.logical(s_sampling))
  {
    ArgumentCheck::addError(msg = "'s_sampling' needs to be of format boolean!",argcheck = Check)
    Method <- 99
  }
  if (!trials%%1==0)
  {
    ArgumentCheck::addError(msg = "'trials' needs to be of format integer!",argcheck = Check)
    Method <- 99
  }
  if (!is.numeric(simplex_struc) || length(simplex_struc)==0)
  {
    ArgumentCheck::addError(msg = "'simplex_struc' is empty or in the wrong format!",argcheck = Check)
    Method <- 99
  }
  ArgumentCheck::finishArgCheck(Check)
  if (Method == 99)
  {
    cat("Fit_DDModel failed")
  }
  else if (Method == 1)
  {
    if (Output == 2)
    {
      ncores <- parallel::detectCores() - 1
      clust <- ParallelLogger::makeCluster(ncores)
      COMP_List <- list()
      for (i in 1:length(data))
      {
        COMP_List[[i]] <-list(model,data[[i]],s_sampling,trials,simplex_struc)
      }
      RES <- ParallelLogger::clusterApply(clust,COMP_List,.Fit_DDModel_rnd)
      ParallelLogger::stopCluster(clust)
      return(RES)
    }
    else if(Output == 1)
    {
      COMP_List <- list(model,data,s_sampling,trials,simplex_struc)
      return(.Fit_DDModel_rnd(COMP_List))
    }
  }
  else if (Method == 2)
  {
    Grid_parts <- list.files(grid_path,full.names = TRUE,pattern = "\\.GRID$")
    if (Output == 2)
    {
      ncores <- parallel::detectCores() - 1
      clust <- ParallelLogger::makeCluster(ncores)
      COMP_List <- list()
      for (i in 1:length(data))
      {
        COMP_List[[i]] <-list(model,data[[i]],Grid_parts,s_sampling,trials,grid_path,simplex_struc)
      }
      RES <- ParallelLogger::clusterApply(clust,COMP_List,.Fit_DDModel_grid)
      ParallelLogger::stopCluster(clust)
      return(RES)
    }
    else if(Output == 1)
    {
      COMP_List <- list(model,data,Grid_parts,s_sampling,trials,grid_path,simplex_struc)
      return(.Fit_DDModel_grid(COMP_List))
    }
  }
  else if (Method == 3)
  {
    if (Output == 2)
    {
      ncores <- parallel::detectCores() - 1
      clust <- ParallelLogger::makeCluster(ncores)
      COMP_List <- list()
      for (i in 1:length(data))
      {
        INP <- t(DDRep_wide(data[[i]]))
        PRE <- as.numeric(predict(DL_model,INP))
        COMP_List[[i]] <-list(model,data[[i]],s_sampling,trials,simplex_struc,PRE)
      }
      RES <- ParallelLogger::clusterApply(clust,COMP_List,.Fit_DDModel_DL)
      ParallelLogger::stopCluster(clust)
      return(RES)
    }
    else if(Output == 1)
    {
      INP <- t(DDRep_wide(data))
      PRE <- as.numeric(predict(DL_model,INP))
      COMP_List <- list(model,data,s_sampling,trials,simplex_struc,PRE)
      return(.Fit_DDModel_DL(COMP_List))
    }
  }
}



