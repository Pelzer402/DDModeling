#' Function to generate a Grid from a given DDModel
#'
#' @name Get_Grid
#' @rdname Get_Grid
#' @slot model \code{DDModel} object
#' @slot path \code{character} that specifies the full path to the directory in which the Grid should be saved
#' @slot name \code{character} that represents the name (and subdirectory in path) of the Grid
#' @return No direct return value inside of the R-session. The calculated Grid will be saved in the specified path!
#' @description  After calling the function the user will be instructed to enter the step sizes corresponding to the parameters listet in the used model.
#'               Step size should allways be of an integer value, as they represent the number of evaluation points per parameter that are used in the grid.
#'               Note that in the given function the evaluation points are allways equally spaced concerning the corresponding parameter domain in the used model.
#'               Therefor, If one would like to specify the used evaluation points it is advised to specify the domain in the model.
Get_Grid <- function(model = NULL, path = NULL, name = NULL){
  Flag <- NULL
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model) || !is(model,"DDModel"))
  {
    ArgumentCheck::addError(msg = "'model' is missing or in the wrong format!",argcheck = Check)
    Flag <- 1
  }
  if (is.null(path) || !is.character(path))
  {
    ArgumentCheck::addError(msg = "'path' is missing or in the wrong format!",argcheck = Check)
    Flag <- 1
  }
  if (is.null(name) || !is.character(name))
  {
    ArgumentCheck::addError(msg = "'name' is missing or in the wrong format!",argcheck = Check)
    Flag <- 1
  }
  ArgumentCheck::finishArgCheck(Check)
  if (is.null(Flag))
  {
    steps <- rep(1,ncol(model@DM))
    cat("Please define the step sizes of all parameters:")
    for (i in 1:length(steps))
    {
      steps[i] <- as.integer(readline(prompt = paste0(colnames(model@DM)[i],": ")))
    }
    grid_path <- paste(path,name,sep = "/")
    dir.create(grid_path,showWarnings = FALSE)
    saveRDS(model,file = paste0(grid_path,"/",name,".Gcfg"))
    ncores <- parallel::detectCores() - 1
    .Get_ParComb_cpp(model,grid_path,name,steps,ncores)
    pc_paths <- list.files(grid_path,full.names = TRUE,pattern = ".ParComb")
    out_paths <- pc_paths
    for (i in 1:length(out_paths))
    {
      out_paths[i] <- RSAGA::set.file.extension(out_paths[i],"GRID")
    }
    COMP_List <- list()
    for (i in 1:length(pc_paths))
    {
      COMP_List[[i]] <- list(model,pc_paths[i],out_paths[i])
    }
    clust <- parallel::makeCluster(ncores)
    parallel::clusterApply(clust, COMP_List, .Get_GRID_cpp)
    unlink(pc_paths)
  }
  else
  {
    return(cat("Get_Grid failed"))
  }
}
