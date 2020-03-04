#' Generates a Grid from a given DDModel
#'
#' @name Generate_GRID
#' @rdname Generate_GRID
#' @param  model \code{DDModel} object
#' @param path \code{character} that specifies the full path to the directory in which the Grid should be saved
#' @param name \code{character} that represents the name (and subdirectory in path) of the Grid
#' @return No direct return value inside of the R-session. The calculated Grid will be saved in the specified path!
#' @description  After calling the function the user will be instructed to enter the step sizes corresponding to the parameters listet in the used model.
#'               Step size should allways be of an integer value, as they represent the number of evaluation points per parameter that are used in the grid.
#'               Note that in the given function the evaluation points are allways equally spaced concerning the corresponding parameter domain in the used model.
#'               Therefor, If one would like to specify the used evaluation points it is advised to specify the domain in the model.
#' @export
Generate_GRID <- function(model = NULL, path = NULL, name = NULL){
  Flag <- NULL
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model) || !methods::is(model,"DDModel"))
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
    pc_paths <- list.files(grid_path,full.names = TRUE,pattern = "\\.ParComb$")
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
    clust <- ParallelLogger::makeCluster(ncores)
    ParallelLogger::clusterApply(clust, COMP_List, .Get_GRID_cpp)
    unlink(pc_paths)
    ParallelLogger::stopCluster(clust)
    return(grid_path)
  }
  else
  {
    return(cat("Generate_GRID failed"))
  }
}

#' Imports a GRID into the R environment
#'
#' @name Import_GRID
#' @rdname Import_GRID
#' @param grid_path  \code{path} to a directory containing a .GRID fileset. If NULL the model will be fitted using 20 randomly drawn startparametersets from the model-DOMAIN.
#' @param to \code{character} that specifies the format to import to (choose from "frame", "keras_data" and "DDRep")
#' @return \code{list} of \code{DDRep}, \code{list} of \code{matrix} or \code{data.frame} dependend on the "to" parameter
#' @export
Import_GRID <- function(grid_path = NULL, to = "frame"){
  files <- list.files(grid_path,full.names = TRUE,pattern = "\\.GRID$")
  model <- readRDS(list.files(grid_path,full.names = TRUE,pattern = "\\.Gcfg$"))
  n_evals <- 0;
  conditions <- names(model@MM)
  CDF_perc <- as.character(model@RF$CDF)
  CAF_perc <- c()
  for (i in 1:(length(model@RF$CAF)-1))
  {
    buff <- (model@RF$CAF[i+1] + model@RF$CAF[i])/2
    CAF_perc <- c(CAF_perc,buff)
  }
  CAF_perc <- as.character(CAF_perc)
  PAR <- colnames(model@DM)
  colsum <- length(conditions)*(length(CDF_perc)*3+length(CAF_perc)*5+2) + length(PAR)
  CN <- c()
  for ( c in 1:length(conditions))
  {
    for (cdf_p in 1:length(CDF_perc))
    {
      CN <- c(CN,paste0(conditions[c],"_CDF_RT_",CDF_perc[cdf_p]))
      CN <- c(CN,paste0(conditions[c],"_CDF_PERC_",CDF_perc[cdf_p]))
      CN <- c(CN,paste0(conditions[c],"_CDF_N_",CDF_perc[cdf_p]))
    }
  }
  for ( c in 1:length(conditions))
  {
    for (caf_p in 1:length(CAF_perc))
    {
      CN <- c(CN,paste0(conditions[c],"_CAF_RT_",CAF_perc[caf_p]))
      CN <- c(CN,paste0(conditions[c],"_CAF_PERC_",CAF_perc[caf_p]))
      CN <- c(CN,paste0(conditions[c],"_CAF_ACC_",CAF_perc[caf_p]))
      CN <- c(CN,paste0(conditions[c],"_CAF_N_corr_",CAF_perc[caf_p]))
      CN <- c(CN,paste0(conditions[c],"_CAF_N_incorr_",CAF_perc[caf_p]))
    }
  }
  for ( c in 1:length(conditions))
  {
    CN <- c(CN,paste0(conditions[c],"_TOTAL_corr"),paste0(conditions[c],"_TOTAL_incorr"))
  }
  CN <- c(CN,PAR)
  IN <- as.data.frame(read.table(file = files[1],header = FALSE,skip = 1))
  for (i in 2:length(files))
  {
    n_evals <- n_evals + as.numeric(read.table(file = files[i],nrows = 1))
    buff <- as.data.frame(read.table(file = files[i],header = FALSE,skip = 1))
    IN <- rbind(IN,buff)
  }
  colnames(IN) <- CN
  if(to == "frame")
  {
    return(IN)
  }
  if(to == "keras_data")
  {
    CN_data <- c()
    for ( c in 1:length(conditions))
    {
      for (cdf_p in 1:length(CDF_perc))
      {
        CN_data <- c(CN_data,paste0(conditions[c],"_CDF_RT_",CDF_perc[cdf_p]))
      }
    }
    for ( c in 1:length(conditions))
    {
      for (caf_p in 1:length(CAF_perc))
      {
        CN_data <- c(CN_data,paste0(conditions[c],"_CAF_RT_",CAF_perc[caf_p]))
        CN_data <- c(CN_data,paste0(conditions[c],"_CAF_ACC_",CAF_perc[caf_p]))
      }
    }
    return(list(INPUT = as.matrix(IN[,CN_data]), OUTPUT = as.matrix(IN[,PAR])))
  }
  if (to == "DDRep")
  {
   return(.GRID_to_DDRep(list(model,files,n_evals)))
  }
}
