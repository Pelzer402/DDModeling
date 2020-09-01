#' Generates a Grid from a given DDModel
#'
#' @name Generate_GRID
#' @rdname Generate_GRID
#' @param  model \code{DDModel} object
#' @param path \code{character} that specifies the full path to the directory in which the Grid should be saved
#' @param name \code{character} that represents the name of the Grid (and indirectly the name of the subdirectory set in the directory specified in path)
#' @param eval_pts \code{numeric} vector containing number of evaluation points of the parameters in the grid
#' @return The path to the grid as a \code{character}
#' @details   Eval_pts specifies the number of evaluation points per parameter that are used in the grid, such the number of integers in
#' eval_pts needs to match the number of parameters in the model. Note that in the given function the evaluation points are allways equally spaced concerning the corresponding parameter domain in the used model.
#' Therefor, if one would like to specify the used evaluation points it is advised to specify the domain in the model.
#' @export
Generate_GRID <- function(model = NULL, path = NULL, name = NULL, eval_pts = NULL) {
  Flag <- NULL
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model) || !methods::is(model, "DDModel")) {
    ArgumentCheck::addError(msg = "'model' is missing or in the wrong format!", argcheck = Check)
    Flag <- 1
  }
  if (is.null(path) || !is.character(path)) {
    ArgumentCheck::addError(msg = "'path' is missing or in the wrong format!", argcheck = Check)
    Flag <- 1
  }
  if (is.null(name) || !is.character(name)) {
    ArgumentCheck::addError(msg = "'name' is missing or in the wrong format!", argcheck = Check)
    Flag <- 1
  }
  if (!all(eval_pts %% 1 == 0) || ncol(model@DM) != length(eval_pts) || is.null(eval_pts)) {
    ArgumentCheck::addError(msg = "'eval_pts' is missing or in the wrong format!", argcheck = Check)
    Flag <- 1
  }
  ArgumentCheck::finishArgCheck(Check)
  if (is.null(Flag)) {
    steps <- eval_pts
    grid_path <- paste(path, name, sep = "/")
    dir.create(grid_path, showWarnings = FALSE)
    saveRDS(model, file = paste0(grid_path, "/", name, ".Gcfg"))
    ncores <- parallel::detectCores() - 1
    .Get_ParComb_cpp(model, grid_path, name, steps, ncores)
    pc_paths <- list.files(grid_path, full.names = TRUE, pattern = "\\.ParComb$")
    out_paths <- pc_paths
    for (i in 1:length(out_paths))
    {
      out_paths[i] <- RSAGA::set.file.extension(out_paths[i], "GRID")
    }
    COMP_List <- list()
    for (i in 1:length(pc_paths))
    {
      COMP_List[[i]] <- list(model, pc_paths[i], out_paths[i])
    }
    clust <- ParallelLogger::makeCluster(ncores)
    ParallelLogger::clusterApply(clust, COMP_List, .Get_GRID_cpp)
    unlink(pc_paths)
    ParallelLogger::stopCluster(clust)
    return(grid_path)
  }
  else {
    return(cat("Generate_GRID failed"))
  }
}

#' Imports a GRID into the R environment
#'
#' @name Import_GRID
#' @rdname Import_GRID
#' @param grid_path  \code{path} to a directory containing a .GRID fileset.
#' @param to \code{character} that specifies the format to import to (choose from "frame", "keras_data" and "DDRep")
#' @details  "frame" will return a \code{data.frame} like structure to the total Grid stored in the grid_path. "keras_data" will return a \code{list}
#' containing an INPUT and OUTPUT element. The INPUT constitutes all relevant values from the CDF/CAF (namely reaction times and accuracy) in a \code{matrix} with each row
#' corresponding to one parameter set in the OUTPUT. "DDRep" will convert the complete Grid to DDReps in a \code{list}.
#' @return \code{list} of \code{DDRep}, \code{list} of \code{matrix} or \code{data.frame} dependend on the "to" parameter
#' @export
Import_GRID <- function(grid_path = NULL, to = "frame") {
  files <- list.files(grid_path, full.names = TRUE, pattern = "\\.GRID$")
  model <- readRDS(list.files(grid_path, full.names = TRUE, pattern = "\\.Gcfg$"))
  n_evals <- 0
  conditions <- names(model@MM)
  CDF_perc <- as.character(model@RF$CDF)
  CAF_perc <- c()
  for (i in 1:(length(model@RF$CAF) - 1))
  {
    buff <- (model@RF$CAF[i + 1] + model@RF$CAF[i]) / 2
    CAF_perc <- c(CAF_perc, buff)
  }
  CAF_perc <- as.character(CAF_perc)
  PAR <- colnames(model@DM)
  colsum <- length(conditions) * (length(CDF_perc) * 3 + length(CAF_perc) * 5 + 2) + length(PAR)
  CN <- c()
  for (c in 1:length(conditions))
  {
    for (cdf_p in 1:length(CDF_perc))
    {
      CN <- c(CN, paste0(conditions[c], "_CDF_RT_", CDF_perc[cdf_p]))
      CN <- c(CN, paste0(conditions[c], "_CDF_N_", CDF_perc[cdf_p]))
    }
  }
  for (c in 1:length(conditions))
  {
    for (caf_p in 1:length(CAF_perc))
    {
      CN <- c(CN, paste0(conditions[c], "_CAF_RT_", CAF_perc[caf_p]))
      CN <- c(CN, paste0(conditions[c], "_CAF_ACC_", CAF_perc[caf_p]))
      CN <- c(CN, paste0(conditions[c], "_CAF_N_corr_", CAF_perc[caf_p]))
      CN <- c(CN, paste0(conditions[c], "_CAF_N_incorr_", CAF_perc[caf_p]))
    }
  }
  CN <- c(CN, PAR)
  n_evals <- 0
  for (i in 1:length(files))
  {
    n_evals <- n_evals + as.numeric(read.table(file = files[i], nrows = 1))
  }
  if (to == "frame") {
    IN <- list()
    for (i in 1:length(files))
    {
      IN[[i]] <- data.table::fread(file = files[i], header = FALSE, skip = 1)
    }
    IN <- data.table::rbindlist(IN)
    colnames(IN) <- CN
    return(as.data.frame(IN))
  }
  if (to == "keras_data") {
    CN_data <- c()
    for (c in 1:length(conditions))
    {
      for (cdf_p in 1:length(CDF_perc))
      {
        CN_data <- c(CN_data, paste0(conditions[c], "_CDF_RT_", CDF_perc[cdf_p]))
      }
    }
    for (c in 1:length(conditions))
    {
      for (caf_p in 1:length(CAF_perc))
      {
        CN_data <- c(CN_data, paste0(conditions[c], "_CAF_RT_", CAF_perc[caf_p]))
        CN_data <- c(CN_data, paste0(conditions[c], "_CAF_ACC_", CAF_perc[caf_p]))
      }
    }
    IN <- list()
    OUT <- list()
    for (i in 1:length(files))
    {
      buff <- data.table::fread(file = files[i], header = FALSE, skip = 1)
      colnames(buff) <- CN
      IN[[i]] <- buff[, ..CN_data]
      OUT[[i]] <- buff[, ..PAR]
      rm(buff)
    }
    IN <- data.table::rbindlist(IN)
    OUT <- data.table::rbindlist(OUT)
    return(list(INPUT = as.matrix(IN), OUTPUT = as.matrix(OUT)))
  }
  if (to == "DDRep") {
    return(.GRID_to_DDRep(list(model, files, n_evals)))
  }
}
