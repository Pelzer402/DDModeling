#' Function to fit  a given DDRep to a given DDModel
#'
#' @name GetGrid
#'
Get_Grid <- function(model = NULL, wd = NULL){
  steps <- rep(1,ncol(model@DM))
  cat("Please define the step sizes of all parameters:")
  for (i in 1:length(steps))
  {
    steps[i] <- as.integer(readline(prompt = paste0(colnames(model@DM)[i],": ")))
  }
  grid_path <- paste(wd,"GRID",sep = "/")
  dir.create(grid_path,showWarnings = FALSE)
  ncores <- parallel::detectCores() - 1
  .Get_ParComb(model,grid_path,steps,ncores)
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
  parallel::clusterApply(clust, COMP_List, .Get_GRID)
  unlink(pc_paths)
}
