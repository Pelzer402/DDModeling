#' Function to simulate a given DDModel
#' @include DDModel.R Sim_DDModel.R
#' @name Sim_DDModel
#'
#' @param model \code{\linkS4class{DDModel}} object to be used in the simulation.
#' @param trials \code{Numeric} specifying the number of trials per condition.
#' @details For now this function will only initialize a given model with randomly drawn parameters from the domains specified in the DDModel object.
#' @examples M1 <- DDModel(model="DSTP",task = "flanker",
#'           CDF_perc = c(0.1,0.3,0.5,0.7,0.9),CAF_perc = c(0.0,0.2,0.4,0.6,0.8,1.0))
#' R1 <- Sim_DDModel(M1,10000)
#' @return \code{DDRep} object.
#' @export
Sim_DDModel <- function(model = NULL, trials = NULL){
  Flag <- 1
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model) || !methods::is(model,"DDModel"))
  {
    ArgumentCheck::addError(msg = "'model' is missing or in the wrong format!",argcheck = Check)
    Flag <- 99
  }
  if (is.null(trials) || !is.numeric(trials))
  {
    ArgumentCheck::addError(msg = "'trials' is missing or in the wrong format!",argcheck = Check)
    Flag <- 99
  }
  ArgumentCheck::finishArgCheck(Check)
  if (Flag == 99)
  {
    return(cat("Sim_DDModel failed"))
  }
  else
  {
    return(.Sim_DDModel(model,trials))
  }
}

