#' Function to simulate a given DDModel by initializing with random parameters
#'
#' @name Sim_DDModel
#'
#' @param model \code{DDModel} Object
#' @param trials \code{Numeric} specifying the number of trials per condition
#'
#' @return \code{DDRep} object
Sim_DDModel <- function(model = NULL, trials = NULL){
  Flag <- 1
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model) || !is(model,"DDModel"))
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

