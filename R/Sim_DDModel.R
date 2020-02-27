#' Function to simulate a given DDModel
#' @include DDModel.R Sim_DDModel.R
#' @name Sim_DDModel
#'
#' @param model \code{\linkS4class{DDModel}} object to be used in the simulation.
#' @param trials \code{Numeric} specifying the number of trials per condition.
#' @param simulations \code{Numeric} specifying the number of simulations
#' @param parameter \code{Data.frame} containing parameters that should be simulated
#' @details If 'parameter' is specified the function will undergo nrow(parameter) times simulation calls.
#' Depending on the input the return value will be either a \code{DDRep}, \code{list} or \code{list} of \code{lists}.
#' @examples # Define a Model
#' M1 <- DDModel(model="DSTP",task = "flanker",
#'           CDF_perc = c(0.1,0.3,0.5,0.7,0.9),CAF_perc = c(0.0,0.2,0.4,0.6,0.8,1.0))
#' # Simulate some Data with randomly generated parameters in the domain of M1
#' R1 <- Sim_DDModel(M1,10000)
#' # If you want to reuse the randomly drawn parameters (or any other) simply insert them
#' R2 <- Sim_DDModel(M1,10000,parameter = R1@PAR)
#' @return \code{DDRep} object or \code{list} of \code{DDRep} depending on the input
#' @export
Sim_DDModel <- function(model = NULL, trials = NULL, simulations = 1, parameter = NULL){
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
  if (!is.numeric(trials))
  {
    ArgumentCheck::addError(msg = "'simulations' is in the format of numeric!",argcheck = Check)
    Flag <- 99
  }
  if (!is.data.frame(parameter) && !is.null(parameter))
  {
    ArgumentCheck::addError(msg = "'parameter' is in the wrong format!",argcheck = Check)
    Flag <- 99
  }
  else if (!is.null(parameter))
  {
    if(all(colnames(model@DM) == colnames(parameter)))
    {
      Flag <- 2
    }
    else
    {
      ArgumentCheck::addError(msg = "'parameter' do not match the definition in the 'model'!",argcheck = Check)
      Flag <- 99
    }
  }
  ArgumentCheck::finishArgCheck(Check)
  if (Flag == 99)
  {
    return(cat("Sim_DDModel failed"))
  }
  else
  {
    if (Flag == 1)
    {
      if (simulations == 1)
      {
        return(.Sim_DDModel_rnd(model,trials))
      }
      else
      {
        OUT <- list()
        for (i in 1:simulations)
        {
          OUT[[i]] <- .Sim_DDModel_rnd(model,trials)
        }
        return(OUT)
      }
    }
    if (Flag == 2)
    {
      if (simulations == 1)
      {
        if (nrow(parameter) == 1)
        {
          par_buff <- as.numeric(parameter[1,])
          return(.Sim_DDModel_par(model,trials,par_buff))
        }
        else
        {
          OUT <- list()
          for (p in 1:nrow(parameter))
          {
            par_buff <- as.numeric(parameter[p,])
            OUT[[p]] <- .Sim_DDModel_par(model,trials,par_buff)
          }
          return(OUT)
        }
      }
      else
      {
        OUT <- list()
        if (nrow(parameter) == 1)
        {
          for (s in 1:simulations)
          {
            par_buff <- as.numeric(parameter[1,])
            OUT[[s]] <- .Sim_DDModel_par(model,trials,par_buff)
          }
          return(OUT)
        }
        else
        {
          for (p in 1:nrow(parameter))
          {
            OUT[[p]] <- list()
            for (s in 1:simulations)
            {
              par_buff <- as.numeric(parameter[p,])
              OUT[[p]][[s]] <- .Sim_DDModel_par(model,trials,par_buff)
            }
          }
          return(OUT)
        }
      }
    }
  }
}

