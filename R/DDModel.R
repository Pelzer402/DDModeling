#' An S4 class to represent a drift diffusion model
#'
#' @name DDModel-class
#' @rdname DDModel-class
#' @slot ID \code{character} that represents the name of the model to be used (i.e. "DSTP","DMC","SSP")
#' @slot MM   \code{list} of \code{matrix} that contain values which map custom parameters to correspondend modelparameters
#' @slot DM  \code{matrix} that contains the domain of all custom parameters (and grid size steps)
#' @slot SP  \code{matrix} that contains a set of simulation-parameters important for simulation
#' @slot RF \code{list} of numeric vectors that contain the percentiles of the representation.
setClass("DDModel",
         slots      = list(
           ID       = "character",
           MM       = "list",
           DM       = "matrix",
           SP       = "matrix",
           RF       = "list"
         )
)


#' Function to generate a DDModel S4 class object
#'
#' @name DDModel
#'
#' @param model \code{character} of the name of the Model to be used (legitamite choices are "DSTP","DMC","SSP")
#' @param task \code{character} specifying a specific predefined modelstructure ("flanker")
#' @param conditions \code{character} vector of the names of conditions
#' @param parameter \code{character} vector of the names of custom parameters
#' @param dt \code{numeric} representing the integration constant of the diffusion process
#' @param sigma \code{numeric} representing the diffusion constant of the diffusion process
#' @param CDF_perc \code{Numeric} vector specifying the CDF percentiles (note: numbers equal to absolut percentiles!)
#' @param CAF_perc \code{Numeric}  vector specifying the CAF percentiles (note: numbers equal to boarders of segments!)
#'
#' @return \code{DDMODEL} object
DDModel <- function(model = NULL,task = NULL,conditions=NULL,parameter=NULL,dt=NULL,sigma = NULL,CDF_perc = NULL, CAF_perc = NULL){
  Flag <- NULL
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model))
  {
    ArgumentCheck::addError(msg = "'model' is missing!",argcheck = Check)
  }
  else
  {
    if (!(model %in% c("DSTP","DMC","SSP")))
      ArgumentCheck::addError(msg = "'model' must be one of 'DSTP', 'DMC' or 'SSP'!",argcheck = Check)
  }
  if (is.null(task))
    ArgumentCheck::addError(msg = "'task' is missing!",argcheck = Check)
  if (!is.null(task))
  {
    if(!(task %in% c("flanker","custom")))
      ArgumentCheck::addError(msg = "'task' must be one of 'flanker' or 'custom'!",argcheck = Check)
    if(task %in% c("flanker")&& ((!is.null(conditions))||(!is.null(parameter))||(!is.null(dt))||(!is.null(sigma))))
    {
      if(!is.null(conditions))
        ArgumentCheck::addMessage(msg = "'conditions' will be discarded because of task specification",argcheck = Check)
      if(!is.null(parameter))
        ArgumentCheck::addMessage(msg = "'parameter' will be discarded because of task specification",argcheck = Check)
      if(!is.null(dt))
        ArgumentCheck::addWarning(msg = "'dt' will be discarded because of task specification",argcheck = Check)
      if(!is.null(sigma))
        ArgumentCheck::addWarning(msg = "'sigma' will be discarded because of task specification",argcheck = Check)
      Flag <- "flanker"
    }
    else
    {
      Flag <- "flanker"
    }
    if(task %in% c("custom"))
    {
      i <- 0
      if(is.null(conditions))
      {
        ArgumentCheck::addError(msg = "'conditions' is missing!",argcheck = Check)
      }
      else if(!is.vector(conditions,mode = "character"))
      {
        ArgumentCheck::addError(msg = "'conditions' is not a character vector!",argcheck = Check)
      }
      else
      {
        i <- i+1
      }
      if(is.null(parameter))
      {
        ArgumentCheck::addError(msg = "'parameter' is missing!",argcheck = Check)
      }
      else if(!is.vector(parameter,mode = "character"))
      {
        ArgumentCheck::addError(msg = "'parameter' is not a character vector!",argcheck = Check)
      }
      else
      {
        i <- i+1
      }
      if(is.null(dt))
      {
        ArgumentCheck::addError(msg = "'dt' is missing!",argcheck = Check)
      }
      else if(!is.numeric(dt))
      {
        ArgumentCheck::addError(msg = "'dt' is not a numeric!",argcheck = Check)
      }
      else
      {
        i <- i+1
      }
      if(is.null(sigma))
      {
        ArgumentCheck::addError(msg = "'sigma' is missing!",argcheck = Check)
      }
      else  if(!is.numeric(sigma))
      {
        ArgumentCheck::addError(msg = "'sigma' is not a numeric!",argcheck = Check)
      }
      else
      {
        i <- i+1
      }
      if(i==4)
        Flag <- "custom"
    }
  }
  ArgumentCheck::finishArgCheck(Check)
  if(is.null(Flag))
    return(cat("DDModel failed"))
  switch(Flag,
         custom ={
           switch(EXPR = model,
                  DSTP={
                    mm <- lapply(1:length(conditions), function(x){ x<-matrix(0,nrow = 7,ncol = length(parameter))
                    rownames(x) <- c("Ter","a","c","mu_RS1","mu_RS2_C","mu_RS2_D","mu_SS")
                    colnames(x) <- parameter
                    return(x)
                    })
                    names(mm) <- conditions
                    dm <- matrix(nrow = 2,ncol = length(parameter))
                    rownames(dm) <- c("Upper_Limit","Lower_Limit")
                    colnames(dm) <- parameter
                  },
                  DMC={
                    mm <- lapply(1:length(conditions), function(x){ x<-matrix(0,nrow = 6,ncol = length(parameter))
                    rownames(x) <- c("Ter","a","zeta","alpha","mu_c","tau")
                    colnames(x) <- parameter
                    return(x)
                    })
                    names(mm) <- conditions
                    dm <- matrix(nrow = 2,ncol = length(parameter))
                    rownames(dm) <- c("Upper_Limit","Lower_Limit")
                    colnames(dm) <- parameter
                  },
                  SSP={
                    mm <- lapply(1:length(conditions), function(x){ x<-matrix(0,nrow = 5,ncol = length(parameter))
                    rownames(x) <- c("Ter","a","P","sda","rd")
                    colnames(x) <- parameter
                    return(x)
                    })
                    names(mm) <- conditions
                    dm <- matrix(nrow = 2,ncol = length(parameter))
                    rownames(dm) <- c("Upper_Limit","Lower_Limit")
                    colnames(dm) <- parameter
                  }
           )
           sp <- as.matrix(data.frame(dt=dt,sigma=sigma))
           return(methods::new("DDModel",ID=model,MM=mm,DM=dm,SP=sp))
         },
         flanker={
           switch(EXPR = model,
                  DSTP={
                    conditions <- c("Cong","Incong")
                    parameter <- c("Ter","a","c","mu_t","mu_f","mu_RS2","mu_SS")
                    dt <- 0.001
                    sigma <- 0.1
                    mm <- lapply(1:length(conditions), function(x){ x<-matrix(0,nrow = 7,ncol = length(parameter))
                    rownames(x) <- c("Ter","a","c","mu_RS1","mu_RS2_C","mu_RS2_D","mu_SS")
                    colnames(x) <- parameter
                    return(x)
                    })
                    names(mm) <- conditions
                    dm <- matrix(nrow = 2,ncol = length(parameter))
                    rownames(dm) <- c("Upper_Limit","Lower_Limit")
                    colnames(dm) <- parameter
                    mm$Cong["Ter",]        <- c(1,0,0,0,0,0,0)
                    mm$Cong["a",]          <- c(0,1,0,0,0,0,0)
                    mm$Cong["c",]          <- c(0,0,1,0,0,0,0)
                    mm$Cong["mu_RS1",]     <- c(0,0,0,1,1,0,0)
                    mm$Cong["mu_RS2_C",]   <- c(0,0,0,0,0,1,0)
                    mm$Cong["mu_RS2_D",]   <- c(0,0,0,0,0,1,0)
                    mm$Cong["mu_SS",]      <- c(0,0,0,0,0,0,1)
                    mm$Incong["Ter",]      <- c(1,0,0,0,0,0,0)
                    mm$Incong["a",]        <- c(0,1,0,0,0,0,0)
                    mm$Incong["c",]        <- c(0,0,1,0,0,0,0)
                    mm$Incong["mu_RS1",]   <- c(0,0,0,1,-1,0,0)
                    mm$Incong["mu_RS2_C",] <- c(0,0,0,0,0,1,0)
                    mm$Incong["mu_RS2_D",] <- c(0,0,0,0,0,-1,0)
                    mm$Incong["mu_SS",]    <- c(0,0,0,0,0,0,1)
                    dm["Upper_Limit",] <- c(0.45,0.38,0.38,0.15,0.25,0.55,1.2)
                    dm["Lower_Limit",] <- c(0.15,0.14,0.14,0.05,0.05,0.25,0.4)
                  },
                  DMC={
                    conditions <- c("Cong","Incong")
                    parameter <- c("Ter","a","zeta","alpha","mu_c","tau")
                    dt <- 0.001
                    sigma <- 4
                    mm <- lapply(1:length(conditions), function(x){ x<-matrix(0,nrow = 6,ncol = length(parameter))
                    rownames(x) <- c("Ter","a","zeta","alpha","mu_c","tau")
                    colnames(x) <- parameter
                    return(x)
                    })
                    names(mm) <- conditions
                    dm <- matrix(nrow = 2,ncol = length(parameter))
                    rownames(dm) <- c("Upper_Limit","Lower_Limit")
                    colnames(dm) <- parameter
                    mm$Cong["Ter",]   <- c(1,0,0,0,0,0)
                    mm$Cong["a",]     <- c(0,1,0,0,0,0)
                    mm$Cong["zeta",]  <- c(0,0,1,0,0,0)
                    mm$Cong["alpha",] <- c(0,0,0,1,0,0)
                    mm$Cong["mu_c",]  <- c(0,0,0,0,1,0)
                    mm$Cong["tau",]   <- c(0,0,0,0,0,1)
                    mm$Incong["Ter",]   <- c(1,0,0,0,0,0)
                    mm$Incong["a",]     <- c(0,1,0,0,0,0)
                    mm$Incong["zeta",]  <- c(0,0,-1,0,0,0)
                    mm$Incong["alpha",] <- c(0,0,0,1,0,0)
                    mm$Incong["mu_c",]  <- c(0,0,0,0,1,0)
                    mm$Incong["tau",]   <- c(0,0,0,0,0,1)
                    dm["Upper_Limit",] <- c(400,160,40,4.5,0.8,120)
                    dm["Lower_Limit",] <- c(270,90,15,1.5,0.2,20)
                  },
                  SSP={
                    conditions <- c("Cong","Incong")
                    parameter <- c("Ter","a","P","sda","rd")
                    dt <- 0.001
                    sigma <- 0.1
                    mm <- lapply(1:length(conditions), function(x){ x<-matrix(0,nrow = 5,ncol = length(parameter))
                    rownames(x) <- c("Ter","a","P","sda","rd")
                    colnames(x) <- parameter
                    return(x)
                    })
                    names(mm) <- conditions
                    dm <- matrix(nrow = 2,ncol = length(parameter))
                    rownames(dm) <- c("Upper_Limit","Lower_Limit")
                    colnames(dm) <- parameter
                    mm$Cong["Ter",]   <- c(1,0,0,0,0)
                    mm$Cong["a",]     <- c(0,1,0,0,0)
                    mm$Cong["P",]     <- c(0,0,1,0,0)
                    mm$Cong["sda",]   <- c(0,0,0,1,0)
                    mm$Cong["rd",]    <- c(0,0,0,0,1)
                    mm$Incong["Ter",] <- c(1,0,0,0,0)
                    mm$Incong["a",]   <- c(0,1,0,0,0)
                    mm$Incong["P",]   <- c(0,0,-1,0,0)
                    mm$Incong["sda",] <- c(0,0,0,1,0)
                    mm$Incong["rd",]   <- c(0,0,0,0,1)
                    dm["Upper_Limit",] <- c(0.45,0.19,0.55,2.6,0.026)
                    dm["Lower_Limit",] <- c(0.15,0.07,0.2,1,0.01)
                  }
           )
           sp <- as.matrix(data.frame(dt=dt,sigma=sigma))
           rf <- list(CDF=CDF_perc,CAF=CAF_perc)
           return(methods::new("DDModel",ID=model,MM=mm,DM=dm,SP=sp,RF = rf))
         }
  )
}

setMethod("show","DDModel",function(object){
  cat("S4 DD.MODEL Object: \n\n",
      "Model: ", object@ID,"\n\n",
      "ModelMatrix: \n")
  print(object@MM)
  cat(" Parameter Domain: \n")
  print(object@DM)
  cat("\n Simulation Parameter: \n")
  print(object@SP)
  cat("\n Form of Representation: \n")
  print(object@RF)
})


