#' \code{DDModel} class definition
#' @name DDModel-class
#' @rdname DDModel-class
#' @slot ID   \code{character} that represents the name of the model to be used (i.e. "DSTP","DMC","SSP")
#' @slot MM   \code{list} of \code{matrix} that contain values which map custom parameters to correspondend modelparameters
#' @slot DM   \code{matrix} that contains the domain of all custom parameters (and grid size steps)
#' @slot SP   \code{matrix} that contains a set of simulation-parameters important for simulation
#' @slot RF   \code{list} of numeric vectors that contain the percentiles of the representation.
setClass("DDModel",
  slots = list(
    ID = "character",
    MM = "list",
    DM = "matrix",
    SP = "matrix",
    RF = "list"
  )
)


#' Constructor for \link{DDModel-class}
#' @include DDModel.R
#' @name DDModel
#' @description Userfriendly function for the construction of a \link{DDModel-class}.
#' @param model       \code{character} of the name of the Model to be used (choices are "DSTP","DMC","SSP")
#' @param task        \code{character} specifying a specific predefined modelstructure ("flanker","custom","RMT_LDT")
#' @param conditions  \code{character} vector of the names of conditions
#' @param parameter   \code{character} vector of the names of custom parameters
#' @param dt          \code{numeric} representing the integration constant of the diffusion process
#' @param sigma       \code{numeric} representing the diffusion constant of the diffusion process
#' @param CDF_perc    \code{Numeric} vector specifying the CDF percentiles (note: numbers equal to absolut percentiles!)
#' @param CAF_perc    \code{Numeric} vector specifying the CAF percentiles (note: numbers equal to boarders of segments!)
#' @return \code{\link{DDModel-class}}.
#' @details The constructor allows for the usage of handy defaults for specific tasks. Choosing task = "flanker" will configure the corresponding model in a way, that is
#' in correspondance to present litretature. Only the type of data representation (i.e. CDF and CAF quantiles) need further specification. If however, one wants to define
#' customized applications of the models, the task = "custom" preset can be choosen. It is important to emphazise that in this case there are several more parameters needed
#' in the constructor (see usage or examples). Moreso after model construction one must specify the parameter domains and model matrizes manually! Be sure to double check
#' your model befor computation. Custom models should only be used, if one knows exactly what he is doing!
#' @examples
#' M_DSTP <- DDModel(
#'   model = "DSTP", task = "flanker",
#'   CDF_perc = c(0.1, 0.3, 0.5, 0.7, 0.9), CAF_perc = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
#' )
#' M_DMC <- DDModel(
#'   model = "DMC", task = "flanker",
#'   CDF_perc = c(0.1, 0.3, 0.5, 0.7, 0.9), CAF_perc = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
#' )
#' M_Custom <- DDModel(
#'   model = "DSTP", task = "custom",
#'   CDF_perc = c(0.1, 0.3, 0.5, 0.7, 0.9), CAF_perc = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
#'   conditions = "single", parameter = c("par_1", "par_2"), dt = 0.001, sigma = 0.1
#' )
#' @export
DDModel <- function(model = NULL, task = NULL, conditions = NULL, parameter = NULL, dt = NULL, sigma = NULL, CDF_perc = NULL, CAF_perc = NULL) {
  Flag <- NULL
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model)) {
    ArgumentCheck::addError(msg = "'model' is missing!", argcheck = Check)
  }
  else {
    if (!(model %in% c("DSTP", "DMC", "SSP", "DDM_classic"))) {
      ArgumentCheck::addError(msg = "'model' must be one of 'DSTP', 'DMC', 'SSP' or 'DDM_classic'!", argcheck = Check)
    }
  }
  if (is.null(task)) {
    ArgumentCheck::addError(msg = "'task' is missing!", argcheck = Check)
  }
  if (!is.null(task)) {
    if (!(task %in% c("flanker", "custom", "RMT_LDT"))) {
      ArgumentCheck::addError(msg = "'task' must be one of 'flanker', 'custom' or 'RMT_LDT'!", argcheck = Check)
    }
    if (task %in% c("flanker", "RMT_LDT") && ((!is.null(conditions)) || (!is.null(parameter)) || (!is.null(dt)) || (!is.null(sigma)))) {
      if (!is.null(conditions)) {
        ArgumentCheck::addMessage(msg = "'conditions' will be discarded because of task specification", argcheck = Check)
      }
      if (!is.null(parameter)) {
        ArgumentCheck::addMessage(msg = "'parameter' will be discarded because of task specification", argcheck = Check)
      }
      if (!is.null(dt)) {
        ArgumentCheck::addWarning(msg = "'dt' will be discarded because of task specification", argcheck = Check)
      }
      if (!is.null(sigma)) {
        ArgumentCheck::addWarning(msg = "'sigma' will be discarded because of task specification", argcheck = Check)
      }
    }
    if (task %in% c("flanker")) {
      Flag <- "flanker"
    }
    if (task %in% c("RMT_LDT")) {
      Flag <- "RMT_LDT"
    }
    if (task %in% c("custom")) {
      i <- 0
      if (is.null(conditions)) {
        ArgumentCheck::addError(msg = "'conditions' is missing!", argcheck = Check)
      }
      else if (!is.vector(conditions, mode = "character")) {
        ArgumentCheck::addError(msg = "'conditions' is not a character vector!", argcheck = Check)
      }
      else {
        i <- i + 1
      }
      if (is.null(parameter)) {
        ArgumentCheck::addError(msg = "'parameter' is missing!", argcheck = Check)
      }
      else if (!is.vector(parameter, mode = "character")) {
        ArgumentCheck::addError(msg = "'parameter' is not a character vector!", argcheck = Check)
      }
      else {
        i <- i + 1
      }
      if (is.null(dt)) {
        ArgumentCheck::addError(msg = "'dt' is missing!", argcheck = Check)
      }
      else if (!is.numeric(dt)) {
        ArgumentCheck::addError(msg = "'dt' is not a numeric!", argcheck = Check)
      }
      else {
        i <- i + 1
      }
      if (is.null(sigma)) {
        ArgumentCheck::addError(msg = "'sigma' is missing!", argcheck = Check)
      }
      else if (!is.numeric(sigma)) {
        ArgumentCheck::addError(msg = "'sigma' is not a numeric!", argcheck = Check)
      }
      else {
        i <- i + 1
      }
      if (i == 4) {
        Flag <- "custom"
      }
    }
  }
  ArgumentCheck::finishArgCheck(Check)
  if (is.null(Flag)) {
    return(cat("DDModel failed"))
  }
  switch(Flag,
    custom = {
      switch(EXPR = model,
        DSTP = {
          mm <- lapply(1:length(conditions), function(x) {
            x <- matrix(0, nrow = 16, ncol = length(parameter))
            rownames(x) <- c("Ter", "a", "c", "mu_RS1", "mu_RS2_C", "mu_RS2_D", "mu_SS", "z_RS", "z_SS", "s_Ter", "s_mu_RS1", "s_mu_RS2_C", "s_mu_RS2_D", "s_SS", "s_z_RS", "s_z_SS")
            colnames(x) <- parameter
            return(x)
          })
          names(mm) <- conditions
          dm <- matrix(nrow = 2, ncol = length(parameter))
          rownames(dm) <- c("Upper_Limit", "Lower_Limit")
          colnames(dm) <- parameter
        },
        DMC = {
          mm <- lapply(1:length(conditions), function(x) {
            x <- matrix(0, nrow = 9, ncol = length(parameter))
            rownames(x) <- c("Ter", "a", "zeta", "alpha", "mu_c", "tau", "z", "s_Ter", "s_z")
            colnames(x) <- parameter
            return(x)
          })
          names(mm) <- conditions
          dm <- matrix(nrow = 2, ncol = length(parameter))
          rownames(dm) <- c("Upper_Limit", "Lower_Limit")
          colnames(dm) <- parameter
        },
        SSP = {
          mm <- lapply(1:length(conditions), function(x) {
            x <- matrix(0, nrow = 8, ncol = length(parameter))
            rownames(x) <- c("Ter", "a", "P", "sda", "rd", "z", "s_Ter", "s_z")
            colnames(x) <- parameter
            return(x)
          })
          names(mm) <- conditions
          dm <- matrix(nrow = 2, ncol = length(parameter))
          rownames(dm) <- c("Upper_Limit", "Lower_Limit")
          colnames(dm) <- parameter
        },
        DDM_classic = {
          mm <- lapply(1:length(conditions), function(x) {
            x <- matrix(0, nrow = 7, ncol = length(parameter))
            rownames(x) <- c("Ter", "a", "mu", "z", "s_Ter", "s_mu", "s_z")
            colnames(x) <- parameter
            return(x)
          })
          names(mm) <- conditions
          dm <- matrix(nrow = 2, ncol = length(parameter))
          rownames(dm) <- c("Upper_Limit", "Lower_Limit")
          colnames(dm) <- parameter
        }
      )
      sp <- as.matrix(data.frame(dt = dt, sigma = sigma))
      rf <- list(CDF = CDF_perc, CAF = CAF_perc)
      return(methods::new("DDModel", ID = model, MM = mm, DM = dm, SP = sp, RF = rf))
    },
    flanker = {
      switch(EXPR = model,
        DSTP = {
          conditions <- c("Cong", "Incong")
          parameter <- c("Ter", "a", "c", "mu_t", "mu_f", "mu_RS2", "mu_SS")
          dt <- 0.001
          sigma <- 0.1
          mm <- lapply(1:length(conditions), function(x) {
            x <- matrix(0, nrow = 16, ncol = length(parameter))
            rownames(x) <- c("Ter", "a", "c", "mu_RS1", "mu_RS2_C", "mu_RS2_D", "mu_SS", "z_RS", "z_SS", "s_Ter", "s_mu_RS1", "s_mu_RS2_C", "s_mu_RS2_D", "s_SS", "s_z_RS", "s_z_SS")
            colnames(x) <- parameter
            return(x)
          })
          names(mm) <- conditions
          dm <- matrix(nrow = 2, ncol = length(parameter))
          rownames(dm) <- c("Upper_Limit", "Lower_Limit")
          colnames(dm) <- parameter
          mm$Cong["Ter", ] <- c(1, 0, 0, 0, 0, 0, 0)
          mm$Cong["a", ] <- c(0, 1, 0, 0, 0, 0, 0)
          mm$Cong["c", ] <- c(0, 0, 1, 0, 0, 0, 0)
          mm$Cong["mu_RS1", ] <- c(0, 0, 0, 1, 1, 0, 0)
          mm$Cong["mu_RS2_C", ] <- c(0, 0, 0, 0, 0, 1, 0)
          mm$Cong["mu_RS2_D", ] <- c(0, 0, 0, 0, 0, 1, 0)
          mm$Cong["mu_SS", ] <- c(0, 0, 0, 0, 0, 0, 1)
          mm$Cong["z_RS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Cong["z_SS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Cong["s_Ter", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Cong["s_mu_RS1", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Cong["s_mu_RS2_C", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Cong["s_mu_RS2_D", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Cong["s_SS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Cong["s_z_RS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Cong["s_z_SS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Incong["Ter", ] <- c(1, 0, 0, 0, 0, 0, 0)
          mm$Incong["a", ] <- c(0, 1, 0, 0, 0, 0, 0)
          mm$Incong["c", ] <- c(0, 0, 1, 0, 0, 0, 0)
          mm$Incong["mu_RS1", ] <- c(0, 0, 0, 1, -1, 0, 0)
          mm$Incong["mu_RS2_C", ] <- c(0, 0, 0, 0, 0, 1, 0)
          mm$Incong["mu_RS2_D", ] <- c(0, 0, 0, 0, 0, -1, 0)
          mm$Incong["mu_SS", ] <- c(0, 0, 0, 0, 0, 0, 1)
          mm$Incong["z_RS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Incong["z_SS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Incong["s_Ter", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Incong["s_mu_RS1", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Incong["s_mu_RS2_C", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Incong["s_mu_RS2_D", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Incong["s_SS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Incong["s_z_RS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          mm$Incong["s_z_SS", ] <- c(0, 0, 0, 0, 0, 0, 0)
          dm["Upper_Limit", ] <- c(0.45, 0.38, 0.38, 0.15, 0.25, 0.55, 1.2)
          dm["Lower_Limit", ] <- c(0.15, 0.14, 0.14, 0.05, 0.05, 0.25, 0.4)
        },
        DMC = {
          conditions <- c("Cong", "Incong")
          parameter <- c("Ter", "a", "zeta", "alpha", "mu_c", "tau")
          dt <- 0.001
          sigma <- 4
          mm <- lapply(1:length(conditions), function(x) {
            x <- matrix(0, nrow = 9, ncol = length(parameter))
            rownames(x) <- c("Ter", "a", "zeta", "alpha", "mu_c", "tau", "z", "s_Ter", "s_z")
            colnames(x) <- parameter
            return(x)
          })
          names(mm) <- conditions
          dm <- matrix(nrow = 2, ncol = length(parameter))
          rownames(dm) <- c("Upper_Limit", "Lower_Limit")
          colnames(dm) <- parameter
          mm$Cong["Ter", ] <- c(1, 0, 0, 0, 0, 0)
          mm$Cong["a", ] <- c(0, 1, 0, 0, 0, 0)
          mm$Cong["zeta", ] <- c(0, 0, 1, 0, 0, 0)
          mm$Cong["alpha", ] <- c(0, 0, 0, 1, 0, 0)
          mm$Cong["mu_c", ] <- c(0, 0, 0, 0, 1, 0)
          mm$Cong["tau", ] <- c(0, 0, 0, 0, 0, 1)
          mm$Cong["z", ] <- c(0, 0, 0, 0, 0, 0)
          mm$Cong["s_Ter", ] <- c(0, 0, 0, 0, 0, 0)
          mm$Cong["s_z", ] <- c(0, 0, 0, 0, 0, 0)
          mm$Incong["Ter", ] <- c(1, 0, 0, 0, 0, 0)
          mm$Incong["a", ] <- c(0, 1, 0, 0, 0, 0)
          mm$Incong["zeta", ] <- c(0, 0, -1, 0, 0, 0)
          mm$Incong["alpha", ] <- c(0, 0, 0, 1, 0, 0)
          mm$Incong["mu_c", ] <- c(0, 0, 0, 0, 1, 0)
          mm$Incong["tau", ] <- c(0, 0, 0, 0, 0, 1)
          mm$Incong["z", ] <- c(0, 0, 0, 0, 0, 0)
          mm$Incong["s_Ter", ] <- c(0, 0, 0, 0, 0, 0)
          mm$Incong["s_z", ] <- c(0, 0, 0, 0, 0, 0)
          dm["Upper_Limit", ] <- c(400, 160, 40, 4.5, 0.8, 120)
          dm["Lower_Limit", ] <- c(270, 90, 15, 1.5, 0.2, 20)
        },
        SSP = {
          conditions <- c("Cong", "Incong")
          parameter <- c("Ter", "a", "P", "sda", "rd")
          dt <- 0.001
          sigma <- 0.1
          mm <- lapply(1:length(conditions), function(x) {
            x <- matrix(0, nrow = 8, ncol = length(parameter))
            rownames(x) <- c("Ter", "a", "P", "sda", "rd", "z", "s_Ter", "s_z")
            colnames(x) <- parameter
            return(x)
          })
          names(mm) <- conditions
          dm <- matrix(nrow = 2, ncol = length(parameter))
          rownames(dm) <- c("Upper_Limit", "Lower_Limit")
          colnames(dm) <- parameter
          mm$Cong["Ter", ] <- c(1, 0, 0, 0, 0)
          mm$Cong["a", ] <- c(0, 1, 0, 0, 0)
          mm$Cong["P", ] <- c(0, 0, 1, 0, 0)
          mm$Cong["sda", ] <- c(0, 0, 0, 1, 0)
          mm$Cong["rd", ] <- c(0, 0, 0, 0, 1)
          mm$Cong["z", ] <- c(0, 0, 0, 0, 0)
          mm$Cong["s_Ter", ] <- c(0, 0, 0, 0, 0)
          mm$Cong["s_z", ] <- c(0, 0, 0, 0, 0)
          mm$Incong["Ter", ] <- c(1, 0, 0, 0, 0)
          mm$Incong["a", ] <- c(0, 1, 0, 0, 0)
          mm$Incong["P", ] <- c(0, 0, -1, 0, 0)
          mm$Incong["sda", ] <- c(0, 0, 0, 1, 0)
          mm$Incong["rd", ] <- c(0, 0, 0, 0, 1)
          mm$Incong["z", ] <- c(0, 0, 0, 0, 0)
          mm$Incong["s_Ter", ] <- c(0, 0, 0, 0, 0)
          mm$Incong["s_z", ] <- c(0, 0, 0, 0, 0)
          dm["Upper_Limit", ] <- c(0.45, 0.19, 0.55, 2.6, 0.026)
          dm["Lower_Limit", ] <- c(0.15, 0.07, 0.2, 1, 0.01)
        }
      )
      sp <- as.matrix(data.frame(dt = dt, sigma = sigma))
      rf <- list(CDF = CDF_perc, CAF = CAF_perc)
      return(methods::new("DDModel", ID = model, MM = mm, DM = dm, SP = sp, RF = rf))
    },
    RMT_LDT = {
      switch(EXPR = model,
        DDM_classic = {
          conditions <- c("New_noWord", "Old_Word")
          parameter <- c("Ter", "a", "mu", "z", "s_Ter", "s_mu", "s_z")
          dt <- 0.001
          sigma <- 0.1
          mm <- lapply(1:length(conditions), function(x) {
            x <- matrix(0, nrow = 7, ncol = length(parameter))
            rownames(x) <- c("Ter", "a", "mu", "z", "s_Ter", "s_mu", "s_z")
            colnames(x) <- parameter
            return(x)
          })
          names(mm) <- conditions
          dm <- matrix(nrow = 2, ncol = length(parameter))
          rownames(dm) <- c("Upper_Limit", "Lower_Limit")
          colnames(dm) <- parameter
          mm$New_noWord["Ter", ] <- c(1, 0, 0, 0, 0, 0, 0)
          mm$New_noWord["a", ] <- c(0, 1, 0, 0, 0, 0, 0)
          mm$New_noWord["mu", ] <- c(0, 0, 1, 0, 0, 0, 0)
          mm$New_noWord["z", ] <- c(0, 0, 0, 1, 0, 0, 0)
          mm$New_noWord["s_Ter", ] <- c(0, 0, 0, 0, 1, 0, 0)
          mm$New_noWord["s_mu", ] <- c(0, 0, 0, 0, 0, 1, 0)
          mm$New_noWord["s_z", ] <- c(0, 0, 0, 0, 0, 0, 1)
          mm$Old_Word["Ter", ] <- c(1, 0, 0, 0, 0, 0, 0)
          mm$Old_Word["a", ] <- c(0, 1, 0, 0, 0, 0, 0)
          mm$Old_Word["mu", ] <- c(0, 0, 1, 0, 0, 0, 0)
          mm$Old_Word["z", ] <- c(0, 0, 0, -1, 0, 0, 0)
          mm$Old_Word["s_Ter", ] <- c(0, 0, 0, 0, 1, 0, 0)
          mm$Old_Word["s_mu", ] <- c(0, 0, 0, 0, 0, 1, 0)
          mm$Old_Word["s_z", ] <- c(0, 0, 0, 0, 0, 0, 1)
          dm["Upper_Limit", ] <- c(0.5, 0.2, 0.4, 0.02, 0.2, 0.1, 0.05)
          dm["Lower_Limit", ] <- c(0.2, 0.05, 0.0, -0.02, 0.0, 0.0, 0.0)
        }
      )
      sp <- as.matrix(data.frame(dt = dt, sigma = sigma))
      rf <- list(CDF = CDF_perc, CAF = CAF_perc)
      return(methods::new("DDModel", ID = model, MM = mm, DM = dm, SP = sp, RF = rf))
    }
  )
}

#' @rdname show-methods
#' @aliases show,DDModel-method
#' @exportMethod show
setMethod("show", "DDModel", function(object) {
  cat(
    "DDModel Object: \n",
    "\nModel: ", object@ID, "\n",
    "\nModelMatrix: \n"
  )
  print(object@MM)
  cat("Parameter Domain: \n")
  print(object@DM)
  cat("\nSimulation Parameter: \n")
  print(object@SP)
  cat("\nForm of Representation: \n")
  print(object@RF)
})
