#' Function to fit a given DDRep to a given DDModel
#' @name Fit_DDModel
#' @param model \code{DDModel} object
#' @param data   \code{DDRep} object or list of \code{DDRep} objects
#' @param DL_model (optional) \code{Model} in the form of a keras neural network model
#' @param DL_scale (optional) \code{list} containing mean and sd of the transformation used in the deep learning model while training.
#' @param grid_path  (optional) \code{path} to a directory containing a .GRID fileset. If NULL the model will be fitted using 20 randomly drawn startparametersets from the model-DOMAIN.
#' @param trials \code{integer} indicating the number of trials used while fitting (s_sampling = FALSE) or the maximum number of trials used while super sampling (s_sampling = TRUE)
#' @param s_sampling \code{bool} indicating super sampling while fitting
#' @param simplex_struc \code{numeric vector} containing the number of simplex iterations per sorting cycle.
#' @param simplex_coef \code{list} containing named parameters for the SIMPLEX
#' @return \code{DDFit} object
#' @details This is a rather complex function that is very flexible depending on your input!
#' \describe{
#'   \item{Random Fitting:}{Exclude all optional parameters}
#'   \item{Grid Fitting:}{Specify a grid_path}
#'   \item{Deep Learning Fitting}{Specify DL_model and DL_scale}
#' }
#' @export
Fit_DDModel <- function(model = NULL, data = NULL,
                        DL_model = NULL, DL_scale = list(INPUT = list(Center = FALSE, Stddev = FALSE), OUTPUT = list(Center = FALSE, Stddev = FALSE)),
                        grid_path = NULL,
                        trials = 10000L, s_sampling = FALSE,
                        simplex_struc = c(0L), simplex_coef = list(alpha = 1.0, beta = 0.5, gamma = 2.0, sigma = 0.5, tol = 0.0000001, nfunc = 120L, nshrink = 3L)) {
  args <- c("model", "data", "DL_model", "DL_scale", "grid_path", "trials", "s_sampling", "simplex_struc", "simplex_coef")
  formats <- list("DDModel", c("DDRep", "list"), NULL, "list", "character", "integer", "logical", "integer", "list")
  codes <- list(
    RND_single  = c(1, 1, 0, 1, 0, 1, 1, 1, 1),
    RND_multi   = c(1, 2, 0, 1, 0, 1, 1, 1, 1),
    GRID_single = c(1, 1, 0, 1, 1, 1, 1, 1, 1),
    GRID_multi  = c(1, 2, 0, 1, 1, 1, 1, 1, 1),
    DL_single   = c(1, 1, 1, 1, 0, 1, 1, 1, 1),
    DL_multi    = c(1, 2, 1, 1, 0, 1, 1, 1, 1)
  )
  names(formats) <- args
  format_checks <- c(rep(NULL, length(args)))
  Check <- ArgumentCheck::newArgCheck()
  for (i in 1:length(args)) {
    if (args[i] == "DL_model") {
      if (is.null(get(args[i]))) {
        format_checks[i] <- 0
      }
      else {
        format_checks[i] <- 1
      }
    }
    else {
      for (f in 1:length(formats[[i]])) {
        if (methods::is(get(args[i]), formats[[i]][f])) {
          format_checks[i] <- f
          switch(i,
            "model" = {
            },
            "data" = {
              if (format_checks[1] == 1) {
                if (format_checks[i] == 1) {
                  if (!identical(data@RF, model@RF)) {
                    ArgumentCheck::addError(msg = "The representation (@RF) of 'data' does not conform to the one specified in 'model'! => Please make sure that both match eachother!", argcheck = Check)
                  }
                }
                else if (format_checks[i] == 2) {
                  error_rep <- c()
                  error_RF <- c()
                  il <- 1
                  while (il <= length(data)) {
                    if (!methods::is(data[[il]], "DDRep")) {
                      error_rep <- c(error_rep, il)
                    }
                    else {
                      if (!identical(data[[il]]@RF, model@RF)) {
                        error_RF <- c(error_RF, il)
                      }
                    }
                    il <- il + 1
                  }
                  if (length(error_rep) > 0) {
                    ArgumentCheck::addError(msg = paste0("The type of your 'data' at indices '", paste(error_rep, collapse = ", "), "' does not match a 'DDRep'! => Please make sure that all your data is resembled as a DDRep!"), argcheck = Check)
                  }
                  if (length(error_RF) > 0) {
                    ArgumentCheck::addError(msg = paste0("The representation (@RF) of your 'data' at indices '", paste(error_RF, collapse = ", "), "' does not conform to the one specified in 'model'! => Please make sure that both match eachother!"), argcheck = Check)
                  }
                  if (length(data) == 0) {
                    ArgumentCheck::addError(msg = "'data' is missing or in the wrong format!", argcheck = Check)
                  }
                }
              }
            },
            "DL_model" = {
            },
            "DL_scale" = {
              if (!names(DL_scale) %in% c("INPUT", "OUTPUT") ||
                !length(names(DL_scale)) == 2 ||
                !names(DL_scale$INPUT) %in% c("Center", "Stddev") ||
                !length(names(DL_scale$INPUT)) == 2 ||
                !names(DL_scale$OUTPUT) %in% c("Center", "Stddev") ||
                !length(names(DL_scale$OUTPUT)) == 2) {
                ArgumentCheck::addError(msg = "'DL_scale' is in the wrong format! Please check the 'Fit_DDModel' help page for more information!", argcheck = Check)
              }
            },
            "grid_path" = {
              if (format_checks[1] == 1) {
                Grid_model <- readRDS(list.files(grid_path, full.names = TRUE, pattern = "\\.Gcfg$"))
                if (!identical(Grid_model, model)) {
                  ArgumentCheck::addError(msg = "The model used for your GRID at 'grid_path' does not conform to the one specified in 'model'! => Please make sure that both match eachother! You may want to import the one in your 'grid_path'.", argcheck = Check)
                }
              }
            },
            "trials" = {
            },
            "s_sampling" = {
            },
            "simplex_struc" = {
              if (!all(simplex_struc >= 0)) {
                ArgumentCheck::addError(msg = "All values specified in 'simplex_struc' must be >0!", argcheck = Check)
              }
              if (length(simplex_struc) > 1) {
                for (it in 1:(length(simplex_struc) - 1))
                {
                  if (simplex_struc[it + 1] > simplex_struc[it]) {
                    ArgumentCheck::addError(msg = "All values specified in 'simplex_struc' must be in an at least descending order!", argcheck = Check)
                    break
                  }
                }
              }
            },
            "simplex_coef" = {
              n_coef <- c("alpha", "beta", "gamma", "sigma", "tol", "nfunc", "nshrink")
              def_coef <- list(alpha = 1.0, beta = 0.5, gamma = 2.0, sigma = 0.5, tol = 0.0000001, nfunc = 120L, nshrink = 3L)
              mod <- n_coef %in% names(simplex_coef)
              for (i in 1:length(mod))
              {
                if (mod[i]) {
                  def_coef[[n_coef[i]]] <- simplex_coef[[n_coef[i]]]
                }
              }
              if (!all(mod)) {
                ArgumentCheck::addMessage(msg = paste0("You did non specify the following SIMPLEX parameter in 'simplex_coef': ", paste(n_coef[!mod], collapse = ", "), ". Default values will be assumed!"), argcheck = Check)
              }
              mod <- names(simplex_coef) %in% n_coef
              if (!all(mod)) {
                ArgumentCheck::addMessage(msg = paste0("The values: ", paste(simplex_coef[!mod], collapse = ", "), " specified in 'simplex_coef' have no meaningfull intepretation for the function and therefore will be discarded!"), argcheck = Check)
              }
              simplex_coef <- def_coef
            }
          )
          break
        }
        else {
          if (f == length(formats[[i]])) {
            switch(i,
              "model" = {
                ArgumentCheck::addError(msg = "'model' is missing or in the wrong format!", argcheck = Check)
              },
              "data" = {
                ArgumentCheck::addError(msg = "'data' is missing or in the wrong format!", argcheck = Check)
              },
              "DL_model" = {
              },
              "DL_scale" = {
              },
              "grid_path" = {
              },
              "trials" = {
                ArgumentCheck::addError(msg = "'trials' is missing or in the wrong format!", argcheck = Check)
              },
              "s_sampling" = {
                ArgumentCheck::addError(msg = "'s_sampling' is missing or in the wrong format!", argcheck = Check)
              },
              "simplex_struc" = {
                ArgumentCheck::addError(msg = "'simplex_struc' is missing or in the wrong format!", argcheck = Check)
              },
              "simplex_coef" = {
                ArgumentCheck::addError(msg = "'simplex_coef' is missing or in the wrong format!", argcheck = Check)
              }
            )
            format_checks[i] <- 0
          }
        }
      }
    }
  }
  ArgumentCheck::finishArgCheck(Check)
  routine <- 99
  for (i in 1:length(codes))
  {
    if (all(codes[[i]] == format_checks)) {
      routine <- i
    }
  }
  switch(routine,
    "RND_singel" = {
      COMP_List <- list(model, data, s_sampling, trials, simplex_struc, simplex_coef)
      return(.Fit_DDModel_rnd(COMP_List))
    },
    "RND_multi" = {
      ncores <- parallel::detectCores() - 1
      clust <- ParallelLogger::makeCluster(ncores)
      COMP_List <- list()
      for (i in 1:length(data))
      {
        COMP_List[[i]] <- list(model, data[[i]], s_sampling, trials, simplex_struc, simplex_coef)
      }
      RES <- ParallelLogger::clusterApply(clust, COMP_List, .Fit_DDModel_rnd)
      ParallelLogger::stopCluster(clust)
      return(RES)
    },
    "GRID_singel" = {
      Grid_parts <- list.files(grid_path, full.names = TRUE, pattern = "\\.GRID$")
      COMP_List <- list(model, data, Grid_parts, s_sampling, trials, grid_path, simplex_struc, simplex_coef)
      return(.Fit_DDModel_grid(COMP_List))
    },
    "GRID_multi" = {
      Grid_parts <- list.files(grid_path, full.names = TRUE, pattern = "\\.GRID$")
      ncores <- parallel::detectCores() - 1
      clust <- ParallelLogger::makeCluster(ncores)
      COMP_List <- list()
      for (i in 1:length(data))
      {
        COMP_List[[i]] <- list(model, data[[i]], Grid_parts, s_sampling, trials, grid_path, simplex_struc, simplex_coef)
      }
      RES <- ParallelLogger::clusterApply(clust, COMP_List, .Fit_DDModel_grid)
      ParallelLogger::stopCluster(clust)
      return(RES)
    },
    "DL_singel" = {
      COMP_List <- list()
      if (is.null(DL_scale)) {
        INP <- t(DDRep_wide(data))
        PRE <- as.numeric(predict(DL_model, INP))
        COMP_List <- list(model, data, s_sampling, trials, simplex_struc, PRE, simplex_coef)
      }
      else {
        INP <- t(DDRep_wide(data))
        INP <- scale(INP, DL_scale$INPUT$Center, DL_scale$INPUT$Stddev)
        PRE <- as.numeric(predict(DL_model, INP))
        PRE <- unscale_scale(IN = PRE, unscaler = DL_scale$OUTPUT)
        COMP_List <- list(model, data, s_sampling, trials, simplex_struc, PRE, simplex_coef)
      }
      return(.Fit_DDModel_DL(COMP_List))
    },
    "DL_multi" = {
      ncores <- parallel::detectCores() - 1
      clust <- ParallelLogger::makeCluster(ncores)
      COMP_List <- list()
      for (i in 1:length(data))
      {
        if (is.null(DL_scale)) {
          INP <- t(DDRep_wide(data[[i]]))
          PRE <- as.numeric(predict(DL_model, INP))
          COMP_List[[i]] <- list(model, data[[i]], s_sampling, trials, simplex_struc, PRE, simplex_coef)
        }
        else {
          INP <- t(DDRep_wide(data[[i]]))
          INP <- scale(INP, DL_scale$INPUT$Center, DL_scale$INPUT$Stddev)
          PRE <- as.numeric(predict(DL_model, INP))
          PRE <- unscale_scale(IN = PRE, unscaler = DL_scale$OUTPUT)
          COMP_List[[i]] <- list(model, data[[i]], s_sampling, trials, simplex_struc, PRE, simplex_coef)
        }
      }
      RES <- ParallelLogger::clusterApply(clust, COMP_List, .Fit_DDModel_DL)
      ParallelLogger::stopCluster(clust)
      return(RES)
    }
  )
}
