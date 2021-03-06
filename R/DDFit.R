#' \code{DDFit} class definition
#' @include DDRep.R DDModel.R DDFitPar.R
#' @name DDFit-class
#' @rdname DDFit-class
#' @slot INP_REP    \code{DDRep} object containing the to be fitted data
#' @slot FIT_REP    \code{DDRep} object containing the fitted data
#' @slot MODEL      \code{DDModel} object containing the model that was used in the fit
#' @slot FIT        \code{DDFitPar} object containing information regarding the fit
setClass("DDFit",
  slots = list(
    INP_REP = "DDRep",
    FIT_REP = "DDRep",
    MODEL = "DDModel",
    FIT = "DDFitPar"
  )
)

#' @rdname summary-methods
#' @aliases summary
#' @exportMethod summary
setMethod("summary", signature("DDFit"), function(object) {
  if (all(object@INP_REP@PAR == 0)) {
    eta_buff <- "Population parameters not known! No calculation possible."
  }
  else {
    ranges <- object@MODEL@DM
    ranges <- ranges["Upper_Limit", ] - ranges["Lower_Limit", ]
    eta_buff <- abs(object@INP_REP@PAR - object@FIT_REP@PAR) / ranges
  }
  OUT <- list(INPUT_Par = object@INP_REP@PAR, FIT_Par = object@FIT_REP@PAR, Eta = eta_buff, mean_Eta = mean(as.numeric(eta_buff)), Fit = object@FIT@FIT_V)
  print(OUT)
})

#' @rdname plot-methods
#' @aliases plot
#' @exportMethod plot
setMethod("plot", signature("DDFit"), function(x,type="simple") {
  switch(type,
         simple = {
           perc <- type <- cond <- acc <- NULL
           CAF_F <- rbindlist(x@FIT_REP@REP$CAF)
           CAF_I <- rbindlist(x@INP_REP@REP$CAF)
           CDF_F <- rbindlist(x@FIT_REP@REP$CDF)
           CDF_I <- rbindlist(x@INP_REP@REP$CDF)
           CDF_F$Source <- "Fit"
           CDF_I$Source <- "Observation"
           CAF_F$Source <- "FIT"
           CAF_I$Source <- "Observation"
           CAF_time <- ggplot2::ggplot() +
             ggplot2::geom_line(data = CAF_F,mapping = ggplot2::aes(x = perc, y = time,color=cond)) +
             ggplot2::geom_point(data = CAF_F,mapping = ggplot2::aes(x = perc, y = time,color=cond, shape = Source)) +
             ggplot2::geom_point(data = CAF_I,mapping = ggplot2::aes(x = perc, y = time,color=cond, shape = Source)) +
             ggplot2::labs(title = "CAF time", x = "percentile", y = "reaction time [ms]", color = "Condition", shape = "Data") +
             ggplot2::theme(legend.position = "bottom")
           L <- cowplot::get_legend(CAF_time)
           CAF_time <- CAF_time + ggplot2::theme(legend.position = "none")
           CAF_acc <- ggplot2::ggplot() +
             ggplot2::geom_line(data = CAF_F,mapping = ggplot2::aes(x = perc, y = acc,color=cond)) +
             ggplot2::geom_point(data = CAF_F,mapping = ggplot2::aes(x = perc, y = acc,color=cond, shape = Source)) +
             ggplot2::geom_point(data = CAF_I,mapping = ggplot2::aes(x = perc, y = acc,color=cond, shape = Source)) +
             ggplot2::labs(title = "CAF acc", x = "percentile", y = "accuracy") +
             ggplot2::theme(legend.position = "none")
           CDF <- ggplot2::ggplot() +
             ggplot2::geom_line(data = CDF_F,mapping = ggplot2::aes(x = perc, y = time,color=cond)) +
             ggplot2::geom_point(data = CDF_F,mapping = ggplot2::aes(x = perc, y = time,color=cond, shape = Source)) +
             ggplot2::geom_point(data = CDF_I,mapping = ggplot2::aes(x = perc, y = time,color=cond, shape = Source)) +
             ggplot2::labs(title = "CDF time", x = "percentile", y = "reaction time [ms]") +
             ggplot2::theme(legend.position = "none")
           C1 <- cowplot::plot_grid(CAF_time, CAF_acc, CDF, nrow = 1, ncol = 3, rel_heights = c(1, 1, 1))
           cat("Fit: ")
           cat(as.numeric(x@FIT@FIT_V))
           cowplot::plot_grid(C1, L, nrow = 2, ncol = 1, rel_heights = c(9, 1))
         },
         research = {
           perc <- type <- cond <- acc <- NULL
           CAF_F <- rbindlist(x@FIT_REP@REP$CAF)
           CAF_I <- rbindlist(x@INP_REP@REP$CAF)
           CDF_F <- rbindlist(x@FIT_REP@REP$CDF)
           CDF_I <- rbindlist(x@INP_REP@REP$CDF)
           CDF_F$Source <- "Fit"
           CDF_I$Source <- "Observation"
           CAF_F$Source <- "FIT"
           CAF_I$Source <- "Observation"
           CAF_P <- ggplot2::ggplot() +
             ggplot2::geom_line(data = CAF_F,mapping = ggplot2::aes(x = time, y = acc,color=cond)) +
             ggplot2::geom_point(data = CAF_F,mapping = ggplot2::aes(x = time, y = acc,color=cond, shape = Source)) +
             ggplot2::geom_point(data = CAF_I,mapping = ggplot2::aes(x = time, y = acc,color=cond, shape = Source)) +
             ggplot2::labs(title = "CAF time", x = "reaction time [ms]", y = "accuracy", color = "Condition", shape = "Data") +
             ggplot2::theme(legend.position = "bottom")
           L <- cowplot::get_legend(CAF_P)
           CAF_P <- CAF_P + ggplot2::theme(legend.position = "none")
           CDF_P <- ggplot2::ggplot() +
             ggplot2::geom_line(data = CDF_F,mapping = ggplot2::aes(x = perc, y = time,color=cond)) +
             ggplot2::geom_point(data = CDF_F,mapping = ggplot2::aes(x = perc, y = time,color=cond, shape = Source)) +
             ggplot2::geom_point(data = CDF_I,mapping = ggplot2::aes(x = perc, y = time,color=cond, shape = Source)) +
             ggplot2::labs(title = "CDF time", x = "percentile", y = "reaction time [ms]") +
             ggplot2::theme(legend.position = "none")
           C1 <- cowplot::plot_grid(CAF_P, CDF_P, nrow = 1, ncol = 2, rel_heights = c(1, 1))
           cat("Fit: ")
           cat(as.numeric(x@FIT@FIT_V))
           cowplot::plot_grid(C1, L, nrow = 2, ncol = 1, rel_heights = c(9, 1))
         })
})
