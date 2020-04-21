#' \code{DDFit} class definition
#' @include DDRep.R DDModel.R DDFitPar.R
#' @name DDFit-class
#' @rdname DDFit-class
#' @slot INP_REP    \code{DDRep} object containing the to be fitted data
#' @slot FIT_REP    \code{DDRep} object containing the fitted data
#' @slot MODEL      \code{DDModel} object containing the model that was used in the fit
#' @slot FIT        \code{DDFitPar} object containing information regarding the fit
setClass("DDFit",
         slots      = list(
           INP_REP      = "DDRep",
           FIT_REP      = "DDRep",
           MODEL        = "DDModel",
           FIT          = "DDFitPar"
         )
)

#' @rdname DDFit-class
#' @aliases summary,DDFit-method
#' @exportMethod summary
setMethod("summary",signature("DDFit"),function(object){
  if (all(object@INP_REP@PAR == 0))
  {
   eta_buff <- "Population parameters not known! No calculation possible."
  }
  else
  {
    ranges <- object@MODEL@DM
    ranges <- ranges["Upper_Limit",] - ranges["Lower_Limit",]
    eta_buff    <- abs(object@INP_REP@PAR - object@FIT_REP@PAR)/ranges
  }
  OUT <- list(INPUT_Par = object@INP_REP@PAR, FIT_Par = object@FIT_REP@PAR, Eta = eta_buff, mean_Eta = mean(as.numeric(eta_buff)), Fit = object@FIT@FIT_V)
  print(OUT)
})

#' @rdname DDFit-class
#' @aliases plot,DDFit-method
#' @exportMethod plot
setMethod("plot",signature("DDFit"),function(x){
  perc <- type <- cond <- acc <- NULL
  CAF <- c()
  CDF <- c()
  for ( i in 1:length(x@MODEL@MM))
  {
    CAF <- rbind(CAF,cbind(x@INP_REP@REP$CAF[[i]],type="TBF"),cbind(x@FIT_REP@REP$CAF[[i]],type="FIT"))
    CDF <- rbind(CDF,cbind(x@INP_REP@REP$CDF[[i]],type="TBF"),cbind(x@FIT_REP@REP$CDF[[i]],type="FIT"))
  }
  CAF_time <- ggplot2::ggplot(CAF,ggplot2::aes(x=perc,y=time,color=type,linetype=cond,shape=type)) +
          ggplot2::geom_line() +
          ggplot2::geom_point() +
          ggplot2::labs(title = "CAF time",x = "percentile",y = "reaction time [ms]", linetype = "Condition", type = "Condition", color = "Data" ) +
          ggplot2::theme(legend.position = "bottom")
  L <- cowplot::get_legend(CAF_time)
  CAF_time <- CAF_time + ggplot2::theme(legend.position = "none")
  CAF_acc <- ggplot2::ggplot(CAF,ggplot2::aes(x=perc,y=acc,color=type,linetype=cond,shape=type)) +
          ggplot2::geom_line() +
          ggplot2::geom_point() +
          ggplot2::labs(title = "CAF acc",x = "percentile",y = "accuracy") +
          ggplot2::theme(legend.position = "none")
  CDF <- ggplot2::ggplot(CDF,ggplot2::aes(x=perc,y=time,color=type,linetype=cond,shape=type)) +
          ggplot2::geom_line() +
          ggplot2::geom_point() +
          ggplot2::labs(title = "CDF time",x = "percentile",y = "reaction time [ms]") +
          ggplot2::theme(legend.position = "none")
  C1 <- cowplot::plot_grid(CAF_time,CAF_acc,CDF,nrow = 1,ncol = 3,rel_heights = c(1,1,1))
  cat("Fit: ")
  cat(as.numeric(x@FIT@FIT_V))
  cowplot::plot_grid(C1,L,nrow=2,ncol=1,rel_heights=c(9,1))
})


