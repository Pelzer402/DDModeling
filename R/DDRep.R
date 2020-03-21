#' \code{DDRep} class definition
#' @name DDRep-class
#' @rdname DDRep-class
#' @slot RAW  \code{list} of data.frames that contain the RAW data (i.e. 3 coloumns: $cond $resp $time)
#' @slot REP  \code{list} of data.frames that contain data representations (CDF and CAF)
#' @slot RF   \code{list} of numeric vectors that contain the percentiles of the representation.
#' @slot PAR  \code{data.frame} containing the parameters associated with the representation
setClass("DDRep",
         slots      = list(
           RAW      = "list",
           REP      = "list",
           RF       = "list",
           PAR      = "data.frame"
         )
)

#' Constructor for \link{DDRep-class}
#' @include DDRep.R
#' @name DDRep
#' @description Userfriendly function for the construction of a \link{DDRep-class}.
#' @param model  \code{\link{DDModel-class}}.
#' @param raw    \code{data.frame} that contains the raw data in a specific format (see details for further information)
#' @return \code{\link{DDRep-class}}.
#' @details Raw data in the DDModeling package is handled in \code{data.frames} containing three coulmns:
#' \describe{
#'   \item{cond}{Factor: Specifying the condition under which a given measure was taken.}
#'   \item{resp}{Numeric: Specifying the binary coded response property (0==false, 1==correct)}
#'   \item{time}{Numeric: Specifying the time (in rounded ms) for a given measure.}
#' }
#' Note: In order for the DDRep function to execute successfully the factor levels of the cond coulmns have to be identical to the condition names specified
#' in the MM slot of the \code{model}!
#' @export
DDRep <- function(model = NULL,raw=NULL){
  Flag <- 1
  Check <- ArgumentCheck::newArgCheck()
  if (!methods::is(model,"DDModel"))
  {
    ArgumentCheck::addError(msg = "'model' is missing or in the wrong format!",argcheck = Check)
    Flag <- 99
  }
  if (!is.data.frame(raw))
  {
    ArgumentCheck::addError(msg = "'raw' is missing or in the wrong format!",argcheck = Check)
    Flag <- 99
  }
  if (Flag == 1)
  {
    if ((length(colnames(raw)) == 3) && all(colnames(raw) %in% c("cond","resp","time")))
    {
      if (is.factor(raw$cond) && is.numeric(raw$resp) && is.numeric(raw$time))
      {
        if ((length(levels(raw$cond)) == length(names(model@MM))) && (all(levels(raw$cond) %in% names(model@MM))))
        {
          raw$cond <- factor(raw$cond,levels = names(model@MM),labels = names(model@MM))
          raw <- split(raw,raw$cond)
          names(raw) <- names(model@MM)
        }
        else
        {
          ArgumentCheck::addError(msg = "The condition names in 'raw' do not match with the one specified in the 'model'!",argcheck = Check)
          Flag <- 99
        }
      }
      else
      {
        ArgumentCheck::addError(msg = "The data formats of the coulmns in 'raw' do not qualify for this computation! See ?DDRep details section!",argcheck = Check)
        Flag <- 99
      }
    }
    else
    {
      ArgumentCheck::addError(msg = "The coulmnnames in 'raw' do not qualify for this computation! See ?DDRep details section!",argcheck = Check)
      Flag <- 99
    }
  }
  ArgumentCheck::finishArgCheck(Check)
  if (Flag == 1)
  {
    return(.DDRep_cpp(model,raw))
  }
  else
  {
    return(cat("DDRep failed"))
  }
}

# Internal converter function
DDRep_wide <- function(rep = NULL){
  OUT <- c()
  for(c in 1:length(rep@RAW))
  {
    OUT <- c(OUT,rep@REP$CDF[[c]]$time)
  }
  for(c in 1:length(rep@RAW))
  {
    for(p in  1:length(rep@REP$CAF[[c]]$time))
    {
      OUT <- c(OUT,rep@REP$CAF[[c]]$time[p])
      OUT <- c(OUT,rep@REP$CAF[[c]]$acc[p])
    }
  }
  return(OUT)
}

#' @export
setMethod("show","DDRep",function(object){
  cat("CDF: \n")
  print(object@REP$CDF)
  cat("\nCAF: \n")
  print(object@REP$CAF)
  cat("\nParameter: \n")
  print(object@PAR)
})

#' @export
setMethod("plot",signature(x="DDRep"),function(x){
  CAF <- c()
  CDF <- c()
  for ( i in 1:length(x@RAW))
  {
    CAF <- rbind(CAF,x@REP$CAF[[i]])
    CDF <- rbind(CDF,x@REP$CDF[[i]])
  }
  CAF_time <- ggplot2::ggplot(CAF,ggplot2::aes(x=perc,y=time,linetype=cond)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(title = "CAF time",x = "percentile",y = "reaction time [ms]", linetype = "Condition") +
    ggplot2::theme(legend.position = "bottom")
  L <- cowplot::get_legend(CAF_time)
  CAF_time <- CAF_time + ggplot2::theme(legend.position = "none")
  CAF_acc <- ggplot2::ggplot(CAF,ggplot2::aes(x=perc,y=acc,linetype=cond)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(title = "CAF acc",x = "percentile",y = "accuracy") +
    ggplot2::theme(legend.position = "none")
  CDF <- ggplot2::ggplot(CDF,ggplot2::aes(x=perc,y=time,linetype=cond)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(title = "CDF time",x = "percentile",y = "reaction time [ms]") +
    ggplot2::theme(legend.position = "none")
  C1 <- cowplot::plot_grid(CAF_time,CAF_acc,CDF,nrow = 1,ncol = 3,rel_heights = c(1,1,1))
  cowplot::plot_grid(C1,L,nrow=2,ncol=1,rel_heights=c(9,1))
})

