#' Function to generate a DDRep from a RAW data format
#'
#' @name Get_DDRep
#'
#' @param data \code{string} or \code{list} of strings to a document containing the RAW data
#'
#' @return \code{list} of \coe{DDRep} object(s)
Get_DDRep <- function(model=NULL,data = NULL){
  Flag <- 1
  Check <- ArgumentCheck::newArgCheck()
  if (is.null(model))
  {
    ArgumentCheck::addError(msg = "Data is missing!",argcheck = Check)
    Flag <- 99
  }
  if (is.data.frame(data))
  {
    ArgumentCheck::addError(msg = "Data needs to be entered in a list()",argcheck = Check)
    Flag <- 99
  }
  if (Flag == 99)
  {
    return(cat("Get_DDREp failed"))
  }
  else
  {
    trans <- list()
    for (i in 1:length(data))
    {
      buff <- split(x = data[[i]],f = data[[i]]$cond)
      trans[[i]] <- .Get_DDRep(model,buff)
    }
    return(trans)
  }
}
