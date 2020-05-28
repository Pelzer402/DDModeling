# Internal rescaling function
unscale_scale <- function(IN = NULL,unscaler = NULL){
  if ((unscaler$Center == FALSE) || (unscaler$Stddev == FALSE))
  {
    return(IN)
  }
  else
  {
    return(t(t(IN)*unscaler$Stddev + unscaler$Center))
  }
}

#' Retrieves scale Information from your data!
#'
#' @name Scale_DL_Data
#' @rdname Scale_DL_Data
#' @param  data \code{List} containing INPUT and OUTPUT element (as returned from Import_GRID)
#' @return DL_scale (\code{List}) for the given data
Scale_DL_Data <- function(data = NULL)
{
  data$INPUT <- scale(data$INPUT)
  data$OUTPUT <- scale(data$OUTPUT)
  out <- list(
    INPUT = list(
      Center = attr(data$INPUT,"scaled:center"),
      Stddev = attr(data$INPUT,"scaled:scale")
    ),
    OUTPUT = list(
      Center = attr(data$OUTPUT,"scaled:center"),
      Stddev = attr(data$OUTPUT,"scaled:scale")
    )
  )
  return(out)
}
