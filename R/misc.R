# Internal rescaling function
unscale_scale <- function(IN = NULL,unscaler = NULL){
  if (is.null(unscaler$Center) & is.null(unscaler$Stddev))
  {
    return(IN)
  }
  else
  {
    return(t(t(IN)*unscaler$Stddev + unscaler$Center))
  }
}
