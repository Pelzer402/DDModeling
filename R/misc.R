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
