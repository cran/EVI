#' Function that calculates the standard deviation
#'
#' Calculates and returns the standard deviation of a vector
#' @param x numeric vector
#'
#'
#' @examples
#' x=c(2,2,1,4,9,10,23,10,9,10,14,12,10)
#' medvol(x)
#'
#' @export
medvol=function(x){
  sdx<-sd(x)
  return(sdx)
}


