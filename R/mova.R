#' Function that calculates the moving average of a time series
#'
#' Calculates the moving average for a time series.
#'
#' @param cases the time series of the newly observed cases per unit of time (ideally per day).
#' @param r_a The window size for the moving average that will be analyzed. If set to 1 the actual observations are analyzed. However, due to the variability of the reported cases between working days and weekends it is recommended that the 7-day moving average is analyzed (i.e. r_a = 7), which is the default for this argument. Users could prefer a longer interval of 14 days or one month (i.e., r_a=14 or 30, respectively).
#'
#' @examples
#' data("Italy")
#' mova(cases=Italy$Cases, r_a=7)
#' mova(cases=Italy$Cases, r_a=14)
#' @export
#'
#'
mova=function(cases, r_a=7){
  ncases=rep(NA, length(cases))
  for (i in 1:length(cases)){
    ncases[i]=mean(cases[((i+1)-min(r_a,i)):i])
  }
  return(ncases)
}


