#' Calculation of the Epidemic Volatility Index
#'
#' Calculates the relative change in the standard deviation between two consecutive rolling windows.
#'
#' @param rollsd numeric vector - returned and stored as roll from the rollsd() function.
#'
#'
#' @examples
#' data("Italy")
#' cases = mova(cases=Italy$Cases)
#' roll = rollsd(cases=cases, lag_t=7)
#' ev = evi(rollsd=roll)
#'
#' @export

evi = function(rollsd) {
  evi=rep(NA,length(rollsd))
  for (i in 2:length(rollsd)){
    evi[i]= (rollsd[i]-rollsd[i-1])/rollsd[i-1]
  }
  evi[is.nan(evi)]<-0
  return(evi)
}



