#' Calculation of the convergence Epidemic Volatility Index
#'
#' This sencondary function produces the convergence Epidemic Volatility Index based on input data.
#'
#' @param cases the time series of the newly observed cases per unit of time (ideally per day).
#' @param lag_n Integer. Restriction of the window size for the rolling window size.
#' @param c_n threshold alpha-level value (0 <= c <= 0.5) for issuing an early warning. If cevi <= c_n an early warning is issued and otherwise is not.
#'
#' For each time point the stored variables are:
#' @return
#' \itemize{
#' \item{Dates: the date for each time point (with origin 01-01-1970).}
#'
#' \item{Days: the serial number for each time point.}
#'
#' \item{EVI: the estimated EVI for each time point.}}
#' @examples
#' cEVI_fun(cases=c(0,0,1,3,4,10,40,90,105,160,210,301,510,670,680,650,670,665),lag_n=3,c_n=0.1)
#' @export
#'
#' @references
#' Pateras K, Meletis E, Denwood M, et al. The convergence epidemic index (cEVI) an early warning tool for identifying waves in an epidemic. Inf Dis Mod, (2023)

cEVI_fun<-function(cases,lag_n,c_n){

  cevi <- rep(NA, length(cases))
  for(k in (lag_n-1):(length(cases)-(lag_n+1))){
    enu=mean(cases[(k+2):(k+lag_n+1)]-cases[(k+1):(k-(lag_n-2))],na.rm = T)
    den1=sd(cases[(k+1):(k-(lag_n-2))])^2/(length(cases[(k+1):(k-(lag_n-2))]))
    den2=sd(cases[(k+2):(k+lag_n+1)])^2/(length(cases[(k):(k+lag_n+1)]))

    # Spectral variances possibly more appropriate but more time consuming
    #den1=spectrum0.ar(cases[(i+1):(i+w_s)])$spec/(length(cases[(i+1):(i+w_s)])) # Spectral variances
    #den2=spectrum0.ar(cases[(i):(i-(w_s-1))])$spec/(length(cases[(i):(i-(w_s-1))]))
    testthat=enu/sqrt(den1+den2)
#    if(test=="ztest"){ # Not large difference between a ztest and a ttest.
#      cevi[k+lag_n+1]=as.numeric((1-pnorm(q = testthat))<=c_n) #*as.numeric(evi[i] >= rate)
#    }
#    if(test=="ttest"){
       cevi[k + lag_n + 1] = as.numeric((1-pt(q = testthat,df = lag_n-1))<=c_n)

#    }
  }
  return(cevi)
}

