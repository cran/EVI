#' Sensitivity-Specificity estimation for each cut-off value and rolling window size
#'
#' This function calculates the sensitivity and the specificity for each cut-off value and rolling window size.
#'
#' @return
#' Returns a list of the estimated Sensitivity, Specificity, apparent and true prevalence for each cut-off value and rolling window size
#'
#' @param evi numeric vector - object (obtained from the evi function and stored as ev) that corresponds to the relative change in the standard deviation.
#' @param cevi numeric vector - object (obtained from the cevi function and stored as cev) that corresponds to either a positive or negative test.
#' @param cases numeric vector with the number of new cases per unit of time (i.e., daily).
#' @param cut threshold value (0 <= c <= 0.5) for issuing an early warning. If evi >= c, an early warning is issued and otherwise is not.
#' @param r Definition for the minimum difference in the mean number of cases, one week before and after each time point that, if present, should be detected. This is the case definition and the default is 0.2 (with 0 <= r <= 1). A value of r=0.2 means that we have a case when the mean number of the newly observed cases in the next 7 days is at least 20% higher than the mean number of the newly observed cases in the past 7 days.
#' @param method either "EVI" or "cEVI", default equals to "EVI".
#'
#' @examples
#' ## EVI Example ##
#' #data("Italy")
#' #cases = mova(cases = Italy$Cases)
#' #roll = rollsd(cases = cases,lag_t = 7)
#' #ev = evi(rollsd = roll)
#' #evifcut(evi = ev, cases = cases, cut = 0.01, r = 0.2)

#' @export
#'
#' @references
#' Kostoulas P, Meletis E, Pateras K, et al. The epidemic volatility index, a novel early warning tool for identifying new waves in an epidemic. Sci Rep 11, 23775 (2021). \doi{10.1038/s41598-021-02622-3}
#' Pateras K, Meletis E, Denwood M, et al. The convergence epidemic index (cEVI) an early warning tool for identifying waves in an epidemic. Inf Dis Mod, (2023)


evifcut=function (evi=NA, cevi=NA, cases, cut=NA, r,method="EVI")
{
  w_s = 7
  ratio = 1/(1 + r)
  test_p = rep(NA, length(cases))
  true_p = rep(NA, length(cases))

  if(method=="EVI"){
    for (i in w_s:(length(cases) - w_s)) {
      if (evi[i] >= cut && cases[i] > mean(cases[i:(i - 6)])) {
        test_p[i] = 1
      }
      else {
        test_p[i] = 0
      }
      if (mean(cases[(i - (w_s - 1)):i],na.rm=T) <= ratio * mean(cases[(i + 1):(i + w_s)],na.rm = T)) {
        true_p[i] = 1
      }
      else {
        true_p[i] = 0
      }
    }
  }

  if(method=="cEVI"){
    for (i in w_s:(length(cases)-w_s)){
      if ((!is.na(cevi[i]) & cevi[i]  == 1) && cases[i]>mean(cases[i:(i-7)])){
        test_p[i] <- 1
      }else{
        test_p[i] <- 0
      }

      cond2<-mean(cases[(i):(i-w_s+3)],na.rm=T) <= ratio * mean(cases[(i+1):(i+w_s)],na.rm=T)
      if (!is.na(cond2) & cond2==TRUE){
        true_p[i] <- 1
      }else{
        true_p[i] <- 0
      }
    }
  }

  sens = length(which(test_p == 1 & true_p == 1))/length(which(true_p == 1))
  spec = length(which(test_p == 0 & true_p == 0))/length(which(true_p == 0))
  ppv = length(which(test_p == 1 & true_p == 1))/length(which(test_p == 1))
  npv = length(which(test_p == 0 & true_p == 0))/length(which(test_p == 0))
  sens[is.nan(sens)] <- 0
  spec[is.nan(spec)] <- 0
  sens[is.nan(sens)] <- 0
  spec[is.nan(spec)] <- 0
  testsin = length(which(test_p == 1))/(length(cases) - w_s)
  prev = length(which(true_p == 1))/(length(cases) - w_s)
  evifcut <- list(sens = sens, spec = spec, testsin = testsin,
                  prev = prev, ppv=ppv, npv=npv)
  return(evifcut)
}
