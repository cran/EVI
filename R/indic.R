#' Issue of an Early Warning
#'
#' This secondary function produces the early warning signal (Index).
#'
#' @return
#' A vector of 0s and 1s is produced, where a 1 (Index = 1) is recorded when an early warning is issued and a 0 (Index = 0) when an early warning is not issued.
#'
#' @param evi numeric vector - object (obtained from the evi function and stored as ev) that corresponds to the relative change in the standard deviation.
#' @param cevi numeric vector - object (obtained from the cevi function and stored as cev) that corresponds to a positive and negative test.
#' @param cut threshold value (0 <= c <= 0.5) for issuing an early warning. If evi >= c an early warning is issued and otherwise is not.
#' @param cases numeric vector with the number of new cases per unit of time (i.e., daily).
#' @param method either "EVI" or "cEVI", default equals to "EVI".
#'
#' @examples
#' ## EVI example ##
#'
#' data("Italy")
#' cases = mova(cases = Italy$Cases, r_a = 7)
#' roll = rollsd(cases = cases, lag_t = 7)
#' ev = evi(rollsd = roll)
#' ind=indic(evi = ev, cut = 0.01, cases = cases)
#'
#' @export
#'
#' @references
#' Kostoulas P, Meletis E, Pateras K, et al. The epidemic volatility index, a novel early warning tool for identifying new waves in an epidemic. Sci Rep 11, 23775 (2021). \doi{10.1038/s41598-021-02622-3}
#' Pateras K, Meletis E, Denwood M, et al. The convergence epidemic index (cEVI) an early warning tool for identifying waves in an epidemic. Inf Dis Mod, (2023)

indic=function (evi=NA, cevi=NA, cut=NA, cases, method="EVI")
{

  if(method=="EVI"){
    ind = rep(NA, length(evi))
    for (i in 3:length(evi)) {
      if (evi[i] >= cut && cases[i] > mean(cases[i:(i - min(7,i))])){
        ind[i] = 1
      }
      else {
        ind[i] = 0
      }
    }
  }

  if(method=="cEVI"){
    ind <- rep(NA,length(cevi))

    for (i in 3:length(cevi)) {
      if (!is.na(cevi[i]) && cevi[i] ==1 && (!is.na(cases[i]) & cases[i] > mean(cases[i:(i - min(7, i))]))){
        ind[i] <- 1
      }else{
        ind[i] <- 0
      }
    }
  }

  return(ind)
}
