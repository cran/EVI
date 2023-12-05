#' A function that defines the true status based on the case definition
#'
#' Status = 1 when the expected rise in the number of cases occurs and Status = 0 when the expected rise in the number of cases does not occur
#'
#' @param cases the time series of the newly observed cases per unit of time (ideally per day).
#' @param r	 Definition for the minimum difference in the mean number of cases, one week before and after each time point that, if present, should be detected. This is the case definition and the default is 0.2 (with 0 <= r <= 1). A value of r=0.2 means that we have a case when the mean number of the newly observed cases in the next 7 days is at least 20% higher than the mean number of the newly observed cases in the past 7 days.
#'
#' @examples
#' data("Italy")
#' cases = mova(cases=Italy$Cases)
#' status = status(cases=cases, r=0.2)
#' @export
#'
#'
#'
status = function(cases, r) {
  ratio = 1/(1+r)
  w_s=7
  status=rep(NA,length(cases))
  status[1]=NA
  for (i in 2:(length(cases)-w_s)){
    if (mean(cases[(i-min((i-1),w_s)):(i-1)])<=ratio*mean((cases[i:(i+min(i,(w_s-1)))])))
    {status[i]=1}
    else
    {status[i]=0}
  }
  return(status)
}
