#' This function is used after first running the deviant function
#'
#' Once the \code{\link[EVI:deviant]{deviant()}} function has been used to analyze the already observed time series,
#' the deviant_update() function is used to obtain the EVI output and early warnings for the new cases that are recorded.
#' After running the deviant_update() function the output of the deviant function is also updated with a new row of data for each newly observed time point.
#'
#' @param all_cases the time series of the previous and newly observed cases per unit of time (ideally per day).
#' @param new_cases the time series of the newly observed cases per unit of time (ideally per day).
#' @param EVI_input the previous EVI function output  
#' @param cum TRUE if the time series is recorded as the cumulative number of the reported cases and FALSE (the default) if newly reported cases per unit of time are recorded.
#' @param r_a The window size for the moving average that will be analyzed. If set to 1 the actual observations are analyzed. However, due to the variability of the reported cases between working days and weekends it is recommended that the 7-day moving average is analyzed (i.e. r_a = 7), which is the default for this argument. Users could prefer a longer interval of 14 days or one month (i.e., r_a=14 or 30, respectively).
#' @param r Definition for the minimum difference in the mean number of cases, one week before and after each time point that, if present, should be detected. This is the case definition and the default is 0.2 (with 0 <= r <= 1). A value of r=0.2 means that we have a case when the mean number of the newly observed cases in the next 7 days is at least 20% higher than the mean number of the newly observed cases in the past 7 days.
#' @param lag_max Integer. Restriction of the maximum window size for the rolling window size. The default is set to one month (lag_max=30) to prevent excess volatility of past epidemic waves from affecting the most recent volatility estimates and the ability of EVI to warn for upcoming waves that may be smaller and of lower volatility than previous ones.
#' @param past Integer. Default the full length of the input data series. Restriction on the historical data that EVI/cEVI will use. This is fixed and equal to the index used in the deviant function. 
#' @param method either "EVI" or "cEVI", default equals to "EVI".

#'
#' @examples
#' ## Single update ##
#' # If we have first observed only the 148 cases from the 
#' # Italian data we run the deviant function on these cases:
#'
#' data("Italy")
#' EVI_output0<-deviant(new_cases=Italy$Cases[20:35], cum=FALSE, r_a=7, r=0.2, lag_max=30)
#' # EVI_output0<-deviant(new_cases=Italy$Cases[1:148], cum=FALSE, r_a=7, r=0.2, lag_max=30)

#' # When the number of cases for the next day is observed we want to obtain 
#' # the estimates for this day without having to reanalyze the entire time series.
#'  
#' # This is done by using the deviant_update function:
#'
#' EVI_output0<-deviant_update(new_cases=Italy$Cases[36], 
#'                             EVI_input=EVI_output0,
#'                             cum=FALSE, r_a=7, r=0.2, lag_max=30)
#'
#' # EVI_output0<-deviant_update(new_cases=Italy$Cases[149], 
#' #                             EVI_input=EVI_output0,
#' #                             cum=FALSE, r_a=7, r=0.2, lag_max=30)
#' # The result of running the deviant_update function is to update the output of 
#' # the deviant_function by adding an additional row with estimates for the new data.
#' 
#' # In this example the EVI_output file will now have 149 rows. If two additional 
#' # days are analyzed two additional rows will be added and so on.
#' 
#' 
#' 
#' # Multiple update steps when new data come using either an EVI or cEVI updater ##
#'data("Italy")
#'EVI_output<-deviant(new_cases=Italy$Cases[1:20], cum=FALSE, 
#'                    r_a=7, r=0.2, lag_max=30,method = "EVI")
#'
#'EVI_output2<-deviant_update(new_cases = c(100,93,80,54,12), 
#'                            EVI_input=EVI_output,method = "EVI")
#'
#'#EVI_output2
#'
#'# Same as above EVI_output2
#'
#'#EVI_output2<-deviant_update(all_cases = c(Italy$Cases[1:20],100,93,80,54,12),
#'#EVI_input=EVI_output,method = "EVI")
#'
#'#EVI_output2
#'
#'#EVI_output3<-deviant_update(new_cases = c(2,2,10,1,0),EVI_input=EVI_output2, 
#'#method = "cEVI")
#'
#'#EVI_output3
#'
#'# Even though EVI and cEVI can be used interchangeably, we suggest users to 
#'# stick to the initial method.
#'
#' @importFrom stats pt
#' @export
#'
#' @references
#' Pateras K, Meletis E, Denwood M, et al. The convergence epidemic index (cEVI) an early warning tool for identifying waves in an epidemic. Inf Dis Mod, (2023). \doi{10.1016/j.idm.2023.05.001}
#' Kostoulas P, Meletis E, Pateras K, et al. The epidemic volatility index, a novel early warning tool for identifying new waves in an epidemic. Sci Rep 11, 23775 (2021). \doi{10.1038/s41598-021-02622-3}


deviant_update=function(all_cases=NA, new_cases=NA,
                        EVI_input, cum = FALSE, r_a=7,
                        r=0.2, lag_max=30, method="EVI",past = 364/2){
  
  
  if(method=="EVI"){
    start_cases=14
    lag_1=7
    c_1=0.01
    w_s =7
  }else if (method=="cEVI"){
    start_cases=18
    lag_1=3
    c_1=0.001
    w_s=7
  }
  
  if (cum == TRUE) new_cases = c(new_cases[1], diff(new_cases))
  
  if ((!is.na(sum(all_cases)) & !is.na(sum(new_cases))) | (is.na(sum(all_cases)) & is.na(sum(new_cases))))
    stop("Please provide either new cases or all cases (old+new) as input")
  
  if(!is.na(sum(all_cases))){
    #calculate the moving average of new confirmed cases
    cases=mova(c(all_cases),r_a)
  }
  if(!is.na(sum(new_cases))){
    #calculate the moving average of new confirmed cases
    cases=mova(c(EVI_input$new_cases,new_cases),r_a)
  }
  
  
  
  roll=rollsd(cases[1:start_cases],lag_1)
  ev=evi(roll)
  if(method=="EVI"){
    ind=indic(evi = ev, cut = c_1, cases = cases[1:start_cases], method = method)
  }else if (method=="cEVI"){
    cevi=cEVI_fun(cases = cases[1:(start_cases)],lag_n = lag_1, c_n = c_1)
    ind=indic(cevi=cevi, cases=cases[1:start_cases], method="cEVI")
  }
  status=status(cases[1:start_cases],r)
  
  #initiate chain for positive predictive value
  ppv=rep(NA, length(cases))
  
  #initiate chain for negative predictive value
  npv=rep(NA, length(cases))
  
  lag_all=rep(NA, start_cases)
  c_all=rep(NA, start_cases)
  
  se_all=rep(NA, start_cases)
  sp_all=rep(NA, start_cases)
  
  
  lag_all[1:start_cases]=lag_1
  c_all[1:start_cases]=c_1
  
  diff= length(cases)-(nrow(EVI_input) +1)
  
  for (i in (nrow(EVI_input)+1): length(cases)){
    case_t=cases[max(1,(i-past)):i]
    lag_s=seq(lag_1,min(lag_max,(length(case_t)-1)), 1)
    c_s=seq(0.01,0.5, 0.01)
    all_lag<-all_cut<-all_se<-all_sp<-NA
    
    if(method=="EVI"){
      
      for (j in lag_s){
        roll_t=rollsd(case_t,j)
        ev_t=evi(roll_t)
        for (l in c_s){
          evicut_t=evifcut(evi = ev_t, cases = case_t, cut = l, r = r,method = "EVI")
          new_j=j
          new_l=l
          new_se=evicut_t$sens
          new_sp=evicut_t$spec
          all_lag[[length(all_lag) + 1]] <- new_j
          all_cut[[length(all_cut) + 1]] <- new_l
          all_se[[length(all_se) + 1]] <- new_se
          all_sp[[length(all_sp) + 1]] <- new_sp
        }
      }
      
    }
    
    if(method=="cEVI"){
      i_n<-max(1,(i-past)):i
      case_t=cases[i_n]
      lag_s=seq(lag_1,min(lag_max,(length(i_n)-length(i_n)/2-4)), 2)
      c_s=seq(0.001,0.5, 0.06)
      all_lag<-all_cut<-all_se<-all_sp<-NA
      
      
      for (l in c_s) {
        for (j in lag_s) {
          # roll_t <- rollsd(case_t,j)
          #  ev_t <- evi(roll_t)
          cevi <- rep(NA, length(case_t))
          for(k in (j+1):(length(case_t)-(j+1))){
            enu=mean(case_t[(k+2):(k+j+1)]-case_t[(k):(k-(j-1))],na.rm = T)
            den1=sd(case_t[(k):(k-(j-1))])^2/(length(case_t[(k):(k-(j-1))]))
            den2=sd(case_t[(k+2):(k+j+1)])^2/(length(case_t[(k+2):(k+j+1)]))
            teststat=enu/sqrt(den1+den2)
            Nn=length((k+1):(k+j))
            cevi[k+j+1]<-as.numeric((1-stats::pt(q = teststat,df = Nn))<=l)
          }
          evicut_t <- evifcut(cevi=cevi,cases = case_t, r = r,method = "cEVI")
          all_lag[[length(all_lag) + 1]] <- j
          all_cut[[length(all_cut) + 1]] <- l
          all_se[[length(all_se) + 1]] <- evicut_t[[1]]
          all_sp[[length(all_sp) + 1]] <- evicut_t[[2]]
        }
      }
    }
    
    sesp=as.data.frame(cbind(all_lag,all_cut,all_se,all_sp))
    
    #Select the row with the right window and cut
    index=which.max(sesp$all_se+sesp$all_sp-1)
    print(i);print(sesp[index,])
    
    #estimate the parameters for the last observed case
    lag_n=sesp$all_lag[index]
    c_n=sesp$all_cut[index]
    

    if (method=="EVI"){
      roll_n=rollsd(case_t,lag_n)
      ev_n=evi(roll_n)
      ind_n=indic(evi = ev_n,cut = c_n, cases = case_t, method=method)
      evicut_n=evifcut(evi = ev_n, cases = case_t, cut = c_n, r = r, method=method)
      roll=c(roll,roll_n[length(ind_n)])
    }else if (method=="cEVI"){
      ev_n=cEVI_fun(cases = case_t,lag_n = lag_n, c_n = c_n) #
      ind_n=indic(cevi = ev_n,cut = c_n , cases = case_t, method=method) #
      evicut_n=evifcut(cevi = ev_n, cases = case_t, r = r, method=method) #
    }
    
    ev=c(ev,ev_n[length(ind_n)])
    ind=c(ind, ind_n[length(ind_n)])
    
    lag_all=c(lag_all,lag_n)
    c_all=c(c_all,c_n)
    
    se_all=c(se_all,all_se[index])
    sp_all=c(sp_all,all_sp[index])
    
    ppv[i]=evicut_n$prev*all_se[index]/
      (evicut_n$prev*all_se[index]+(1-evicut_n$prev)*(1-all_sp[index]))
    
    npv[i]=(1-evicut_n$prev)*all_sp[index]/
      ((1-evicut_n$prev)*all_sp[index]+evicut_n$prev*(1-all_se[index]))
  }
  
  Days=((length(cases)-diff):length(cases))
  EVI=ev[((length(ev)-diff):length(ev))]
  Cases=cases[((length(cases)-diff):length(cases))]
  Index=ind[((length(ind)-diff):length(ind))]
  ppv=ppv[((length(ppv)-diff):length(ppv))]
  npv=npv[((length(npv)-diff):length(npv))]
  lag_all=lag_all[((length(lag_all)-diff):length(lag_all))]
  c_all=c_all[((length(c_all)-diff):length(c_all))]
  se_all=se_all[((length(se_all)-diff):length(se_all))]
  sp_all=sp_all[((length(sp_all)-diff):length(sp_all))]
  
  
  EVI_out_add=as.data.frame(cbind(Days, EVI, Cases, Index, ppv, npv,
                                  lag_all, c_all, se_all, sp_all,new_cases))
  
  EVI_input=rbind(EVI_input,EVI_out_add)
  
  EVI_output<-(EVI_input)
  
  return(EVI_output)
  
}
