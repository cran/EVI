#' Combining multiple epidemic volatility indices into a multi-early warning index plot
#'
#' This function returns an overlap graph containing two indices.
#' 1. If index number 1 results in an early warning,then a yellow point will be shown.
#' 2. If index number 2 results in an early warning, then a orange point will be shown.
#' 3. If both indexes result in early warnings, then a red point will be shown.
#'
#' The evirlap function can literally combine any two indices with only 2 requirements. The two vectors (Index1, Index2) to be of equal length and to contain either 0s (no warnings) or 1s (warnings).
#'
#' For each time point the stored variables are:
#'@return \itemize{
#' \item{evirlap: A graph with two EVI indexes shown simulteounsly.}}

#' @param Index1 e.g. output of the \code{\link[EVI:deviant]{deviant()}} function
#' @param Index2 e.g. output of the \code{\link[EVI:deviant]{deviant()}} function (e.g. method = cEVI)
#' @param ln  TRUE or FALSE; If TRUE (the default) the output of the graph will be presented on the logarithmic scale. IF FALSE the output data will be presented on the original scale.
#' @param type By default, points are plotted on EVI graphs. In cases where, changes are very sudden or data sparsely available, type="l" introduces lines on top of points for the "EVI" type of graph.
#' @param Index1.lab Text for index 1, "EVI1"
#' @param Index2.lab Text for index 2, "EVI2"
#' @param Index3.lab Text for index conjuction, "EVI-"
#' @param Index.country Label for country/area/region under analysis
#' @param size.index Size of the points, default=1
#'
#' @examples
#'
#' # Run basic functions to acquire early warnings results
#' #data('Italy')
#' #EVI_output=deviant(new_cases=Italy$Cases, cum=FALSE, method="EVI")
#' #cEVI_output=deviant(new_cases=Italy$Cases, cum=FALSE, method="cEVI")
#'
#' # Plot the EVI combination graph
#' #evirlap(Index1=EVI_output, Index2=cEVI_output, ln=TRUE)
#' #evirlap(Index1=EVI_output, Index2=cEVI_output, ln=TRUE, type="l") # For the line EVI plot
#'
#' # Austria first 150 case example
#' #data('Austria')
#' #tmp_EVI_at=deviant(new_cases = Austria$Cases)
#' #tmp_cEVI_at=deviant(new_cases = Austria$Cases,lag_max = 40,method = "cEVI")
#' #evirlap(Index1 = tmp_EVI_at,Index2 = tmp_cEVI_at)
#'
#' @export
#'
#' @import ggplot2
#' @import cowplot
#'
#'
#' @references
#' Pateras K, Meletis E, Denwood M, et al. The convergence epidemic index (cEVI) an early warning tool for identifying waves in an epidemic. Inf Dis Mod, (2023)


evirlap <- function(Index1,Index2, ln=T, type="p",size.index=1,
                            Index1.lab="EVI1",Index2.lab="EVI2",Index3.lab="EVI-",Index.country=NULL) {

  Index1$Index=Index1$Index
  Index2$Index=Index2$Index*2
  Index=Index1
  Index$Index=Index$Index+Index2$Index
  if(length(table(Index$Index))<3)
    Index$Index[1:3]<-1:3
  Index$cases_1=Index$Cases*Index$Index
  Index$cases_1[Index$cases_1 == 0] <- NA
  Index$cases_0=Index$Cases*(1-Index$Index)
  Index$cases_0[Index$cases_0 == 0] <- NA

  Index$npv=Index$npv*(1-Index$Index)
  Index$npv[Index$npv == 0] <- NA
  Index$ppv=Index$ppv*Index$Index
  Index$ppv[Index$ppv == 0] <- NA
  Index$variable<-"x"
  Index$Index[is.na(Index$Index)]<-0
  Index$Index<-factor(Index$Index,labels = c("No warning",paste(Index1.lab,"alone"),paste(Index2.lab,"alone"), Index3.lab))
  if (ln==F) {
    sp3<-ggplot(Index, aes_string(x="Days",group="variable"))+
      list(
        geom_point(aes_string(y=("Cases"), color="Index"), size=size.index),
        #scale_color_manual(values=c("grey69", "yellow3", "orange3", "red4")),
        scale_colour_grey(start = 1,end = 0),
        scale_color_manual(values=c("grey69", "yellow3", "orange3", "red4")),
        labs(title = paste0("Graph combining outputs ",Index1.lab,", ", Index2.lab," and ", Index3.lab," - ",Index.country), y = "Cases", x="Days"),
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size=8),
              legend.key.height = unit(0, 'cm')),
        if (type=="l")  geom_path(aes_string(y="Cases",colour="Index"),size=size.index)
      )
  }

  if (ln==T) {
    sp3<-ggplot(Index, aes_string(x="Days",group="variable"))+
      list(
        geom_point(aes_string(y="log(Cases)", color="Index"), size=size.index),
        #scale_color_manual(values=c("grey69", "yellow3", "orange3", "red4")),
        scale_colour_grey(start = 1,end = 0),
        labs(title = paste0("Graph combining outputs ",Index1.lab,", ", Index2.lab," and ", Index3.lab," - ",Index.country), y = "log(Cases)", x="Days"),
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size=8),
              legend.key.height = unit(0, 'cm')),
        if (type=="l")  geom_path(aes_string(y="log(Cases)",colour="Index"), size=size.index)
      )
  }
  print(sp3)

}
