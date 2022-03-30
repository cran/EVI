## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(EVI)

## ---- results='hide'----------------------------------------------------------
data("Italy")

## ---- results="hide"----------------------------------------------------------
deviant(new_cases = Italy$Cases[1:49], cum = FALSE, r_a = 7, r = 0.2, lag_max = 30)

## ---- results="hide"----------------------------------------------------------
deviant_update(new_cases=Italy$Cases[1:50], cum=FALSE, r_a=7, r=0.2, lag_max=30)

## ---- results="hide", include=FALSE-------------------------------------------
evi.graphs(EVI_output=EVI_output, graph="EVI", ln=T) # confirmed cases plotted on the logarithmic scale
evi.graphs(EVI_output=EVI_output, graph="PPV", ln=F) # Positive predictive value plotted on the original scale
evi.graphs(EVI_output=EVI_output, graph="NPV", ln=T) # Negative predictive value plotted on the logarithmic scale
evi.graphs(EVI_output=EVI_output, graph="EVI", ln=T, type="l") # confirmed cases plotted on the logarithmic scale with lines on top of points

## ---- results="hide", eval=FALSE----------------------------------------------
#  deviant(new_cases = Italy$Cases, cum = FALSE, r_a = 7, r = 0.2, lag_max = 30)

## -----------------------------------------------------------------------------
citation("EVI")

## -----------------------------------------------------------------------------
sessionInfo()

