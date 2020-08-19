# ABOUT THIS SCRIT ================================================================================
# R script for calculating trends in gridded geospatial time series (i.e. stacks of rasters). 
# The trend test involves, on a per-pixel basis, pre-whitening the time series, then testing 
# for trend significance and magnitude using the Mann-Kendall and Theil-Sen approaches, respectively. 
# Reference for trend test:
# Zhang, X., L. A. Vincent, W. Hogg, and A. Niitsoo. 2000. Temperature and precipitation trends in Canada during the 20th century. Atmosphere-Ocean 38:395-429.
# Author: Logan Berner, NAU
# Date: 2018-02-16

# CUSTOM FUNCTION ================================================================================
require(raster)
require(zyp)
require(Kendall)

# base function
z.trend <- function(x){
  if (all(is.na(x))){
    out <- c(NA,NA,NA,NA,NA)
    out
  } else{
    trnd <- zyp.zhang(x)
    out <- c(trnd['trend'], trnd['trendp'], trnd['intercept'], trnd['sig'], trnd['tau'])
    out}
}

# wrapper function that also returns trend in pcnt / yr relative to first year in time series
stack.trend <- function(x, p=0.05){
  trend <- calc(x, fun=z.trend)
  names(trend) <- c('trend','trendp','intercept','sig','tau')
  
  # pcnt change / yr
  slp.pcnt <- (trend$trendp/trend$intercept * 100) / nlayers(x) # total change over time period / value during first year / n years in period
  
  #output
  out <- stack(trend, slp.pcnt)
  names(out) <- c('trend','total.change','intercept','pval','tau','trend.pcnt')
  out
}

# END SCRIPT ================================================================================