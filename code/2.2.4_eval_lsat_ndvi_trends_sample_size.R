# ABOUT THIS SCRIPT  ==============================================================================================================
# This R script assesses the effects of sample size on estimates of NDVI trends in the Arctic
# AUTHOR: LOGAN BERNER, NAU
# DATE: 2020-04-16

# SET UP WORKSPACE ==============================================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
library(dplyr)
require(tidyr)
library(data.table)
library(reshape2)
library(zyp)
library(R.utils)

args <- commandArgs(TRUE)
i = as.numeric(args[1])
# i = 1

setwd('/projects/above_gedi/lberner/arctic_greening/')

# WRAPPER FUNCTIONS FOR COMUPTING AND SUMMARIZING TRENDS IN NDVI ACROSS SITES ==================================================
calc.trends <- function(x,y){
  xx <- zyp.yuepilon(y,x) ## note the order of x and y are switched in this call!!!
  return(data.table(slope=round(xx['trend'],5), int=round(xx['intercept'],5), tau=round(xx['tau'],3), pval=round(xx['sig'],4)))
}


# READ IN LANDSAT VI TIME SERIES  ==============================================================================================================
site.dt <- fread('output/tundra_site_conditions.csv')

files <- list.files('output/ndvi_max_timeseries/', full.names = T)
site.ts.dt <- fread(files[i], fill=T)
site.ts.dt <- site.ts.dt[year >= 2000][first.yr <= 2001]
site.ts.dt[, year.rsc := year - 2000]
setnames(site.ts.dt, 'ndvi.max','ndvi')

# reformat i for sorting
if (i < 10){
  i <- paste0('000',i)
} else if (i < 100){
  i <- paste0('00',i)
} else if (i < 1000){
  i <- paste0('0',i)
}

# create output lists
biome.trnd.list <- list()
site.trnd.list <- list()

# specify sample sizes to test
sample.sizes <- c(seq(100,1000,100), seq(2000,40000,1000))

# iteratie through sample sizes
cnt = 1
for (j in sample.sizes){
  
  # select data from given number of sample sites
  sites <- unique(site.ts.dt$site)
  site.ts.ss.dt <- site.ts.dt[site %in% sites[sample(length(sites), j)]]

  # normalize NDVI data
  site.ts.ss.dt <- site.ts.ss.dt[, ':='(ndvi.avg = mean(ndvi, na.rm=T), ndvi.sd = sd(ndvi, na.rm=T)), by = c('site')]
  site.ts.ss.dt <- site.ts.ss.dt[, ndvi.anom := ndvi - ndvi.avg]
  site.ts.ss.dt <- site.ts.ss.dt[, ndvi.zanom := ndvi.anom / ndvi.sd]
  
  # compute trend in mean biome ndvi and store to list
  biome.ts.ss.dt <- site.ts.ss.dt[, .(ndvi.avg=mean(ndvi, na.rm = T), ndvi.sd=sd(ndvi, na.rm = T) , n.sites = .N,
                                           ndvi.anom.avg=mean(ndvi.anom, na.rm = T), ndvi.anom.sd=sd(ndvi.anom, na.rm = T),
                                           ndvi.zanom.avg=mean(ndvi.zanom, na.rm = T), ndvi.zanom.sd=sd(ndvi.zanom, na.rm = T),
                                           lat.zone = 'biome'), 
                                       by = c('year')]
  
  biome.trnd.ss.dt <- biome.ts.ss.dt %>% group_by() %>% do(out=calc.trends(x=.$year, y=.$ndvi.avg)) %>% unnest(cols=c(out)) %>% data.table()
  biome.trnd.ss.dt[, ':='(rep = i, sample.size = j)]
  biome.trnd.list[[cnt]] <- biome.trnd.ss.dt
  
  # compute trend for each site, categorize, compute frequency, and store to list
  site.trnd.ss.dt <- site.ts.ss.dt %>% group_by(site) %>% do(out=calc.trends(x=.$year, y=.$ndvi)) %>% unnest(cols=c(out)) %>% data.table()
  site.trnd.ss.dt <- site.trnd.ss.dt[is.na(slope) == F]
  site.trnd.ss.dt[, trend.cat := cut(slope, c(-Inf, 0, Inf), c('browning','greening'))]
  site.trnd.ss.dt[pval >= 0.10001, trend.cat := 'insig']
  site.trnd.freq.ss.dt <- site.trnd.ss.dt[, .(n.sites = .N), by=c('trend.cat')] # frequency
  site.trnd.freq.ss.dt[, n.sites.tot := sum(n.sites)][, pcnt.sites := (n.sites / n.sites.tot) * 100]
  site.trnd.freq.ss.dt[, ':='(rep = i, sample.size = j)]
  site.trnd.list[[cnt]] <- site.trnd.freq.ss.dt
  
  # status
  cnt = cnt + 1
  print(paste0('finished sample size ', j))  
}  

# combine output into one data table
biome.trnd.smry.dt <- rbindlist(biome.trnd.list)
site.trnd.smry.dt <- rbindlist(site.trnd.list)

# write out
mkdirs('output/lsat_sample_size_test/trend_freq_mcreps')
mkdirs('output/lsat_sample_size_test/trend_mag_mcreps')

fwrite(site.trnd.smry.dt, paste0('output/lsat_sample_size_test/trend_freq_mcreps/site_trend_freq_rep_',i,'.csv'))
fwrite(biome.trnd.smry.dt, paste0('output/lsat_sample_size_test/trend_mag_mcreps/biome_trend_rep_',i,'.csv'))

print("All done!!")
# END SCRIPT ================================================================================================================