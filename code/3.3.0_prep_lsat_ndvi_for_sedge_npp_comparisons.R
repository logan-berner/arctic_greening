# ABOUT THIS SCRIPT ===================================================================================================
# This R script prepares Landsat NDVI data for comparison with estimates of graminoid productivity on Bylot Island 
# Author: Logan Berner, NAU
# Date: 2019-01-15

# SET UP WORKSPACE ===================================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(dplyr)
require(data.table)
require(tidyr)
require(lattice)
require(plotrix)
require(ranger)

setwd('/projects/arctic/users/lberner/arctic_greening/')
source('/home/lb968/code/arctic_greening/0.1_fun_lsat_tools.R')

lsat <- fread('output/lsat_site_extractions/lsat_bylot_bcv_subsite_samples_buf100m.csv', header = T, fill=T)
lsat

lsat5.xcal.rf.files <- list.files('output/xcal_ndvi/', pattern = 'LT05_xcal_rf.RData', full.names = T)
lsat8.xcal.rf.files <- list.files('output/xcal_ndvi/', pattern = 'LC08_xcal_rf.RData', full.names = T)

# ESTIMATE LANDSAT MAXIMUM SUMMER ndvi AT EACH OF THE SHRUB SITES --------------------------------------------------------
lsat <- lsat_general_prep(lsat)
lsat <- lsat_qaqc_flags(lsat)
lsat <- lsat_ngb_mean(lsat)
lsat <- lsat_spec_index(lsat, 'ndvi')
lsat <- lsat[year >= 1988]

# MONTE CARLO LOOP
mc.reps <- 1000
lsat.max.list <- list()
for (i in 1:mc.reps){
  
  # permute reflectance
  LT05.red.scalers <- runif(nrow(lsat[satellite == 'LT05']), -0.07, 0.07)
  LT05.nir.scalers <- runif(nrow(lsat[satellite == 'LT05']), -0.07, 0.07)
  LE07.red.scalers <- runif(nrow(lsat[satellite == 'LE07']), -0.05, 0.05)
  LE07.nir.scalers <- runif(nrow(lsat[satellite == 'LE07']), -0.05, 0.05)
  LC08.red.scalers <- runif(nrow(lsat[satellite == 'LC08']), -0.03, 0.03)
  LC08.nir.scalers <- runif(nrow(lsat[satellite == 'LC08']), -0.03, 0.03)
  
  lsat[satellite == 'LT05', nir := nir + nir * LT05.nir.scalers]
  lsat[satellite == 'LT05', red := red + red * LT05.red.scalers]
  lsat[satellite == 'LE07', nir := nir + nir * LE07.nir.scalers]
  lsat[satellite == 'LE07', red := red + red * LE07.red.scalers]
  lsat[satellite == 'LC08', nir := nir + nir * LC08.nir.scalers]
  lsat[satellite == 'LC08', red := red + red * LC08.red.scalers]
  
  # apply xcal models
  lsat5.xcal.rf <- readRDS(lsat5.xcal.rf.files[i])
  lsat8.xcal.rf <- readRDS(lsat8.xcal.rf.files[i])
  lsat <- lsat[, ndvi.xcal := ndvi]
  lsat5 <- lsat[satellite == 'LT05']
  lsat8 <- lsat[satellite == 'LC08']
  lsat <- lsat[satellite == 'LT05', ndvi.xcal := predict(lsat5.xcal.rf, data = lsat5)$predictions]
  lsat <- lsat[satellite == 'LC08', ndvi.xcal := predict(lsat8.xcal.rf, data = lsat8)$predictions]
  lsat <- lsat[, rep:= i]

  # apply phenology model
  rep.spar <- 0.7 + runif(1,0.-0.02,0.02)
  lsat.pheno <- lsat_pheno(lsat, 'ndvi.xcal', window.yrs = 17, window.min.obs = 20, spar = rep.spar, spl.fit.outfile = NA)
  
  # growing season max ndvi
  lsat.max <- lsat_pheno_max(lsat.pheno, vi = 'ndvi.xcal', min.frac.of.max = 0.75)
  lsat.max <- lsat.max[order(site, year)]
  
  
  # full factorial
  full.fac <- lsat.max %>% expand(site, year = 1988:2017)
  lsat.max <- left_join(full.fac, lsat.max, by = c('site','year'))
  lsat.max <- data.table(lsat.max)
  lsat.max <- lsat.max[, rep := i]
  setnames(lsat.max, 'ndvi.xcal.max.pred','ndvi')
  
  # add multi-year means and lags
  lsat.max <- lsat.max %>% group_by(site) %>% mutate(ndvi.2yr.mean = (ndvi + lag(ndvi,1))/2, ndvi.3yr.mean = (ndvi + lag(ndvi,1) + lag(ndvi,2))/3)
  lsat.max <- lsat.max %>% mutate(ndvi.lag1 = lag(ndvi,1), ndvi.lag2 = lag(ndvi, 2))
  
  # store to iteration to list
  lsat.max.list[[i]] <- lsat.max
  
  # status
  print(paste0('Finished ', i / mc.reps))
}

lsat.max.mc <- rbindlist(lsat.max.list)

unique(lsat.max.mc$site)

fwrite(lsat.max.mc, 'output/site_comparisons/lsat_ndvi_max_for_sedge_npp_subsites_mcreps.csv')

# END SCRIPT ===================================================================================================