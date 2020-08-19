# ABOUT THIS SCRIPT ===================================================================================================
# This R script prepares Landsat NDVI data for comparison with estimates of gross primary productivity from Arctic flux towers. 
# Author: Logan Berner, NAU
# Date: 2019-01-15

# SET UP WORKSPACE ===================================================================================================
rm(list=ls())
require(lattice)
require(data.table)
require(ranger)

setwd('/projects/above_gedi/lberner/arctic_greening/')
source('/home/lb968/code/arctic_greening/0.1_fun_lsat_tools.R')

lsat <- fread('output/lsat_site_extractions/lsat_fluxnet_tower_samples_20190116.csv')

lsat5.xcal.rf.files <- list.files('output/xcal_ndvi/', pattern = 'LT05_xcal_rf.RData', full.names = T)
lsat8.xcal.rf.files <- list.files('output/xcal_ndvi/', pattern = 'LC08_xcal_rf.RData', full.names = T)


# rename several sites based on updated FLUXNET naming convention ======================================================
lsat[site == 'DK-NuF', site := 'GL-NuF']
lsat[site == 'DK-ZaF', site := 'GL-ZaF']
lsat[site == 'DK-ZaH', site := 'GL-ZaH']
lsat[site == 'NO-Adv', site := 'SJ-Adv']
lsat[site == 'NO-Blv', site := 'SJ-Blv']

# ESTIMATE LANDSAT NDVI MAX ============================================================================================
lsat[,ORIG_FID := NULL]
lsat <- lsat_general_prep(lsat)
lsat <- lsat_qaqc_flags(lsat)
lsat <- lsat_ngb_mean(lsat)
lsat <- lsat_spec_index(lsat, 'ndvi')

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
  
  # apply phenology model
  rep.spar <- 0.7 + runif(1,0.-0.02,0.02)
  lsat.pheno <- lsat_pheno(lsat, 'ndvi.xcal', window.yrs = 17, window.min.obs = 20, spar = rep.spar, spl.fit.outfile = NA)
  
  # growing season max ndvi
  lsat.max <- lsat_pheno_max(lsat.pheno, vi = 'ndvi.xcal', min.frac.of.max = 0.75)
  lsat.max <- lsat.max[order(site, year)]
  
  lsat.max <- lsat.max[, rep := i]
  setnames(lsat.max, 'ndvi.xcal.max.pred','ndvi')
  
  # store to iteration to list
  lsat.max.list[[i]] <- lsat.max
  
  # status
  print(paste0('Finished ', i / mc.reps))
}

lsat.max.mc <- rbindlist(lsat.max.list)


fwrite(lsat.max.mc, 'output/site_comparisons/lsat_ndvi_max_for_flux_tower_sites_mcreps.csv')


# END SCRIPT ============================================================================================ 