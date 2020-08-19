# ABOUT THIS SCRIPT ===================================================================================================
# This R script prepares Landsat NDVI data for comparison with shrub ring width measuremts from sites around the Arctic. 
# Author: Logan Berner, NAU
# Date: 2019-01-15

# SET UP WORKSPACE ===================================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.6.2'))
require(dplyr)
require(dplR)
require(data.table)
require(lattice)
require(plotrix)
require(ranger)
require(R.utils)

setwd('/projects/arctic/users/lberner/arctic_greening/')
source('/home/lb968/code/arctic_greening/0.1_fun_lsat_tools.R')

args <- commandArgs(TRUE)
i = as.numeric(args[1])

# LOAD RAW LANDSAT DATA, CLEANUP, AND WRITE OUT ===================================================================================================
# DO THIS STEP-WISE SINCE IT TAKES ~35 GB RAM AND AN HOUR TO RUN
# lsat <- fread('output/lsat_site_extractions/lsat_tundra_shrub_cluster_sites_buf100m_20190507.csv', fill=T)
# lsat <- lsat_general_prep(lsat)
# lsat <- lsat_qaqc_flags(lsat)
# lsat <- lsat_ngb_mean(lsat)
# lsat <- lsat_spec_index(lsat, 'ndvi')
# lsat <- lsat[ndvi > 0.20]
# lsat <- lsat[year <= 2017]
# lsat <- lsat[, n.obs := .N, by = site][n.obs > 20]
# fwrite(lsat, 'output/lsat_site_extractions/lsat_tundra_shrub_cluster_sites_buf100m_cleaned.csv')

# LOAD DATA SETS ===================================================================================================
sites.dt <- fread('data/field_data/tundra_shrub_growth/tundra_shrub_cluster_site_coords_20190507.csv', fill=T)
lsat <- fread('output/lsat_site_extractions/lsat_tundra_shrub_cluster_sites_buf100m_cleaned.csv')

lsat5.xcal.rf.files <- list.files('output/xcal_ndvi/', pattern = 'LT05_xcal_rf.RData', full.names = T)
lsat8.xcal.rf.files <- list.files('output/xcal_ndvi/', pattern = 'LC08_xcal_rf.RData', full.names = T)

lsat5.xcal.rf <- readRDS(lsat5.xcal.rf.files[i])
lsat8.xcal.rf <- readRDS(lsat8.xcal.rf.files[i])

# reformat i for sorting
if (i < 10){
  i <- paste0('000',i)
} else if (i < 100){
  i <- paste0('00',i)
} else if (i < 1000){
  i <- paste0('0',i)
}

# ESTIMATE LANDSAT MAXIMUM SUMMER NDVI FOR EACH SHRUB SITE CLUSTER IN MONTE CARLO FRAMEWORK ============================
  
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
lsat <- lsat[, ndvi.xcal := ndvi]
lsat5 <- lsat[satellite == 'LT05']
lsat8 <- lsat[satellite == 'LC08']
lsat <- lsat[satellite == 'LT05', ndvi.xcal := predict(lsat5.xcal.rf, data = lsat5)$predictions]
lsat <- lsat[satellite == 'LC08', ndvi.xcal := predict(lsat8.xcal.rf, data = lsat8)$predictions]

# apply phenology model
rep.spar <- 0.7 + runif(1,0.-0.02,0.02)
lsat.pheno <- lsat_pheno(lsat, 'ndvi.xcal', window.yrs = 17, window.min.obs = 20, spar = rep.spar, spl.fit.outfile = NA)

# calc growing season max ndvi
lsat.max <- lsat_pheno_max(lsat.pheno, vi = 'ndvi.xcal', min.frac.of.max = 0.75)
lsat.max <- lsat.max[order(site, year)]
setnames(lsat.max, 'ndvi.xcal.max.pred','ndvi')

# detrend max ndvi
lsat.max <- lsat.max[, n.yrs := .N, by = 'site']
lsat.max <- lsat.max[, spl.fit := ffcsaps(ndvi, year, nyrs = runif(1, 15, 25), f = runif(1, 0.45, 0.55)), by = 'site']
lsat.max <- lsat.max[, ndvi.dt := ndvi / spl.fit]

# add meta data
lsat.max <- sites.dt[lsat.max, on = 'site']

# The Sagwan_N Salix site is sandwitched between a road and river (both w/in 100 m of the point) so use NDVI from ~500 m away at the Alder site 
sag <- lsat.max[cluster == 'Sagwon_N' & genus == 'Alnus'][, genus := 'Salix']
lsat.max <- lsat.max[(cluster == 'Sagwon_N' & genus == 'Salix') == F]
lsat.max <- rbind(lsat.max, sag)

# summarize ndvi across clusters of shrub ring-width sites 
lsat.max.clst <- lsat.max[, .(ndvi.med = median(ndvi, na.rm = T), ndvi.dt.med = median(ndvi.dt, na.rm = T), 
                              n.obs.avg = mean(n.obs), n.sample.sites = .N),
                          by = c('cluster','year','genus')]

# write out iteration
mkdirs('output/lsat_site_extractions/lsat_shrub_site_mcreps/')
lsat.max.clst <- lsat.max.clst[, rep := i]

fwrite(lsat.max.clst, paste0('output/lsat_site_extractions/lsat_shrub_site_mcreps/lsat_ndvi_max_for_shrub_rw_clusters_rep_',i,'.csv'))

# END SCRIPT  ===================================================================================================
