# ABOUT THIS SCRIPT =========================================================================
# This R script takes Landsat data downloaded from GEE and then post-processess the data to 
# estimate NDVImax for 50,000 sites randomly located throughtout the Arctic tundra biome. 
# Date: 2018-10-17
# Author: Logan Berner, Northern Arizona University

# SET UP WORKSPACE ============================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(tidyr)
require(data.table)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ranger)
require(R.utils)
setwd('/projects/above_gedi/lberner/arctic_greening/')
source('/home/lb968/code/arctic_greening/0.1_fun_lsat_tools.R')

args <- commandArgs(TRUE)
i = as.numeric(args[1])
# i = 1

# reformat i for sorting
if (i < 10){
  i <- paste0('000',i)
} else if (i < 100){
  i <- paste0('00',i)
} else if (i < 1000){
  i <- paste0('0',i)
}

# # RANDOMLY SELECT 50K SAMPLING SITES  ================================================================  
# lsat <- fread('data/lsat_samples/lsat_tundra_samples.csv')
# lsat <- lsat[year <= 2016]
# lsat <- lsat[ndvi > 0]
# lsat.sites <- unique(lsat$site)
# lsat.sites.50k <- lsat.sites[sample(1:length(lsat.sites), size = 50000, replace = F)]
# lsat <- lsat[lsat$site %in% lsat.sites.50k]
# length(unique(lsat$site)) # n sites
# # nrow(lsat)
# fwrite(lsat, 'data/lsat_samples/lsat_tundra_samples_n50000.csv')

# LOAD LANDSAT DATA  ================================================================  
lsat <- fread('data/lsat_samples/lsat_tundra_samples_n50000.csv')
lsat.sites <- unique(lsat$site)
length(lsat.sites)
# lsat.sites.50k <- lsat.sites[sample(1:length(lsat.sites), size = 1000, replace = F)]
# lsat <- lsat[lsat$site %in% lsat.sites.50k]


# PERMUTE REFLECTANCE AND CALC NDVI ================================================================
# sensor specific uncertanity from Markham and Helder (2012, RSE) and Markham et al. (2014; RS)
# mkdirs('output/reflectance_scalars')
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


# calculate NDVI
lsat <- lsat_spec_index(lsat, 'ndvi')

# CROSS-CALIBRATE NDVI AMONG LANDSAT SENSORS ================================================================
print('starting xcal...')
mkdirs('output/xcal_ndvi/')
lsat <- lsat_xcal_rf(lsat, 'ndvi', doy.rng = 152:243, min.obs = 5, frac.eval = 0.33, outfile.prefix = paste0('ndvi_rep_',i), outdir = 'output/xcal_ndvi')


# FILTER OUT BARREN SITES AND SITES WITH TOO FEW YEARS OF OBSERVATIONS  ================================================================ 
# subset barren sites
lsat <- lsat[, ndvi.xcal.avg := mean(ndvi.xcal, na.rm = T), by = site]
lsat.barren <- lsat[ndvi.xcal.avg < 0.10]
# fwrite(lsat.barren, 'output/tundra_barren_sites.csv')

length(unique(lsat.barren$site)) # number of barren sites
lsat <- lsat[ndvi.xcal.avg >= 0.10][ndvi.xcal >= 0.10]
length(unique(lsat$site)) # after excluding barren

# identify sites with too short a record 
lsat <- lsat[, n.yrs := length(unique(year)), by = site]
lsat.shortrecords <- lsat[n.yrs < 10]
lsat <- lsat[n.yrs >= 10]
length(unique(lsat$site)) # after excluding short records

# identify sites with too few observations
lsat <- lsat[, n.obs := .N, by = site]
lsat <- lsat[n.obs > 30]
length(unique(lsat$site)) # after excluding short records


# ESTIMATE ANNUAL MAXIMUM SUMMER NDVI ================================================================ 
print('starting pheno...')
# apply phenological correction
mkdirs('output/pheno_curves')
mkdirs('output/pheno_timeseries')

# derive phenological curves for each site
rep.spar <- 0.7 + runif(1,0.-0.02,0.02)
spl.outfile <- paste0('output/pheno_curves/pheno_curves_rep_',i,'.csv')
lsat.pheno <- lsat_pheno(lsat, 'ndvi.xcal', window.yrs = 17, window.min.obs = 20, spar = rep.spar, spl.fit.outfile = spl.outfile)
fwrite(lsat.pheno, paste0('output/pheno_timeseries/tundra_landsat_ndvi_pheno_corrected_timeseries_rep_',i,'.csv'))

# compute max NDVI
lsat.max <- lsat_pheno_max(lsat.pheno, vi = 'ndvi.xcal', min.frac.of.max = 0.75)

# evaluate estimates of max NDVI
lsat.max.eval <- lsat_pheno_max_eval(lsat.pheno, vi = 'ndvi.xcal', min.frac.of.max = 0.75, min.obs = 11, reps = 10, 
                                     outdir = 'output/pheno_max_eval/', outfile.suffix = paste0('rep_',i))

# rename some of the NDVI colunms
colnames(lsat.max) <- gsub('.pred','',gsub('.pred.min','.lower',gsub('.pred.max','.upper',gsub('.xcal','',colnames(lsat.max)))))

# DETERMINE OBS DENSITY ==============================================================================================================
# determine first year of obs and number of years with obs
lsat.max <- lsat.max[, first.yr := min(year, na.rm=T), by = site]
lsat.max <- lsat.max[, n.yr.obs := length(unique(year)), by = site]
length(unique(lsat.max$site))

# convert implicit NA to explicit NA (i.e., ensure all site x year combinations)
full.fac <- lsat.max %>% expand(site, year = 1982:2016)
full.fac <- data.table(full.fac)
lsat.max <- lsat.max[full.fac, on = c('site','year')]

length(unique(lsat.max$site))

# WRITE OUT TIME SERIES ============================================================================================================== 
mkdirs('output/ndvi_max_timeseries')
fwrite(lsat.max, paste0('output/ndvi_max_timeseries/tundra_site_lsat_ndvi_max_timeseries_rep_',i,'.csv'))

# END SCRIPT  ============================================================================================================== 
print('DONE!!!')