# ABOUT THIS SCRIPT ====================================================================================================
# This R script computes correlations between time series of NDVI and SWI for sites and zones in the Arctic.
# Author: Logan Berner, Northern Arizona University
# DATE: 2020-04-02

# SET UP WORKSPACE ====================================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(data.table)
require(dplyr)
require(tidyr)
require(raster)
require(ggplot2)
require(R.utils)

args <- commandArgs(TRUE)
i = as.numeric(args[1])
# i = 2

setwd('/projects/arctic/users/lberner/arctic_greening/')

# READ IN LANDSAT NDVI AND SWI TIME SERIES, PREP FOR CORRELATIONS ====================================================================================================
clim.dt <- fread(list.files('output/swi_site_timeseries/mc_reps/', full.names = T)[i])
lsat.dt <- fread(list.files('output/ndvi_max_timeseries/', full.names = T)[i])
cond.dt <- fread('output/tundra_site_conditions.csv')

# rename NDVI variable and drop NAs 
setnames(lsat.dt, 'ndvi.max','ndvi')
lsat.dt <- lsat.dt[!is.na(ndvi)]

# compute 2-yr and 3-yr average SWI
clim.dt[, swi.2yr.avg := movingFun(swi, 2, mean, type = 'to', na.rm = T), by = 'site']
clim.dt[, swi.3yr.avg := movingFun(swi, 3, mean, type = 'to', na.rm = T), by = 'site']

# combine lsat and clim data sets
site.dt <- clim.dt[lsat.dt, on = c('site','year')]
site.dt <- site.dt[!is.na(swi)]
cond.dt <- cond.dt[,c('latitude','longitude','lat.zone'):=NULL]
site.dt <- cond.dt[site.dt, on = c('site')]
site.dt <- site.dt[lat.zone != '']
  
# split data into two time periods, then recombine
site.gte1985.dt <- site.dt[, first.yr := min(year, na.rm=T), by = site][first.yr <= 1986][year >= 1985]
site.gte1985.dt <- site.gte1985.dt[, n.yrs := .N, by = c('site')][n.yrs >= 10]
site.gte1985.dt <- site.gte1985.dt[, period := '1985-2016']

site.gte2000.dt <- site.dt[, first.yr := min(year, na.rm=T), by = site][first.yr <= 2001][year >= 2000]
site.gte2000.dt <- site.gte2000.dt[, n.yrs := .N, by = c('site')][n.yrs >= 10]
site.gte2000.dt <- site.gte2000.dt[, period := '2000-2016']
  
site.dt <- rbind(site.gte1985.dt, site.gte2000.dt)
  
# mean-center each time series
site.dt <- site.dt[, ndvi.anom := scale(ndvi, center = T, scale = F), by = c('site', 'period')]
site.dt <- site.dt[, swi.anom := scale(swi, center = T, scale = F), by = c('site', 'period')]
site.dt <- site.dt[, swi.2yr.avg.anom := scale(swi.2yr.avg, center = T, scale = F), by = c('site', 'period')]

# detrend time series
site.dt[, ndvi.anom.dt := lm(ndvi.anom ~ year)$resid, by=c('site','period')]
site.dt[, swi.anom.dt := lm(swi.anom ~ year)$resid, by=c('site','period')]
site.dt[, swi.2yr.avg.anom.dt := lm(swi.2yr.avg.anom ~ year)$resid, by=c('site','period')]

# reformat i for naming output
if (i < 10){
  i <- paste0('000',i)
} else if (i < 100){
  i <- paste0('00',i)
} else if (i < 1000){
  i <- paste0('0',i)
}

print('finished data prep...')  


# CORRELATE NDVI AND SWI FOR EACH SAMPLING SITE =======================================================================

# compute correlations
swi.anom.site.cor.dt <- site.dt[, .(r = cor.test(ndvi.anom, swi.anom, method = 'spearman')$estimate,
                                    p = cor.test(ndvi.anom, swi.anom, method = 'spearman')$p.val,
                                    r.2yr = cor.test(ndvi.anom, swi.2yr.avg.anom, method = 'spearman')$estimate,
                                    p.2yr = cor.test(ndvi.anom, swi.2yr.avg.anom, method = 'spearman')$p.val,
                                    r.dt = cor.test(ndvi.anom.dt, swi.anom.dt, method = 'spearman')$estimate,
                                    p.dt = cor.test(ndvi.anom.dt, swi.anom.dt, method = 'spearman')$p.val,
                                    r.2yr.dt = cor.test(ndvi.anom.dt, swi.2yr.avg.anom.dt, method = 'spearman')$estimate,
                                    p.2yr.dt = cor.test(ndvi.anom.dt, swi.2yr.avg.anom.dt, method = 'spearman')$p.val), 
                                by = c('site','lat.zone','clim.pxl','period')]

# recast into long format with each SWI variable as a row 
site.cors.dt <- melt.data.table(swi.anom.site.cor.dt, id.vars = c('site','lat.zone','clim.pxl','period'), measure.vars = c('r','r.2yr','r.dt','r.2yr.dt'), variable.name = 'swi.var', value.name = 'r')
site.cors.p.dt <- melt.data.table(swi.anom.site.cor.dt, id.vars = c('site','lat.zone','clim.pxl','period'), measure.vars = c('p','p.2yr','p.dt','p.2yr.dt'), value.name = 'p')
site.cors.dt$p <- site.cors.p.dt$p

site.cors.dt[swi.var == 'r', swi.var := 'swi']
site.cors.dt[swi.var == 'r.dt', swi.var := 'swi.dt']
site.cors.dt[swi.var == 'r.2yr', swi.var := 'swi.2yr']
site.cors.dt[swi.var == 'r.2yr.dt', swi.var := 'swi.2yr.dt']

# categorize trends
site.cors.dt[, sig := cut(p, c(-Inf, 0.0500, 0.1000, Inf), c('sig.p5','sig.p10','insig'))]
site.cors.dt[, r.dir := cut(r, c(-Inf, 0, Inf), c('negative','positive'))]
site.cors.dt[, cor.cat := paste(r.dir, sig, sep = '.')]
site.cors.dt[grep('insig', cor.cat), cor.cat := 'insig']

# write out correlations
site.cors.dt$rep <- i

mkdirs('output/lsat_site_swi_cors/mc_reps')
fwrite(site.cors.dt, paste0('output/lsat_site_swi_cors/mc_reps/lsat_ndvi_swi_correlations_by_site_rep_',i,'.csv'))

print('finished site correlations')


# COMPUTE PROPORTION OF SITES IN EACH BIOCLIMATIC ZONE THAT HAVE SPECIFIC NDVI - SWI CORRELATION CATEGORIES ====================================================

# zones
cor.freq.by.zone.dt <- site.cors.dt[ cor.cat != 'NA.NA' | lat.zone != '', .(n.sites = .N), by = c('swi.var', 'period', 'lat.zone', 'cor.cat')]
cor.freq.by.zone.dt <- cor.freq.by.zone.dt[, n.sites.zone := sum(n.sites), by = c('swi.var', 'period', 'lat.zone')][, pcnt.sites := n.sites / n.sites.zone * 100]
cor.freq.by.zone.dt <- cor.freq.by.zone.dt[, cor.cat := factor(cor.cat, levels = c('negative.sig.p5','negative.sig.p10','insig','positive.sig.p10','positive.sig.p5'),
                                                               labels = c('negative (p<0.05)','negative (p<0.10)','no cor.','positive (p<0.10)','positive (p<0.05)'))]
cor.freq.by.zone.dt <- cor.freq.by.zone.dt[order(swi.var, period, lat.zone, cor.cat)]  
cor.freq.by.zone.dt <- cor.freq.by.zone.dt[, pcnt.position := cumsum(pcnt.sites)-pcnt.sites/2, by = c('swi.var', 'period', 'lat.zone')]
cor.freq.by.zone.dt <- cor.freq.by.zone.dt[, cnt.position := cumsum(n.sites)-n.sites/2, by = c('swi.var', 'period', 'lat.zone')]

# biome
cor.freq.by.biome.dt <- site.cors.dt[ cor.cat != 'NA.NA', .(n.sites = .N), by = c('swi.var', 'period', 'cor.cat')]
cor.freq.by.biome.dt <- cor.freq.by.biome.dt[, n.sites.zone := sum(n.sites), by = c('swi.var', 'period')][, pcnt.sites := n.sites / n.sites.zone * 100]
cor.freq.by.biome.dt <- cor.freq.by.biome.dt[, cor.cat := factor(cor.cat, levels = c('negative.sig.p5','negative.sig.p10','insig','positive.sig.p10','positive.sig.p5'),
                                                                 labels = c('negative (p<0.05)','negative (p<0.10)','no cor.','positive (p<0.10)','positive (p<0.05)'))]
cor.freq.by.biome.dt <- cor.freq.by.biome.dt[order(swi.var, period, cor.cat)]  
cor.freq.by.biome.dt <- cor.freq.by.biome.dt[, pcnt.position := cumsum(pcnt.sites)-pcnt.sites/2, by = c('swi.var', 'period')]
cor.freq.by.biome.dt <- cor.freq.by.biome.dt[, cnt.position := cumsum(n.sites)-n.sites/2, by = c('swi.var', 'period')]
cor.freq.by.biome.dt <- cor.freq.by.biome.dt[, lat.zone := 'Arctic']

# combine zones and biome
cor.freq.by.zone.dt <- data.table(rbind(cor.freq.by.zone.dt, cor.freq.by.biome.dt))

# write out
cor.freq.by.zone.dt$rep <- i
mkdirs('output/lsat_zonal_freq_swi_cors/mc_reps/')
fwrite(cor.freq.by.zone.dt, paste0('output/lsat_zonal_freq_swi_cors/mc_reps/lsat_ndvi_swi_correlations_pcnt_cat_by_zone_rep_',i, '.csv'))

print('finished cor freq summary')  


# COMPUTE AND MAP OUT MEAN NDVI - SWI CORRELATION FOR EACH CLIMATE GRID CELL =======================================================================

mkdirs('output/lsat_gridded_swi_cors/mc_reps/')

arctic.zones.r <- raster('data/gis_data/arctic_zones/arctic_oroarctic_lat_zones_laea_50km.tif')
arctic.tmplt.r <- arctic.zones.r

periods <- unique(site.cors.dt$period)
swi.vars <- unique(site.cors.dt$swi.var)

for (j in swi.vars){
  for (k in periods){
    
    # select cols
    rep.dt <- site.cors.dt[swi.var == j & period == k]
    
    # correlate ndvi and swi 
    pxl.cor.dt <- rep.dt[, .(r = mean(r, na.rm=T)), by = 'clim.pxl']
    
    # flush template raster with NAs
    arctic.tmplt.r[arctic.tmplt.r >= 0] <- NA
    
    # map mean correlations
    cor.avg.r <- arctic.tmplt.r
    cor.avg.r[pxl.cor.dt$clim.pxl] <- pxl.cor.dt$r
    
    # write out correlation rasters
    #outname <- paste0('output/lsat_gridded_swi_cors/mc_reps/arctic_lsat_ndvi_', gsub('\\.','_',cor.types.dt$swi.var[j]), '_cor_', gsub('-','to',k),'_rep_',i,'.tif')
    outname <- paste0('output/lsat_gridded_swi_cors/mc_reps/arctic_mean_lsat_ndvi_', gsub('\\.','_',j), '_cor_', gsub('-','to',k),'_rep_',i,'.tif')
    writeRaster(cor.avg.r, outname, overwrite=T)
  }
}

print('finished gridding correlations')  


# COMPUTE ZONAL NDVI - SWI CORRELATIONS ======================================================================================
# mean climate data by pixel for locations with NDVI data
clim.pxl.dt <- site.dt[, .(swi = mean(swi, na.rm=T), swi.2yr.avg = mean(swi.2yr.avg, na.rm=T)), by = c('clim.pxl','lat.zone','year','period')]
clim.pxl.dt <- clim.pxl.dt[, swi.anom := scale(swi, center = T, scale = F), by = c('clim.pxl','lat.zone','period')]
clim.pxl.dt <- clim.pxl.dt[, swi.2yr.avg.anom := scale(swi.2yr.avg, center = T, scale = F), by = c('clim.pxl','lat.zone','period')]
clim.pxl.dt <- clim.pxl.dt[, swi.anom.dt := lm(swi.anom ~ year)$resid, by = c('clim.pxl','lat.zone','period')]
clim.pxl.dt <- clim.pxl.dt[, swi.2yr.avg.anom.dt := lm(swi.2yr.avg.anom ~ year)$resid, by = c('clim.pxl','lat.zone','period')]

swi.latzone.dt <- clim.pxl.dt[, .(swi.anom = mean(swi.anom, na.rm=T), 
                                  swi.anom.dt = mean(swi.anom.dt, na.rm=T),
                                  swi.2yr.avg.anom = mean(swi.2yr.avg.anom, na.rm=T), 
                                  swi.2yr.avg.anom.dt = mean(swi.2yr.avg.anom.dt, na.rm=T)), 
                              by = c('lat.zone', 'year','period')][order(lat.zone,year)]
swi.biome.dt <- clim.pxl.dt[, .(swi.anom = mean(swi.anom, na.rm=T), 
                                swi.anom.dt = mean(swi.anom.dt, na.rm=T),
                                swi.2yr.avg.anom = mean(swi.2yr.avg.anom, na.rm=T), 
                                swi.2yr.avg.anom.dt = mean(swi.2yr.avg.anom.dt, na.rm=T)), 
                            by = c('year','period')][order(year)][,lat.zone:='Arctic']
swi.zonal.dt <- rbind(swi.latzone.dt, swi.biome.dt)

latzone.dt <- site.dt[, .(ndvi.anom = mean(ndvi.anom, na.rm=T), ndvi.anom.dt = mean(ndvi.anom.dt, na.rm=T)), by = c('lat.zone', 'year','period')][order(lat.zone,year)]
biome.dt <- site.dt[, .(ndvi.anom = mean(ndvi.anom, na.rm=T), ndvi.anom.dt = mean(ndvi.anom.dt, na.rm=T)), by = c('year','period')][order(year)][,lat.zone:='Arctic']
zonal.dt <- rbind(latzone.dt, biome.dt)

zonal.dt <- swi.zonal.dt[zonal.dt, on = c('lat.zone','year','period')]
zonal.dt <- zonal.dt[order(lat.zone,period,year)]

# write out zonal avg time series
mkdirs('output/lsat_zonal_avg_timeseries/mc_reps')
zonal.dt$rep <- i
fwrite(zonal.dt, paste0('output/lsat_zonal_avg_timeseries/mc_reps/lsat_zonal_avg_timeseries_rep_',i,'.csv'))

# compute correlation
zonal.cor.dt <- zonal.dt[, .(r = round(cor.test(ndvi.anom, swi.anom, method = 'spearman')$estimate,2),
                             p = round(cor.test(ndvi.anom, swi.anom, method = 'spearman')$p.val,5),
                             r.dt = round(cor.test(ndvi.anom.dt, swi.anom.dt, method = 'spearman')$estimate,2),
                             p.dt = round(cor.test(ndvi.anom.dt, swi.anom.dt, method = 'spearman')$p.val,5),
                             r.2yr = round(cor.test(ndvi.anom, swi.2yr.avg.anom, method = 'spearman')$estimate,2),
                             p.2yr = round(cor.test(ndvi.anom, swi.2yr.avg.anom, method = 'spearman')$p.val,5),
                             r.2yr.dt = round(cor.test(ndvi.anom.dt, swi.2yr.avg.anom.dt, method = 'spearman')$estimate,2),
                             p.2yr.dt = round(cor.test(ndvi.anom.dt, swi.2yr.avg.anom.dt, method = 'spearman')$p.val,5)), 
                         by=c('lat.zone','period')]

# write out table of correlation coefficients
mkdirs('output/lsat_zonal_avg_swi_cors/mc_reps')
zonal.cor.dt$rep <- i
fwrite(zonal.cor.dt, paste0('output/lsat_zonal_avg_swi_cors/mc_reps/Landsat_NDVI_clim_cors_zonal_rep_',i,'.csv'))

print('finished zonal correlations... All done!')  
# write out zonal time series
# mkdirs('output/swi_zonal_avg_anom_timeseries/mc_reps')
# fwrite(zonal.dt, paste0('output/swi_zonal_avg_anom_timeseries/mc_reps/zonal_mean_swi_ndvi_anom_timeseries_1985to2016_rep_',i,'.csv'))

# quick figure
# ggplot(zonal.dt, aes(swi.anom, ndvi.anom)) + geom_point() + facet_grid(lat.zone ~ period, scales = 'free')
# ggplot(zonal.dt, aes(swi.anom.dt, ndvi.anom.dt)) + geom_point() + facet_grid(lat.zone ~ period, scales = 'free_y')
# ggplot(zonal.dt, aes(swi.2yr.avg.anom.dt, ndvi.anom.dt)) + geom_point() + facet_grid(lat.zone ~ period, scales = 'free_y')

# END SCRIPT ==================================================================================================

# cols <- c('site','year','clim.pxl',cor.types.dt$ndvi.var[j], cor.types.dt$swi.var[j])
# rep.dt <- site.dt[period == k, ..cols]
# names(rep.dt) <- c('site','year','clim.pxl','ndvi','swi')
# # summarize ndvi and swi across clim pixels
# pxl.avg.dt <- rep.dt[, ":="(ndvi.avg = mean(ndvi, na.rm=T), swi.avg = mean(swi, na.rm=T)), by = c('clim.pxl','year')]
# # correlate ndvi and swi 
# pxl.cor.dt <- pxl.avg.dt[, .(r = cor.test(ndvi.avg, swi.avg, method = 'spearman')$estimate), by = 'clim.pxl']