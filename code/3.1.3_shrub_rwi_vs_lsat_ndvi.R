# ABOUT THIS SCRIPT ===================================================================================================
# This R script computes correlation between shrub radial growth and Landsat NDVI time series for sites across the Arctic in a Monte Carlo framework
# Author: Logan Berner, NAU
# Date: 220-04-14

# SET UP WORKSPACE ===================================================================================================
rm(list=ls())
require(data.table)
require(lattice)
require(dplyr)
setwd('/projects/arctic/users/lberner/arctic_greening/')


# LOAD DATA SETS ===================================================================================================
crn.files <- list.files('output/site_comparisons/shrub_chrons_mcreps/', full.names = T)
lsat.files <- list.files('output/lsat_site_extractions/lsat_shrub_site_mcreps/', full.names = T)
sites.dt <- fread('data/field_data/tundra_shrub_growth/tundra_shrub_cluster_site_coords_20190507.csv')

# COMPUTE CORRELATIONS BETWEEN LANDSAT NDVI AND SHRUB RWI ================================================================
mc.reps <- 1000
cor.mc.list <- list()
bio.mc.list <- list()

for (i in 1:mc.reps){
  # get data sets
  crn.dt <- fread(crn.files[i])
  lsat.dt <- fread(lsat.files[i])
  
  # merge data into one data table
  clst.dt <- lsat.dt[crn.dt, on = c('cluster','genus','year')]
  
  # filter years with too few obs
  clst.dt <- clst.dt[is.na(ndvi.med) == F]
  clst.dt <- clst.dt[, n.yrs := .N, by = c('cluster','genus')]
  clst.dt <- clst.dt[n.yrs >= 10] # require at least 10 years of overlap
  
  # bootstrap years used at each site
  # clst.dt <- clst.dt[, .SD[sample(.N,.N, replace = T)], by = c('cluster','genus')]
  clst.dt <- clst.dt[, .SD[sample(.N,.N*0.9)], by = c('cluster','genus')]
  
  # compute correlations
  cor.dt <- clst.dt[, .(r = cor.test(ndvi.dt.med, rwi.med, method = 'spearman')$estimate,
                        rbar.avg = mean(rbar.avg)),
                    by = c('cluster','genus')]
  
  # save correlations and time series lists
  cor.dt[, rep := i]
  cor.mc.list[[i]] <- cor.dt
  
  clst.dt[, rep := i]
  bio.mc.list[[i]] <- clst.dt
  
  # status 
  print(i/mc.reps)
}

cor.mc.dt <- rbindlist(cor.mc.list)
bio.mc.dt <- rbindlist(bio.mc.list)

# SUMMARIZE CORRELATIONS ACROSS MONTE CARLO ITERATIONS ===================================================================

# median correlation for ecah site 
site.cor.smry.dt <- cor.mc.dt[, .(r=median(r), r.q025=quantile(r,0.025), r.q975=quantile(r,0.975), 
                                  rbar=median(rbar.avg), rbar.q025=quantile(rbar.avg,0.025), rbar.q975=quantile(rbar.avg,0.975)), by = c('cluster','genus')]
site.cor.smry.dt[order(r)]
fwrite(sit.cor.smry.fancy.dt, 'output/site_comparisons/lsat_ndvi_cor_with_shrub_growth_smry_fancy.csv')

# fancy output table
sit.cor.smry.fancy.dt <- site.cor.smry.dt
sit.cor.smry.fancy.dt[, ':='(r = paste0(sprintf('%.2f', r)), r.95 = paste0(sprintf('%.2f', r), ' [', sprintf('%.2f', r.q025), ',', sprintf('%.2f', r.q975), ']'),
                             rbar = paste0(sprintf('%.2f', rbar), ' [', sprintf('%.2f', rbar.q025), ',', sprintf('%.2f', rbar.q975), ']'))]
sit.cor.smry.fancy.dt[, c('r.q025','r.q975','rbar.q025','rbar.q975') := NULL]
sit.cor.smry.fancy.dt$country <- sites.dt$country[match(sit.cor.smry.fancy.dt$cluster, sites.dt$cluster)]
sit.cor.smry.fancy.dt$lat <- round(sites.dt$lat[match(sit.cor.smry.fancy.dt$cluster, sites.dt$cluster)],3)
sit.cor.smry.fancy.dt$lon <- round(sites.dt$lon[match(sit.cor.smry.fancy.dt$cluster, sites.dt$cluster)],3)
sit.cor.smry.fancy.dt <- sit.cor.smry.fancy.dt %>% arrange(country,cluster,genus) %>% dplyr::select(country, cluster, lat, lon, genus, rbar, everything()) # rbar coef
sit.cor.smry.fancy.dt

fwrite(sit.cor.smry.fancy.dt, 'output/site_comparisons/lsat_ndvi_cor_with_shrub_growth_smry.csv')

# median correlation across all simulations 
cor.mc.dt[, .(r=median(r)), by = c('rep')][, .(r=median(r), r.q025=quantile(r,0.025), r.q975=quantile(r,0.975))]

# number ad % of sites with consisitent positive correalitons
sum(site.cor.smry.dt$r.q025 > 0, na.rm = T) 
sum(site.cor.smry.dt$r.q025 > 0, na.rm = T) / nrow(site.cor.smry.dt) * 100


# QUICK SUMMARY OF NDVI-SWI CORRELATION ===================================================================================================
median(site.cor.smry.dt$r)
min(site.cor.smry.dt$r)
max(site.cor.smry.dt$r)


bwplot(r ~ genus, site.cor.smry.dt)
site.cor.smry.dt %>% group_by(genus) %>% summarise(cor.med = median(r))


# PLOT HISTOGRAM OF RWI - NDVI CORRELATION COEFFICIENTS ===================================================================================================
my.pch.cex=1.5
my.cex.axis=1.5

my.ylab  <- "Number of shrub chronologies"
my.xlab <- expression("Landsat NDVI"[max]~'- shrub RWI correlation (r'[s]~')')

jpeg('figures/Landsat_NDVI_vs_shrubRWI_cor_hist.jpg', 6, 5, units = 'in', res = 400)
par.top <- par(mar=c(5,5.5,1,1))
hist(site.cor.smry.dt$r, xaxt='n', yaxt='n', xlim=c(-0.3,1), ylim=c(0,6), breaks = seq(-0.3,0.9,0.1),
     xlab = '', ylab='',cex = my.cex.axis, col = 'gray50', main='')
axis(1, seq(-0.3,1,0.2), cex.axis=my.cex.axis)
axis(2, seq(0,6,1), las=2, cex.axis=my.cex.axis, labels = T)
mtext(side = 1, line = 3.5, cex = my.pch.cex, my.xlab)
mtext(side = 2, line = 4.0, cex = my.pch.cex, my.ylab)
dev.off()

# PLOT MEAN NDVI AND RWI TIME SERIES FOR EACH SITE USING DOUBLE Y AXES ===================================================================================================
ndvi.plot <- xyplot(ndvi.dt.avg ~ year | cluster, clst.dt, type='b', ylab = 'NDVI (detrended)', xlim=c(1980,2017), scales=list(y=list(relation='free'))) #, col='green')
rwi.plot <- xyplot(rwi.avg ~ year | cluster, clst.dt, type='b', ylab = 'RWI (unitless)', xlim=c(1980,2017), scales=list(y=list(relation='free'))) #, col='brown')
ndvi.rwi.ts.plot <- doubleYScale(ndvi.plot, rwi.plot, col='bl')
ndvi.rwi.ts.plot

jpeg('figures/shrub_rwi_lsat_ndvi_timeseries.jpg', 20, 10, res=400, units = 'in')
ndvi.rwi.ts.plot
dev.off()
print(rwi.plot)


# SCATTER PLOTS OF NDVI VS RWI ===================================================================================================
jpeg('figures/lsat_ndvi_vs_shrub_rwi_by_site.jpg', 10, 8, res=400, units = 'in')
xyplot(ndvi ~ rwi | cluster, bio.smry.dt)
dev.off()


# SCATTER PLOTS OF NDVI VS RWI FOR SELECT SITES ===================================================================================================
exmpl.dt <- bio.smry.dt[cluster == 'Cherskii']

jpeg('figures/lsat_ndvi_vs_shrub_rwi_at_Cherskii.jpg', 3, 3, res=500, units = 'in')
xyplot(ndvi.dt.avg ~ rwi.avg, exmpl.dt, xlab='Shrub ring width index', ylab='Landsat NDVI (detrended)', pch = 16)
dev.off()

# END SCRIPT ===================================================================================================
