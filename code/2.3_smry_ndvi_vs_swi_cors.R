# ABOUT THIS SCRIPT  ==============================================================================================================
# This R script summarizes NDVI - SWI correlations for sites to zones across the Arctic 
# AUTHOR: LOGAN BERNER, NAU
# DATE: 2020-04-27

# SET UP WORKSPACE ==============================================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(data.table)
require(dplyr)
require(ggplot2)
require(raster)

setwd('/projects/arctic/users/lberner/arctic_greening/')

# READ IN FILES
site.cors.dt <- do.call("rbind", lapply(list.files('output/lsat_site_swi_cors/mc_reps/', full.names = T), fread))
zonal.avg.cors.dt <- do.call("rbind", lapply(list.files('output/lsat_zonal_avg_swi_cors/mc_reps/', full.names = T), fread))
zonal.freq.cors.dt <- do.call("rbind", lapply(list.files('output/lsat_zonal_freq_swi_cors/mc_reps/', full.names = T), fread))
cor.raster.files <- list.files('output/lsat_gridded_swi_cors/mc_reps/', full.names = T)


# MEAN CORRELATION AMONG SITES SPLIT SEVERAL WAYS -------------------------------------------------
site.cors.dt[cor.cat == 'positive.sig.p5' | cor.cat == 'positive.sig.p10', cor.cat := 'positive']
site.cors.dt[cor.cat == 'negative.sig.p5' | cor.cat == 'negative.sig.p10', cor.cat := 'negative']

# mean correlation by SWI variable and NDVI-SWI correlation category
site.cor.avg.dt <- site.cors.dt[, .(r=median(r, na.rm=T)), by = c('period','swi.var','cor.cat','site')][, .(r.avg=mean(r, na.rm=T), r.sd=sd(r,na.rm=T)), by = c('period','swi.var','cor.cat')]
fwrite(site.cor.avg.dt, 'output/lsat_site_swi_cors/lsat_ndvi_swi_mean_cor_by_cor_swi_cats.csv')

# mean correlation by SWI variable
site.cor.avg.dt <- site.cors.dt[, .(r=median(r, na.rm=T)), by = c('period','swi.var','cor.cat','site')][, .(r.avg=mean(r, na.rm=T), r.sd=sd(r,na.rm=T)), by = c('period','swi.var')]
fwrite(site.cor.avg.dt, 'output/lsat_site_swi_cors/lsat_ndvi_swi_mean_cor_by_swi_cats.csv')


# CORRELATION BETWEEN ZONAL AVERAGE NDVI AND SWI -------------------------------------------------------------------
zonal.avg.cors.smry.dt <- zonal.avg.cors.dt[, .(r=median(r), r.q005=quantile(r,0.005), r.q025=quantile(r,0.025), r.q975=quantile(r,0.975),
                                                r.dt=median(r.dt), r.dt.q005=quantile(r.dt,0.005), r.dt.q025=quantile(r.dt,0.025), r.dt.q975=quantile(r.dt,0.975),
                                                r.2yr=median(r.2yr), r.2yr.q025=quantile(r.2yr,0.025), r.2yr.q975=quantile(r.2yr,0.975),
                                                r.2yr.dt=median(r.2yr.dt), r.2yr.dt.q005=quantile(r.2yr.dt,0.005), r.2yr.dt.q025=quantile(r.2yr.dt,0.025), r.2yr.dt.q975=quantile(r.2yr.dt,0.975)),
                                             by = c('period','lat.zone')]

# write out
zonal.avg.cors.smry.dt
fwrite(zonal.avg.cors.smry.dt, 'output/lsat_zonal_avg_swi_cors/lsat_zonal_avg_swi_cors_summary.csv')

# fancify table
zonal.avg.cors.smry.fancy.dt <- zonal.avg.cors.smry.dt[, .(r = paste0(sprintf('%.2f', r),' [', sprintf('%.2f', r.q025),',', sprintf('%.2f', r.q975),']'),
                                                           r.dt = paste0(sprintf('%.2f', r.dt),' [', sprintf('%.2f', r.dt.q025),',', sprintf('%.2f', r.dt.q975),']'),
                                                           r.2yr = paste0(sprintf('%.2f', r.2yr),' [', sprintf('%.2f', r.2yr.q025),',',sprintf('%.2f', r.2yr.q975),']'),
                                                           r.2yr.dt = paste0(sprintf('%.2f', r.2yr.dt),' [',sprintf('%.2f', r.2yr.dt.q025),',',sprintf('%.2f', r.2yr.dt.q975),']')), 
                                                       by=c('period','lat.zone')]

zonal.avg.cors.smry.fancy.dt[, lat.zone := recode(zonal.avg.cors.smry.fancy.dt$lat.zone, High = "High Arctic", Low = "Low Arctic", Oro = "Oro Arctic", biome = "Arctic")]
zonal.avg.cors.smry.fancy.dt[, lat.zone := factor(zonal.avg.cors.smry.fancy.dt$lat.zone, levels = c('Arctic','High Arctic','Low Arctic','Oro Arctic'))]
zonal.avg.cors.smry.fancy.dt <- zonal.avg.cors.smry.fancy.dt[order(period, lat.zone)]

# write out fancy table
zonal.avg.cors.smry.fancy.dt
fwrite(zonal.avg.cors.smry.fancy.dt, 'output/lsat_zonal_avg_swi_cors/lsat_zonal_avg_swi_cors_summary_fancy.csv')


# FREQ OF NDVI TRENDS BY ZONE ------------------------------------------------------------------------
zonal.freq.cors.dt <- zonal.freq.cors.dt[cor.cat != '']
zonal.freq.cors.dt[cor.cat == 'positive (p<0.05', cor.cat := 'positive (p<0.05)'] # fix typo
zonal.freq.cor.smry.dt <- zonal.freq.cors.dt[, .(n.sites=as.numeric(median(n.sites)), n.sites.q025=quantile(n.sites,0.025), n.sites.q975=quantile(n.sites,0.975),
                                                    pcnt.sites=median(pcnt.sites), pcnt.sites.q025=quantile(pcnt.sites,0.025), pcnt.sites.q975=quantile(pcnt.sites,0.975)), 
                                                by = c('swi.var','period','lat.zone','cor.cat')]

zonal.freq.cor.smry.dt
fwrite(zonal.freq.cor.smry.dt, 'output/lsat_zonal_freq_swi_cors/lsat_zonal_freq_swi_cors_summary.csv')

# aggregate trends with P < 0.1 and P < 0.05 and recompute freq
zonal.freq.cors.dt[, cor.cat.p10 := cor.cat]
zonal.freq.cors.dt[cor.cat.p10 == 'negative (p<0.05)' | cor.cat.p10 == 'negative (p<0.10)', cor.cat.p10 := 'negative']
zonal.freq.cors.dt[cor.cat.p10 == 'positive (p<0.05)' | cor.cat.p10 == 'positive (p<0.10)', cor.cat.p10 := 'positive']


zonal.freq.cors.p10.dt <- zonal.freq.cors.dt[, .(n.sites = sum(n.sites), n.sites.zone=first(n.sites.zone)), by = c('swi.var','period','lat.zone','cor.cat.p10','rep')]
zonal.freq.cors.p10.dt <- zonal.freq.cors.p10.dt[, pcnt.sites := n.sites/n.sites.zone*100]

zonal.freq.cors.p10.smry.dt <- zonal.freq.cors.p10.dt[, .(n.sites=median(n.sites), n.sites.q025=quantile(n.sites,0.025), n.sites.q975=quantile(n.sites,0.975),
                                                            pcnt.sites=median(pcnt.sites), pcnt.sites.q025=quantile(pcnt.sites,0.025), pcnt.sites.q975=quantile(pcnt.sites,0.975)), 
                                                        by = c('swi.var','period','lat.zone','cor.cat.p10')]
zonal.freq.cors.p10.smry.dt

fwrite(zonal.freq.cors.p10.smry.dt, 'output/lsat_zonal_freq_swi_cors/lsat_zonal_freq_swi_cors_p10_summary.csv')



# GRIDDED NDVI TRENDS ------------------------------------------------------------
# swi.vars <- c('swi_cor','swi_dt_cor','swi_2yr_cor','swi_2yr_dt_cor')
swi.vars <- c('swi_dt_cor','swi_2yr_cor','swi_2yr_dt_cor')
periods <- c('1985to2016','2000to2016')

for (j in swi.vars){
  for (k in periods){
    # identify files and load raster stack
    cor.files <- grep(k, grep(j, cor.raster.files, value = T), value = T)
    stk <- stack(cor.files)
    
    # summarize across raster stacks
    stk.med <- calc(stk, fun=median)
    writeRaster(stk.med, paste0('output/lsat_gridded_swi_cors/arctic_lsat_ndvi_', j, '_', k, '_med_50km_laea.tif'), overwrite=T)
    
#     stk.p025 <- calc(stk, fun=function(x){quantile(x, 0.025, na.rm = T)})
#     writeRaster(stk.p025, paste0('output/lsat_gridded_swi_cors/arctic_lsat_ndvi_', j, '_', k, '_p025_50km_laea.tif'), overwrite=T)
#     
#     stk.p975 <- calc(stk, fun=function(x){quantile(x, 0.975, na.rm = T)})
#     writeRaster(stk.p975, paste0('output/lsat_gridded_swi_cors/arctic_lsat_ndvi_', j, '_', k, '_p975_50km_laea.tif'), overwrite=T)
#     
    # status report
    print(paste0('Finished: ', j, ' from ', k))
  }
}

# QUICK FIGURES
# ggplot(zonal.avg.cors.dt, aes(dataset, r)) + geom_boxplot() + facet_grid(period ~ lat.zone)
# ggplot(zonal.avg.cors.dt, aes(dataset, r.2yr.dt)) + geom_boxplot() + facet_grid(period ~ lat.zone)
# ggplot(zonal.avg.cors.dt, aes(r)) + geom_histogram() + facet_grid(period ~ lat.zone)

# END SCRIPTS --------------------------------------


## HOLDING 
# zonal.avg.cors.smry.dt <- zonal.avg.cors.dt[, .(r=median(r), r.q025=quantile(r,0.025), r.q975=quantile(r,0.975),
#                                                 p=median(p), p.q95=quantile(p,0.95),
#                                                 r.dt=median(r.dt), r.dt.q025=quantile(r.dt,0.025), r.dt.q975=quantile(r.dt,0.975),
#                                                 p.dt=median(p.dt), p.dt.q95=quantile(p.dt,0.95),
#                                                 r.2yr=median(r.2yr), r.yr.q025=quantile(r.2yr,0.025), r.2yr.q975=quantile(r.2yr,0.975),
#                                                 p.2yr=median(p.2yr), p.q95=quantile(p.2yr,0.95),
#                                                 r.2yr.dt=median(r.2yr.dt), r.2yr.dt.q025=quantile(r.2yr.dt,0.025), r.2yr.dt.q975=quantile(r.2yr.dt,0.975),
#                                                 p.2yr.dt=median(p.2yr.dt), p.2yr.dt.q95=quantile(p.2yr.dt,0.95)),
#                                             by = c('period','lat.zone')]
