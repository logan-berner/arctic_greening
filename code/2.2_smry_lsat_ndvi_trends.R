# ABOUT THIS SCRIPT  ==============================================================================================================
# This R script summarizes Landsat cross-sensor calibration random forest models
# AUTHOR: LOGAN BERNER, NAU
# DATE: 2020-02-25

# SET UP WORKSPACE ==============================================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(data.table)
require(ggplot2)
require(raster)
require(dplyr)

setwd('/projects/arctic/users/lberner/arctic_greening/')

# READ IN FILES
site.trnds.dt <- do.call("rbind", lapply(list.files('output/lsat_site_trends/mc_reps/', full.names = T), fread))
zonal.avg.trnds.dt <- do.call("rbind", lapply(list.files('output/lsat_zonal_avg_trends/mc_reps/', full.names = T), fread))
zonal.freq.trnds.dt <- do.call("rbind", lapply(list.files('output/lsat_zonal_freq_trends/mc_reps/', full.names = T), fread))
grn.gte1985.raster.files <- list.files('output/lsat_gridded_trends/mc_reps/', pattern = 'greening_p10_1985to2016', full.names = T)
grn.gte2000.raster.files <- list.files('output/lsat_gridded_trends/mc_reps/', pattern = 'greening_p10_2000to2016', full.names = T)
brn.gte1985.raster.files <- list.files('output/lsat_gridded_trends/mc_reps/', pattern = 'browning_p10_1985to2016', full.names = T)
brn.gte2000.raster.files <- list.files('output/lsat_gridded_trends/mc_reps/', pattern = 'browning_p10_2000to2016', full.names = T)
site.cnt.gte1985.raster.files <- list.files('output/lsat_gridded_trends/mc_reps/', pattern = 'site_cnt_1985to2016', full.names = T)
site.cnt.gte2000.raster.files <- list.files('output/lsat_gridded_trends/mc_reps/', pattern = 'site_cnt_2000to2016', full.names = T)

# SITE TRENDS NDVI -------------------------------------------------------------------
site.trnds.smry.dt <- site.trnds.dt[, .(int=median(int, na.rm = T), int.q025=quantile(int,0.025, na.rm = T), int.q975=quantile(int,0.975, na.rm = T),
                                        ndvi.change=median(total.change, na.rm = T), ndvi.change.q025=quantile(total.change,0.025, na.rm = T), ndvi.change.q975=quantile(total.change,0.975, na.rm = T),
                                        ndvi.change.pcnt=median(total.change.pcnt, na.rm = T), ndvi.change.pcnt.q025=quantile(total.change.pcnt,0.025, na.rm = T), ndvi.change.pcnt.q975=quantile(total.change.pcnt,0.975, na.rm = T)),
                                    by = c('site','period')]

fwrite(site.trnds.smry.dt, 'output/lsat_site_trends/lsat_ndvi_site_trend_summary.csv')


# TRENDS IN ZONAL AVERAGE NDVI -------------------------------------------------------------------
zonal.avg.trnds.smry.dt <- zonal.avg.trnds.dt[, .(#int=median(int), int.q025=quantile(int,0.025), int.q975=quantile(int,0.975),
                                                  slope=median(slope), slope.q025=quantile(slope,0.025), slope.q975=quantile(slope,0.975),
                                                  tau=median(tau), tau.q025=quantile(tau,0.025), tau.q975=quantile(tau,0.975),
                                                  pval=median(pval), pval.q025=quantile(pval,0.025), pval.q975=quantile(pval,0.975),
                                                  total.change=median(total.change), total.change.q025=quantile(total.change,0.025), total.change.q975=quantile(total.change,0.975),
                                                  total.change.pcnt=median(total.change.pcnt), total.change.pcnt.q025=quantile(total.change.pcnt,0.025), total.change.pcnt.q975=quantile(total.change.pcnt,0.975)),
                                              by = c('period','lat.zone')]

zonal.avg.trnds.smry.dt
fwrite(zonal.avg.trnds.smry.dt, 'output/lsat_zonal_avg_trends/lsat_zonal_avg_trend_summary.csv')

# pretty output table 
zonal.avg.trnds.smry.fancy.dt <- zonal.avg.trnds.smry.dt[, .(tau = paste0(sprintf("%.2f", tau),' [',sprintf("%.2f", tau.q025),',',sprintf("%.2f", tau.q975),']'),
                                                             total.change = paste0(sprintf("%.3f", total.change),' [',sprintf("%.3f", total.change.q025),',',sprintf("%.3f", total.change.q975),']'),
                                                             total.change.pcnt = paste0(sprintf("%.1f", total.change.pcnt),' [',sprintf("%.1f", total.change.pcnt.q025),',',sprintf("%.1f", total.change.pcnt.q975),']'),
                                                             slope = paste0(sprintf("%.4f", slope),' [',sprintf("%.4f", slope.q025),',',sprintf("%.4f", slope.q975),']')), 
                                                         by=c('period','lat.zone')]

zonal.avg.trnds.smry.fancy.dt[, lat.zone := recode(lat.zone, High = "High Arctic", Low = "Low Arctic", Oro = "Oro Arctic", biome = "Arctic")]
zonal.avg.trnds.smry.fancy.dt[, lat.zone := factor(lat.zone, levels = c('Arctic','High Arctic','Low Arctic','Oro Arctic'))]
zonal.avg.trnds.smry.fancy.dt <- zonal.avg.trnds.smry.fancy.dt[order(period, lat.zone)]

zonal.avg.trnds.smry.fancy.dt
fwrite(zonal.avg.trnds.smry.fancy.dt, 'output/lsat_zonal_avg_trends/lsat_zonal_avg_trend_summary_fancy.csv')


# FREQ OF NDVI TRENDS BY ZONE ------------------------------------------------------------------------
zonal.freq.trnds.smry.dt <- zonal.freq.trnds.dt[, .(n.sites=median(n.sites), n.sites.q025=quantile(n.sites,0.025), n.sites.q975=quantile(n.sites,0.975),
                                                    pcnt.sites=median(pcnt.sites), pcnt.sites.q025=quantile(pcnt.sites,0.025), pcnt.sites.q975=quantile(pcnt.sites,0.975)), 
                                                by = c('period','lat.zone','trend.cat')]

zonal.freq.trnds.smry.dt
fwrite(zonal.freq.trnds.smry.dt, 'output/lsat_zonal_freq_trends/lsat_zonal_freq_trend_summary.csv')

# aggregate trends with P < 0.1 and P < 0.05 and recompute freq
zonal.freq.trnds.dt$slope.cat <- NA_character_
zonal.freq.trnds.dt[trend.cat == 'browning (p<0.05)' | trend.cat == 'browning (p<0.10)', slope.cat := 'browning']
zonal.freq.trnds.dt[trend.cat == 'greening (p<0.05)' | trend.cat == 'greening (p<0.10)', slope.cat := 'greening']
zonal.freq.trnds.dt[trend.cat == 'no trend', slope.cat := 'no trend']

zonal.freq.trnds.p10.dt <- zonal.freq.trnds.dt[, .(n.sites = sum(n.sites),n.sites.zone=first(n.sites.zone)), by = c('period','lat.zone','slope.cat','rep')]
zonal.freq.trnds.p10.dt <- zonal.freq.trnds.p10.dt[, pcnt.sites := n.sites/n.sites.zone*100]

## calc ratio of greening to browning
g2b.dt <- dcast(zonal.freq.trnds.p10.dt, period + lat.zone + rep ~ slope.cat, value.var = 'pcnt.sites')
g2b.dt <- g2b.dt[, g2b := greening / browning][, c('browning','greening','no trend') := NULL]
zonal.freq.trnds.p10.dt <- g2b.dt[zonal.freq.trnds.p10.dt, on = c('period','lat.zone','rep')]

zonal.freq.trnds.p10.smry.dt <- zonal.freq.trnds.p10.dt[, .(n.sites=median(n.sites), n.sites.q025=quantile(n.sites,0.025), n.sites.q975=quantile(n.sites,0.975),
                                                    pcnt.sites=median(pcnt.sites), pcnt.sites.q025=quantile(pcnt.sites,0.025), pcnt.sites.q975=quantile(pcnt.sites,0.975), 
                                                    g2b=median(g2b), g2b.q025=quantile(g2b,0.025), g2b.q975=quantile(g2b,0.975)),
                                                by = c('period','lat.zone','slope.cat')]
zonal.freq.trnds.p10.smry.dt <- cbind(zonal.freq.trnds.p10.smry.dt[, 1:3], round(zonal.freq.trnds.p10.smry.dt[,-c(1:3)],1)) # round numeric cols

fwrite(zonal.freq.trnds.p10.smry.dt, 'output/lsat_zonal_freq_trends/lsat_zonal_freq_trend_p10_summary.csv')

# fancy wide-format table
zonal.freq.trnds.p10.smry.dt[, ':='(n.sites = sum(n.sites), n.sites.q025 = sum(n.sites.q025), n.sites.q975 = sum(n.sites.q975)), by = c('period','lat.zone')]
zonal.freq.trnds.p10.smry.fancy.dt <- zonal.freq.trnds.p10.smry.dt[, .(n.sites = paste0(sprintf("%.0f", n.sites),' [',sprintf("%.0f", n.sites.q025),', ',sprintf("%.0f", n.sites.q975),']'),
                                                                       pcnt.sites = paste0(sprintf("%.1f", pcnt.sites),' [',sprintf("%.1f", pcnt.sites.q025),', ',sprintf("%.1f", pcnt.sites.q975),']'),
                                                                       g2b = paste0(sprintf("%.1f", g2b),' [',sprintf("%.1f", g2b.q025),', ',sprintf("%.1f", g2b.q975),']')), 
                                                                   by=c('period','lat.zone','slope.cat')]

zonal.freq.trnds.p10.smry.fancy.dt <- zonal.freq.trnds.p10.smry.fancy.dt[order(period,lat.zone)]
zonal.freq.trnds.p10.smry.fancy.wide.dt <- dcast(zonal.freq.trnds.p10.smry.fancy.dt, period + lat.zone + n.sites ~ slope.cat, value.var = c('pcnt.sites','g2b'))
fwrite(zonal.freq.trnds.p10.smry.fancy.wide.dt, 'output/lsat_zonal_freq_trends/lsat_zonal_freq_trend_p10_summary_fancy_wide.csv')

# GRIDDED NDVI TRENDS ------------------------------------------------------------
print('starting gridded greening trends...')
# greening 1985 to 2016 ---
stk <- stack(grn.gte1985.raster.files)

stk.med <- calc(stk, fun=median)
writeRaster(stk.med, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_greening_p10_1985to2016_med_50km_laea.tif', overwrite=T)

stk.p025 <- calc(stk, fun=function(x){quantile(x, 0.025, na.rm = T)})
writeRaster(stk.p025, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_greening_p10_1985to2016_p025_50km_laea.tif', overwrite=T)

stk.p975 <- calc(stk, fun=function(x){quantile(x, 0.975, na.rm = T)})
writeRaster(stk.p975, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_greening_p10_1985to2016_p975_50km_laea.tif', overwrite=T)

rm(stk)

# greening 2000 to 2016 ---
stk <- stack(grn.gte2000.raster.files)

stk.med <- calc(stk, fun=median)
writeRaster(stk.med, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_greening_p10_2000to2016_med_50km_laea.tif', overwrite=T)

stk.p025 <- calc(stk, fun=function(x){quantile(x, 0.025, na.rm = T)})
writeRaster(stk.p025, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_greening_p10_2000to2016_p025_50km_laea.tif', overwrite=T)

stk.p975 <- calc(stk, fun=function(x){quantile(x, 0.975, na.rm = T)})
writeRaster(stk.p975, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_greening_p10_2000to2016_p975_50km_laea.tif', overwrite=T)

rm(stk)

# browning 1985 to 2016 --- 
print('starting gridded browning trends...')
stk <- stack(brn.gte1985.raster.files)

stk.med <- calc(stk, fun=median)
writeRaster(stk.med, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_browning_p10_1985to2016_med_50km_laea.tif', overwrite=T)

stk.p025 <- calc(stk, fun=function(x){quantile(x, 0.025, na.rm = T)})
writeRaster(stk.p025, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_browning_p10_1985to2016_p025_50km_laea.tif', overwrite=T)

stk.p975 <- calc(stk, fun=function(x){quantile(x, 0.975, na.rm = T)})
writeRaster(stk.p975, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_browning_p10_1985to2016_p975_50km_laea.tif', overwrite=T)

plot(stk.med)

rm(stk)

# browning 2000 to 2016 ---
stk <- stack(brn.gte2000.raster.files)

stk.med <- calc(stk, fun=median)
writeRaster(stk.med, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_browning_p10_2000to2016_med_50km_laea.tif', overwrite=T)

stk.p025 <- calc(stk, fun=function(x){quantile(x, 0.025, na.rm = T)})
writeRaster(stk.p025, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_browning_p10_2000to2016_p025_50km_laea.tif', overwrite=T)

stk.p975 <- calc(stk, fun=function(x){quantile(x, 0.975, na.rm = T)})
writeRaster(stk.p975, 'output/lsat_gridded_trends/arctic_lsat_ndvi_pcnt_browning_p10_2000to2016_p975_50km_laea.tif', overwrite=T)

plot(stk.med)

rm(stk)

# stack count 1985 to 2016 ---
print('starting gridded site counts...')

stk <- stack(site.cnt.gte1985.raster.files)

stk.med <- calc(stk, fun=median)
writeRaster(stk.med, 'output/lsat_gridded_trends/arctic_lsat_site_cnt_1985to2016_med_50km_laea.tif', overwrite=T)

# stack count 1985 to 2016 ---
stk <- stack(site.cnt.gte2000.raster.files)

stk.med <- calc(stk, fun=median)
writeRaster(stk.med, 'output/lsat_gridded_trends/arctic_lsat_site_cnt_2000to2016_med_50km_laea.tif', overwrite=T)

plot(stk.med)

rm(stk)

print('All done!')
# END SCRIPTS --------------------------------------
# 
# ggplot(zonal.avg.trnds.dt, aes(total.change.pcnt)) + geom_histogram() + facet_grid(lat.zone ~ period)
