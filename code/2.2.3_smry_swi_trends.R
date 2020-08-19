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

setwd('/projects/arctic/users/lberner/arctic_greening/')

# READ IN FILES
zonal.avg.trnds.dt <- do.call("rbind", lapply(list.files('output/swi_zonal_trends/mc_reps/', full.names = T), fread))
site.trnds.dt <- do.call("rbind", lapply(list.files('output/swi_site_trends/mc_reps/', full.names = T), fread))

trnd.gte1985.raster.files <- list.files('output/swi_gridded_trends/mc_reps/1985to2016/', full.names = T)
trnd.gte2000.raster.files <- list.files('output/swi_gridded_trends/mc_reps/2000to2016/', full.names = T)

# SITE TRENDS -------------------------------------------------------------------
site.trnds.smry.dt <- site.trnds.dt[, .(
  swi.slope=median(slope, na.rm = T), swi.slope.q025=quantile(slope,0.025, na.rm = T), swi.slope.q975=quantile(slope,0.975, na.rm = T),
  swi.int=median(int, na.rm = T), swi.int.q025=quantile(int,0.025, na.rm = T), swi.int.q975=quantile(int,0.975, na.rm = T),
  swi.change=median(total.change, na.rm = T), swi.change.q025=quantile(total.change,0.025, na.rm = T), swi.change.q975=quantile(total.change,0.975, na.rm = T),
                                        swi.change.pcnt=median(total.change.pcnt, na.rm = T), swi.change.pcnt.q025=quantile(total.change.pcnt,0.025, na.rm = T), swi.change.pcnt.q975=quantile(total.change.pcnt,0.975, na.rm = T)),
                                    by = c('site','period')]

fwrite(site.trnds.smry.dt, 'output/swi_site_trends/lsat_swi_site_trend_summary.csv')


# TRENDS IN ZONAL AVERAGE SWI -------------------------------------------------------------------
zonal.avg.trnds.smry.dt <- zonal.avg.trnds.dt[, .(#int=median(int), int.q025=quantile(int,0.025), int.q975=quantile(int,0.975),
                                                  slope=median(slope), slope.q025=quantile(slope,0.025), slope.q975=quantile(slope,0.975),
                                                  tau=median(tau), tau.q025=quantile(tau,0.025), tau.q975=quantile(tau,0.975),
                                                  pval=median(pval), pval.q025=quantile(pval,0.025), pval.q975=quantile(pval,0.975),
                                                  total.change=median(total.change), total.change.q025=quantile(total.change,0.025), total.change.q975=quantile(total.change,0.975),
                                                  total.change.pcnt=median(total.change.pcnt), total.change.pcnt.q025=quantile(total.change.pcnt,0.025), total.change.pcnt.q975=quantile(total.change.pcnt,0.975)),
                                              by = c('period','lat.zone')]

zonal.avg.trnds.smry.dt
fwrite(zonal.avg.trnds.smry.dt, 'output/swi_zonal_trends/swi_zonal_avg_trend_summary.csv')

# pretty output table 
zonal.avg.trnds.smry.fancy.dt <- zonal.avg.trnds.smry.dt[, .(tau = paste0(sprintf("%.2f", tau),' [',sprintf("%.2f", tau.q025),',',sprintf("%.2f", tau.q975),']'),
                                                             total.change = paste0(sprintf("%.1f", total.change),' [',sprintf("%.1f", total.change.q025),',',sprintf("%.1f", total.change.q975),']'),
                                                             total.change.pcnt = paste0(sprintf("%.1f", total.change.pcnt),' [',sprintf("%.1f", total.change.pcnt.q025),',',sprintf("%.1f", total.change.pcnt.q975),']'),
                                                             slope = paste0(sprintf("%.2f", slope),' [',sprintf("%.2f", slope.q025),',',sprintf("%.2f", slope.q975),']')), 
                                                         by=c('period','lat.zone')]

zonal.avg.trnds.smry.fancy.dt[, lat.zone := factor(lat.zone, levels = c('Arctic','High Arctic','Low Arctic','Oro Arctic'))]
zonal.avg.trnds.smry.fancy.dt <- zonal.avg.trnds.smry.fancy.dt[order(period, lat.zone)]

zonal.avg.trnds.smry.fancy.dt
fwrite(zonal.avg.trnds.smry.fancy.dt, 'output/swi_zonal_trends/swi_zonal_avg_trend_summary_fancy.csv')


# GRIDDED SWI TRENDS ------------------------------------------------------------

# trends 1985 to 2016 ---
print("Starting gridded trends 1985-2016...")
stk <- stack(trnd.gte1985.raster.files, bands = 1)

stk.med <- calc(stk, fun=median)
plot(stk.med)
writeRaster(stk.med, 'output/swi_gridded_trends/arctic_swi_ensemble_trend_1985to2016_med_50km_laea.tif', overwrite=T)

stk.p025 <- calc(stk, fun=function(x){quantile(x, 0.025, na.rm = T)})
writeRaster(stk.p025, 'output/swi_gridded_trends/arctic_swi_ensemble_trend_1985to2016_p025_50km_laea.tif', overwrite=T)

stk.p975 <- calc(stk, fun=function(x){quantile(x, 0.975, na.rm = T)})
writeRaster(stk.p975, 'output/swi_gridded_trends/arctic_swi_ensemble_trend_1985to2016_p975_50km_laea.tif', overwrite=T)

rm(stk)

# trends 2000 to 2016 ---
print("Starting gridded trends 2000-2016...")
stk <- stack(trnd.gte2000.raster.files, bands = 1)

stk.med <- calc(stk, fun=median)
writeRaster(stk.med, 'output/swi_gridded_trends/arctic_swi_ensemble_trend_2000to2016_med_50km_laea.tif', overwrite=T)

stk.p025 <- calc(stk, fun=function(x){quantile(x, 0.025, na.rm = T)})
writeRaster(stk.p025, 'output/swi_gridded_trends/arctic_swi_ensemble_trend_2000to2016_p025_50km_laea.tif', overwrite=T)

stk.p975 <- calc(stk, fun=function(x){quantile(x, 0.975, na.rm = T)})
writeRaster(stk.p975, 'output/swi_gridded_trends/arctic_swi_ensemble_trend_2000to2016_p975_50km_laea.tif', overwrite=T)

rm(stk)

print("All done!")
# END SCRIPTS --------------------------------------
