# This R script generates a table summarizing NDVI trends, SWI trends, and NDVI-SWI correlations by Arctic zone 
# 2019-11-18

rm(list=ls())
require(data.table)
setwd('C:/Users/lb968/Google Drive/research/nau/nsf_arctic/arctic_greening/')
# setwd('C:/Users/Logan/Google Drive/research/nau/nsf_arctic/arctic_greening/')

swi.trend.dt <- fread('output/swi_zonal_mean_trend_summary.csv')
ndvi.avg.trend.dt <- fread('output/Landsat_NDVImax_mean_anom_tundra_zone_avg_trends.csv')
# ndvi.anom.trend.dt <- fread('output/Landsat_NDVImax_mean_anom_tundra_zone_anom_trends.csv')
ndvi.swi.cor.dt <- fread('output/Landsat_NDVI_clim_cors_zonal.csv')


# NDVI trend =============================================
ndvi.avg.trend.dt <- ndvi.avg.trend.dt %>% arrange(period, lat.zone)

# absolute change in NDVI
# ndvi.abs.chng.dt <- ndvi.anom.trend.dt[, .(total.change, total.change.lbound, total.change.ubound)]
ndvi.abs.chng.dt <- ndvi.avg.trend.dt[, .(total.change, total.change.lbound, total.change.ubound)]
ndvi.abs.chng.dt <- round(ndvi.abs.chng.dt,3)
ndvi.abs.chng.vec <- paste(format(ndvi.abs.chng.dt$total.change, digits=3), ' (', format(ndvi.abs.chng.dt$total.change.lbound, digits=3), ', ', format(ndvi.abs.chng.dt$total.change.ubound, digits=3), ')', sep='')

# % change in NDVI
ndvi.pcnt.chng.dt <- ndvi.avg.trend.dt[, .(total.change.pcnt, total.change.lbound.pcnt, total.change.ubound.pcnt)]
ndvi.pcnt.chng.dt <- round(ndvi.pcnt.chng.dt,1)
ndvi.pcnt.chng.vec <- paste(format(ndvi.pcnt.chng.dt$total.change.pcnt, digits=2), ' (', format(ndvi.pcnt.chng.dt$total.change.lbound.pcnt, digits=1), ', ', format(ndvi.pcnt.chng.dt$total.change.ubound.pcnt, digits=2), ')', sep='')
ndvi.pcnt.chng.vec

# trend stats
# ndvi.tau = round(ndvi.anom.trend.dt$tau, 2)
# ndvi.p = round(ndvi.anom.trend.dt$pval, 3)
ndvi.tau = round(ndvi.avg.trend.dt$tau, 2)
ndvi.p = round(ndvi.avg.trend.dt$pval, 3)

ndvi.p[ndvi.p <= 0.001] <- '<0.001'
ndvi.tau.p.vec <- paste0(ndvi.tau, ' (', ndvi.p, ')')


# SWI trend ===============================================
# abs change in SWI
swi.abs.chng.dt <- round(swi.trend.dt[, .(total.change, total.change.lbound, total.change.ubound)],1)
swi.abs.chng.vec <- paste(format(swi.abs.chng.dt$total.change, digits = 2), ' (', format(swi.abs.chng.dt$total.change.lbound, digits = 2), ', ', format(swi.abs.chng.dt$total.change.ubound, digits=2), ')', sep='')

# % change in SWI
swi.pcnt.chng.dt <- round(swi.trend.dt[, .(total.change.pcnt, total.change.lbound.pcnt, total.change.ubound.pcnt)])
swi.pcnt.chng.vec <- paste(swi.pcnt.chng.dt$total.change.pcnt, ' (', swi.pcnt.chng.dt$total.change.lbound.pcnt, ', ', swi.pcnt.chng.dt$total.change.ubound.pcnt, ')', sep='')

# trend stats
swi.t = round(swi.trend.dt$tau,2)
swi.p = round(swi.trend.dt$pval,3)
swi.p[swi.p <= 0.001] <- '<0.001'
swi.tau.p.vec <- paste0(swi.t, ' (', swi.p, ')')


# NDVI - SWI CORRELATION =====================================
ndvi.swi.cor.dt <- ndvi.swi.cor.dt[swi.var != 'swi.anom.dt']
ndvi.swi.cor.dt$p <- round(ndvi.swi.cor.dt$p, 3)
ndvi.swi.cor.dt$p[ndvi.swi.cor.dt$p <= 0.001] <- '<0.001'
ndvi.swi.r.p.vec <- paste0(ndvi.swi.cor.dt$r, ' (', ndvi.swi.cor.dt$p, ')')

# COMBINE ===================================================

smry.dt <- cbind(swi.trend.dt[, .(period, lat.zone)], ndvi.abs.chng.vec, ndvi.pcnt.chng.vec, ndvi.t = ndvi.tau.p.vec,
                    swi.abs.chng.vec, swi.pcnt.chng.vec, swi.t = swi.tau.p.vec, cor.r.p = ndvi.swi.r.p.vec)

fwrite(smry.dt, 'output/summary_table_ndvi_swi_zonal_trends_cors.csv')
