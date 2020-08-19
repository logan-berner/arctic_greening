#-------------------------------------------------------------------------------------------------------------------
# PLOT LANDSAT NDVI TIME SERIES (BY ZONE), SWI TIME SERIES, AND MEAN ARCTIC NDVI VS SWI ANOMALIES, FOR 1985/2000 TO 2016
# AUTHOR: LOGAN BERNER, NAU
# DATE: 2020-04-09
#-------------------------------------------------------------------------------------------------------------------
rm(list=ls())
require(data.table)
require(dplyr)
require(lattice)
require(grid)
require(zyp)
require(reshape2)

setwd('/projects/arctic/users/lberner/arctic_greening/')

# LOAD DATA SETS ===============================================================================================

# zonal average ndvi time series
ndvi.zonal.yrly.mc.dt <- do.call("rbind", lapply(list.files('output/lsat_zonal_avg_timeseries/mc_reps/', full.names = T), fread))

# zonal average temperature data (MC reps)
swi.zonal.yrly.mc.dt <- do.call("rbind", lapply(list.files('output/swi_zonal_avg_timeseries/mc_reps/', full.names = T), fread))

# zonal correlations between NDVI and SWI anomalies
zonal.cor.dt <- fread('output/lsat_zonal_avg_swi_cors/lsat_zonal_avg_swi_cors_summary_fancy.csv')


# COMPUTE ZONAL AVERAGE NDVI TIME SERIES ACROSS MONTE CARLO ITERATION ===================================================================
ndvi.zonal.yrly.mc.dt <- ndvi.zonal.yrly.mc.dt[, lat.zone := factor(lat.zone, levels = c('High Arctic','Low Arctic','Oro Arctic','Arctic'))]
ndvi.zonal.yrly.dt <- ndvi.zonal.yrly.mc.dt[, .(ndvi.anom.med = median(ndvi.anom, na.rm = T), ndvi.anom.q025 = quantile(ndvi.anom, 0.025, na.rm = T), ndvi.anom.q975 = quantile(ndvi.anom, 0.975, na.rm=T),
                                                swi.anom.med = median(swi.anom, na.rm = T), swi.anom.q025 = quantile(swi.anom, 0.025, na.rm = T), swi.anom.q975 = quantile(swi.anom, 0.975, na.rm=T),
                                                swi.2yr.anom.med = median(swi.2yr.avg.anom, na.rm = T), swi.2yr.anom.q025 = quantile(swi.2yr.avg.anom, 0.025, na.rm = T), swi.2yr.anom.q975 = quantile(swi.2yr.avg.anom, 0.975, na.rm=T)), 
                                            by = c('period','lat.zone','year')]

ndvi.zonal.yrly.gte1985.dt <- ndvi.zonal.yrly.dt[period == '1985-2016']
ndvi.zonal.yrly.gte2000.dt <- ndvi.zonal.yrly.dt[period == '2000-2016']


# COMPUTE ZONAL AVERAGE SWI TIME SERIES ACROSS MONTE CARLO ITERATION ===================================================================
swi.zonal.yrly.mc.dt <- swi.zonal.yrly.mc.dt[, lat.zone := factor(lat.zone, levels = c('High Arctic','Low Arctic','Oro Arctic','Arctic'))]

swi.zonal.yrly.gte1985.mc.dt <-  swi.zonal.yrly.mc.dt[year >= 1985][, swi.anom := swi.avg - mean(swi.avg), by = 'lat.zone']
swi.zonal.yrly.gte1985.dt <- swi.zonal.yrly.gte1985.mc.dt[, .(swi.anom.med = median(swi.anom, na.rm = T), swi.anom.q025 = quantile(swi.anom, 0.025, na.rm = T), swi.anom.q975 = quantile(swi.anom, 0.975, na.rm = T)), by = c('lat.zone','year')]

swi.zonal.yrly.gte2000.mc.dt <-  swi.zonal.yrly.mc.dt[year >= 2000][, swi.anom := swi.avg - mean(swi.avg), by = 'lat.zone']
swi.zonal.yrly.gte2000.dt <- swi.zonal.yrly.gte2000.mc.dt[, .(swi.anom.med = median(swi.anom), swi.anom.q025 = quantile(swi.anom, 0.025), swi.anom.q975 = quantile(swi.anom, 0.975)), by = c('lat.zone','year')]

# PLOT DETAILS =====================================================================================================================

# correlation R values
r.lab.gte1985 <- bquote('r'[s]~' = ' ~ .(zonal.cor.dt[period == '1985-2016' & lat.zone == "Arctic"]$r))
r.lab.gte2000 <- bquote('r'[s]~' = ' ~ .(zonal.cor.dt[period == '2000-2016' & lat.zone == "Arctic"]$r))

# labels
lab.year <- 'Year'
lab.zone <- 'Arctic zone'
lab.swi <- expression('Mean SWI anom. ('*degree~'C)')
# lab.ndvi <- expression('Mean NDVI'[max]~' anom. (unitless)')
lab.ndvi <- expression('Mean NDVI anom. (unitless)')

# CEX 
my.cex=2.1
my.cex.adj=1

# colors
decade.cols <- c('skyblue4','skyblue2','darkorange1','firebrick2')
zone.cols <- rev(c('black','#d7191c','orange','#2b83ba'))

floor_decade <- function(value){ return(value - value %% 10) }

ndvi.zonal.yrly.gte1985.dt[, decade := floor_decade(year)]
ndvi.zonal.yrly.gte1985.dt[, decade.col := 'text']
ndvi.zonal.yrly.gte1985.dt[decade == 1980, decade.col := decade.cols[1]]
ndvi.zonal.yrly.gte1985.dt[decade == 1990, decade.col := decade.cols[2]]
ndvi.zonal.yrly.gte1985.dt[decade == 2000, decade.col := decade.cols[3]]
ndvi.zonal.yrly.gte1985.dt[decade == 2010, decade.col := decade.cols[4]]

ndvi.zonal.yrly.gte2000.dt[, decade := floor_decade(year)]
ndvi.zonal.yrly.gte2000.dt[, decade.col := 'text']
ndvi.zonal.yrly.gte2000.dt[decade == 2000, decade.col := decade.cols[3]]
ndvi.zonal.yrly.gte2000.dt[decade == 2010, decade.col := decade.cols[4]]

# keys 
zones <- c('Arctic','Oro Arctic','Low Arctic','High Arctic')
zone.key <- list(title='Arctic zone', corner=c(0.85,0.1),text=list(rev(zones)), 
                 lines=list(col=zone.cols,lwd=c(rep(1.5,3),3)), cex=1.6)

decade.key = list(title='Decade', corner=c(0.85,0.025), cex=1.6,
                  text=list(rev(c('1980s','1990s','2000s','2010s')), col = rev(decade.cols)))




# PLOT MEAN NDVI ANOMALIES FOR EACH ZONE ======================================================================================

panel.bands <- function(x,y,upper,lower,fill,col,subscripts, ...){
  upper <- upper[subscripts]
  lower <- lower[subscripts]
  panel.polygon(c(x,rev(x)), c(upper, rev(lower)), col = fill, border = F, ...)
}

plot.ndvi.ts.gte1985 <- xyplot(ndvi.anom.med ~ year, ndvi.zonal.yrly.gte1985.dt, groups = lat.zone, xlim = c(1984, 2017), ylim = c(-0.068,0.03),
                               type='l', lwd = c(1,1,1,2), key = zone.key,
                               ylab=list(lab.ndvi, cex=my.cex*my.cex.adj), xlab=list(lab.year, cex=my.cex*my.cex.adj), 
                               scales=list(cex=my.cex, tck = c(1,0), x=list(at=seq(1985,2015,10))),
                               upper = ndvi.zonal.yrly.gte1985.dt$ndvi.anom.q975, lower = ndvi.zonal.yrly.gte1985.dt$ndvi.anom.q025,
                               panel = function(x,y,...){ 
                                 panel.superpose(x,y, panel.groups = panel.bands, fill = adjustcolor(zone.cols, alpha.f = 0.3), ...)
                                 panel.abline(h = 0, lty=3, col = 'gray50')
                                 panel.xyplot(x,y, col = zone.cols, ...)}
                               )

plot.ndvi.ts.gte2000 <- xyplot(ndvi.anom.med ~ year, ndvi.zonal.yrly.gte2000.dt, groups = lat.zone, xlim = c(1999, 2017), ylim = c(-0.021,0.021),
                               type='l', lwd = c(1,1,1,2), 
                               ylab=list(lab.ndvi, cex=my.cex*my.cex.adj), xlab=list(lab.year, cex=my.cex*my.cex.adj), 
                               scales=list(cex=my.cex, tck = c(1,0), x=list(at=seq(2000,2015,5)), y=list(at=seq(-0.01,0.01,0.01))),
                               upper = ndvi.zonal.yrly.gte2000.dt$ndvi.anom.q975, lower = ndvi.zonal.yrly.gte2000.dt$ndvi.anom.q025,
                               panel = function(x,y,...){ 
                                 panel.superpose(x,y, panel.groups = panel.bands, fill = adjustcolor(zone.cols, alpha.f = 0.3), ...)
                                 panel.abline(h = 0, lty=3, col = 'gray50')
                                 panel.xyplot(x,y, col = zone.cols, ...)}
                               )

# PLOT MEAN ARCTIC SWI FROM 1985/2000 TO 2016 ===================================================================================

plot.swi.ts.gte1985 <- xyplot(swi.anom.med ~ year, swi.zonal.yrly.gte1985.dt, groups = lat.zone, xlim = c(1984, 2017), ylim = c(-6.4,6.4), type='l', lwd = c(1,1,1,2), 
                              ylab=list(lab.swi, cex=my.cex*my.cex.adj), xlab=list(lab.year, cex=my.cex*my.cex.adj), 
                              scales=list(cex=my.cex, tck = c(1,0), x=list(at=seq(1985,2015,10)), y=list(at=seq(-6,6,3))),
                              upper = swi.zonal.yrly.gte1985.dt$swi.anom.q975, lower = swi.zonal.yrly.gte1985.dt$swi.anom.q025,
                              panel = function(x,y,...){ 
                                panel.superpose(x,y, panel.groups = panel.bands, fill = adjustcolor(zone.cols, alpha.f = 0.3), ...)
                                panel.abline(h = 0, lty=3, col = 'gray50')
                                panel.xyplot(x,y, col = zone.cols, ...)}
                              )


plot.swi.ts.gte2000 <- xyplot(swi.anom.med ~ year, swi.zonal.yrly.gte2000.dt, groups = lat.zone, xlim = c(1999, 2017), ylim = c(-4.5,4.5), type='l', lwd = c(1,1,1,2), 
                               ylab=list(lab.swi, cex=my.cex*my.cex.adj), xlab=list(lab.year, cex=my.cex*my.cex.adj), 
                               scales=list(cex=my.cex, tck = c(1,0), x=list(at=seq(2000,2015,5))),
                               upper = swi.zonal.yrly.gte2000.dt$swi.anom.q975, lower = swi.zonal.yrly.gte2000.dt$swi.anom.q025,
                               panel = function(x,y,...){ 
                                 panel.superpose(x,y, panel.groups = panel.bands, fill = adjustcolor(zone.cols, alpha.f = 0.3), ...)
                                 panel.abline(h = 0, lty=3, col = 'gray50')
                                 panel.xyplot(x,y, col = zone.cols, ...)}
                              )

# PLOT MEAN ARCTIC NDVI VS SWI ANOMALIES 1985/2000 TO 2016 ================================================================

plot.ndvi.swi.gte1985 <- xyplot(ndvi.anom.med ~ swi.anom.med, ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic'], xlim = c(-6.4, 6.4), ylim = c(-0.068,0.03), type='p', pch = 21,
                                ylab=list(lab.ndvi, cex=my.cex*my.cex.adj), xlab=list(lab.swi, cex=my.cex*my.cex.adj),
                                fill = 'black', cex = 0.5, col = 'black',
                                scales=list(cex=my.cex, tck = c(1,0), x=list(at=seq(-6,6,3)), y=list(at=seq(-0.06,0.02,0.02))),
                                par.settings = list(layout.widths = list(ylab.axis.padding = 3)),
                                panel = function(x,y,...){
                                  panel.arrows(x0=ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$swi.anom.q025, 
                                               x1=ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$swi.anom.q975, 
                                               y0 = ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$ndvi.anom.med,
                                               y1 = ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$ndvi.anom.med,
                                               type = 'open', length=0, col = 'black')
                                  panel.arrows(x0=ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$swi.anom.med, 
                                               x1=ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$swi.anom.med, 
                                               y0 = ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$ndvi.anom.q025,
                                               y1 = ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$ndvi.anom.q975,
                                               type = 'open', length=0, col = 'black')
                                  panel.xyplot(x,y,...)
                                  panel.text(ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$swi.anom.med[c(8,32)], 
                                             ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$ndvi.anom.med[c(8,32)]+0.004, 
                                             ndvi.zonal.yrly.gte1985.dt[lat.zone == 'Arctic']$year[c(8,32)], cex = 1.25) 
                                  panel.text(0, -0.062, r.lab.gte1985, cex = my.cex*0.75)
                                })


plot.ndvi.swi.gte2000 <- xyplot(ndvi.anom.med ~ swi.anom.med, ndvi.zonal.yrly.gte2000.dt[lat.zone == 'Arctic'], xlim = c(-4.5, 4.5), ylim = c(-0.014,0.014), type='p', pch = 21,
                                ylab=list(lab.ndvi, cex=my.cex*my.cex.adj), xlab=list(lab.swi, cex=my.cex*my.cex.adj),
                                fill = 'black', cex = 0.5, col = 'black',
                                scales=list(cex=my.cex, tck = c(1,0), x=list(at=seq(-4,4,2)), y=list(at=seq(-0.01,0.01,0.01))),
                                par.settings = list(layout.widths = list(ylab.axis.padding = 3)),
                                panel = function(x,y,...){
                                  panel.arrows(x0=ndvi.zonal.yrly.gte2000.dt[lat.zone == 'Arctic']$swi.anom.q025, 
                                               x1=ndvi.zonal.yrly.gte2000.dt[lat.zone == 'Arctic']$swi.anom.q975, 
                                               y0 = ndvi.zonal.yrly.gte2000.dt[lat.zone == 'Arctic']$ndvi.anom.med,
                                               y1 = ndvi.zonal.yrly.gte2000.dt[lat.zone == 'Arctic']$ndvi.anom.med,
                                               type = 'open', length=0, col = 'black')
                                  panel.arrows(x0=ndvi.zonal.yrly.gte2000.dt[lat.zone == 'Arctic']$swi.anom.med, 
                                               x1=ndvi.zonal.yrly.gte2000.dt[lat.zone == 'Arctic']$swi.anom.med, 
                                               y0 = ndvi.zonal.yrly.gte2000.dt[lat.zone == 'Arctic']$ndvi.anom.q025,
                                               y1 = ndvi.zonal.yrly.gte2000.dt[lat.zone == 'Arctic']$ndvi.anom.q975,
                                               type = 'open', length=0, col = 'black')
                                  panel.xyplot(x,y,...)
                                  panel.text(0, -0.0125, r.lab.gte2000, cex = my.cex*0.75)
                                })


# ASSEMBLE 6 PANEL FIGURE -- NDVI zonal time series, SWI time series, NDVI vs SWI ============================================

pdf(file = 'figures/Landsat_NDVI_SWI_timeseries_cors.pdf', width = 16, height = 10)
# jpeg('figures/Landsat_NDVI_SWI_timeseries_cors.jpg', 16,10, units = 'in', res = 500)

# NDVI time series 
print(plot.ndvi.ts.gte1985, position=c(0.0,0.50,0.33,1), more=T)
print(plot.ndvi.ts.gte2000, position=c(0.0,0,0.33,0.50), more=T)

# SWI time series
print(plot.swi.ts.gte1985, position=c(0.33,0.50,0.66,1), more=T)
print(plot.swi.ts.gte2000, position=c(0.33,0.00,0.66,0.5), more=T)

# NDVI vs SWI
print(plot.ndvi.swi.gte1985, position=c(0.66,0.5,1.0,1.0), more=T)
print(plot.ndvi.swi.gte2000, position=c(0.66,0.0,1.0,0.50), more=T)

# labels
grid.text(expression(bold("(a)")), .12, 0.93, gp=gpar(fontsize=25))
grid.text(expression(bold("(b)")), .42, 0.93, gp=gpar(fontsize=25))
grid.text(expression(bold("(c)")), .785, 0.93, gp=gpar(fontsize=25))
grid.text(expression(bold("(d)")), .12, 0.43, gp=gpar(fontsize=25))
grid.text(expression(bold("(e)")), .42, 0.43, gp=gpar(fontsize=25))
grid.text(expression(bold("(f)")), .785, 0.43, gp=gpar(fontsize=25))

dev.off()

# END SCRIPT =====================================================================