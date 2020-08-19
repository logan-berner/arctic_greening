# ABOUT THIS SCRIPT =================================================================================================
# This R script summarizes the frequency and co-occurrence greening, warming, and ndvi-swi correlations for sites in the Arctic. 
# Author: Logan Berner, NAU
# Date: 2019-10-10

# SET UP WORKSPACE =================================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(lattice)
require(grid)
require(data.table)
require(dplyr)
setwd('/projects/arctic/users/lberner/arctic_greening/')


# LOAD DATA SETS AND COMBINE INTO ONE DATA TABLE =================================================================================================
# read in data sets and combine into one data table
swi.cor.files <- list.files('output/lsat_site_swi_cors/mc_reps/', full.names = T)
ndvi.trend.files <- list.files('output/lsat_site_trends/mc_reps/', full.names = T)  
swi.trend.files <- list.files('output/swi_site_trends/mc_reps/', full.names = T)

mc.reps <- 1000
i=1

for (i in mc.reps){
  ndvi.trend.dt <- fread(ndvi.trend.files[i])
  ndvi.swi.cor.dt <- fread(swi.cor.files[i])
  swi.trend.dt <- fread(swi.trend.files[i])
  
  # select variables of interest
  ndvi.trend.dt <- ndvi.trend.dt[, c('period','site','trend.cat')]
  ndvi.swi.cor.dt <- ndvi.swi.cor.dt[swi.var == 'swi', c('period','lat.zone','site','cor.cat')]
  swi.trend.dt <- swi.trend.dt[, c('period','site','swi.trend.cat')]
  
  # combine into one data table
  site.dt <- swi.trend.dt[ndvi.trend.dt, on = c('period','site')]
  site.dt <- ndvi.swi.cor.dt[site.dt, on = c('period','site')]
  
  # type cast
  site.dt[, cor.cat := as.character(cor.cat)]
  site.dt[, swi.trend.cat := as.character(swi.trend.cat)]
  site.dt[, trend.cat := as.character(trend.cat)]
  
  # reclassify trends and correlations
  site.dt[cor.cat == 'positive.sig.p5' | cor.cat == 'positive.sig.p10', cor.cat := 'positive']
  site.dt[cor.cat == 'negative.sig.p5' | cor.cat == 'negative.sig.p10', cor.cat := 'negative']
  site.dt[cor.cat == 'insig', cor.cat := 'none']
  site.dt[, cor.cat := factor(cor.cat, levels = c('negative','none','positive'))]
          
  site.dt[swi.trend.cat == 'warming.sig.p5' | swi.trend.cat == 'warming.sig.p10', swi.trend.cat := 'warming']
  site.dt[swi.trend.cat == 'cooling.sig.p5' | swi.trend.cat == 'cooling.sig.p10', swi.trend.cat := 'cooling']
  site.dt[swi.trend.cat == 'warming.insig' | swi.trend.cat == 'cooling.insig', swi.trend.cat := 'none']
  site.dt[, swi.trend.cat := factor(swi.trend.cat, levels = c('cooling','none','warming'))]
  
  site.dt[trend.cat == 'greening.sig.p5' | trend.cat == 'greening.sig.p10', trend.cat := 'greening']
  site.dt[trend.cat == 'browning.sig.p5' | trend.cat == 'browning.sig.p10', trend.cat := 'browning']
  site.dt[trend.cat == 'insig', trend.cat := 'none']
  site.dt[, trend.cat := factor(trend.cat, levels = c('browning','none','greening'))]
  
  
  # FREQUENCY OF NDVI - SWI CORRELATIONS BY NDVI TREND ===========================================================================
  cor.by.ndvi.trend.dt <- site.dt[is.na(cor.cat)==F, .(n.sites=.N), by = c('period','trend.cat','cor.cat')]
  cor.by.ndvi.trend.dt <- cor.by.ndvi.trend.dt[order(period,trend.cat,cor.cat)]
  cor.by.ndvi.trend.dt[, n.sites.trend := sum(n.sites), by = c('period','trend.cat')][, pcnt.sites := n.sites / n.sites.trend * 100]
  cor.by.ndvi.trend.dt
}



# # COMPUTE FREQUERCY OF NDVI TRENDS BY WARMING TRENDS =============================================================================
# # by zone
# ndvi.by.swi.cats.zonal.dt <- site.dt %>% group_by(period, lat.zone, swi.trend.cat, trend.cat) %>% summarize(n.ndvi.cat = n()) %>% 
#   group_by(period, lat.zone, swi.trend.cat) %>% mutate(n.tot = sum(n.ndvi.cat), pcnt.sites = n.ndvi.cat/n.tot*100)
# 
# 
# # by biome
# ndvi.by.swi.cats.biome.dt <- site.dt %>% group_by(period, swi.trend.cat, trend.cat) %>% summarize(n.ndvi.cat = n()) %>% 
#   group_by(period, swi.trend.cat) %>% mutate(n.tot = sum(n.ndvi.cat), pcnt.sites = n.ndvi.cat/n.tot*100) %>%
#   mutate(lat.zone = 'Arctic')
# 
# ndvi.by.swi.cat.dt <- rbind(ndvi.by.swi.cats.biome.dt, ndvi.by.swi.cats.zonal.dt)
# 
# ndvi.by.warming.cat.dt <- ndvi.by.swi.cat.dt %>% filter(swi.trend.cat == 'warming')
# # ndvi.by.warming.cat.dt <- ndvi.by.warming.cat.dt %>% group_by(period, lat.zone) %>% mutate(pcnt.position = cumsum(pcnt.sites), pcnt.sites.lab = paste(round(pcnt.sites),'%',sep=''))
# # ndvi.by.warming.cat.dt
# 
# ndvi.by.warming.cat.gte1985.dt <- ndvi.by.warming.cat.dt %>% filter(period == "1985-2016")
# ndvi.by.warming.cat.gte2000.dt <- ndvi.by.warming.cat.dt %>% filter(period == "2000-2016")
# 
# ndvi.by.warming.cat.gte1985.dt
# ndvi.by.warming.cat.gte2000.dt
# 
# # figure
# xlab.year <- 'Year'
# xlab.zone <- 'Arctic zone'
# ylab.pcnt.zone <- expression('Trend in NDVI'[max]~' (% warming of sites)')
# my.cex=1.75
# my.cex.adj=1
# trend.cols <- c('lightsalmon4','gray50','springgreen4')
# 
# 
# plot.pcnt.trend.cat.by.warming.gte1985 <- barchart(pcnt.sites ~ lat.zone | period, ndvi.by.warming.cat.gte1985.dt, groups = trend.cat, stack =T, col = trend.cols, 
#                                                    ylab=list(ylab.pcnt.zone, cex=my.cex*my.cex.adj), xlab=list(xlab.zone, cex=my.cex*my.cex.adj), scales=list(cex=my.cex, tck = c(1,0)),
#                                                    panel=function(x,y,groups,subscripts,...){
#                                                      panel.barchart(x,y,groups=groups,subscripts=subscripts,...)
#                                                      panel.abline(v = 3.5, lty=3, lwd=2)
#                                                    })
# plot.pcnt.trend.cat.by.warming.gte1985 <- plot.pcnt.trend.cat.by.warming.gte1985 + layer(panel.text(x, ndvi.by.warming.cat.gte1985.dt$pcnt.position-2, ndvi.by.warming.cat.gte1985.dt$pcnt.sites.lab, 
#                                                                                                     data = ndvi.by.warming.cat.gte1985.dt, cex=0.85, col='white'))
# plot.pcnt.trend.cat.by.warming.gte1985 <- plot.pcnt.trend.cat.by.warming.gte1985 + layer(panel.text(x, 104, ndvi.by.warming.cat.gte1985.dt$n.tot, 
#                                                                                                     data = ndvi.by.warming.cat.gte1985.dt, cex=0.85, col='black'))
# 
# plot.pcnt.trend.cat.by.warming.gte2000 <- barchart(pcnt.sites ~ lat.zone | period, ndvi.by.warming.cat.gte2000.dt, groups = trend.cat, stack =T, col = trend.cols, 
#                                                    ylab=list(ylab.pcnt.zone, cex=my.cex*my.cex.adj), xlab=list(xlab.zone, cex=my.cex*my.cex.adj), scales=list(cex=my.cex, tck = c(1,0)),
#                                                    panel=function(x,y,groups,subscripts,...){
#                                                      panel.barchart(x,y,groups=groups,subscripts=subscripts,...)
#                                                      panel.abline(v = 3.5, lty=3, lwd=2)
#                                                    })
# plot.pcnt.trend.cat.by.warming.gte2000 <- plot.pcnt.trend.cat.by.warming.gte2000 + layer(panel.text(x, ndvi.by.warming.cat.gte2000.dt$pcnt.position-2, ndvi.by.warming.cat.gte2000.dt$pcnt.sites.lab, 
#                                                                                                     data = ndvi.by.warming.cat.gte2000.dt, cex=0.85, col='white'))
# plot.pcnt.trend.cat.by.warming.gte2000 <- plot.pcnt.trend.cat.by.warming.gte2000 + layer(panel.text(x, 104, ndvi.by.warming.cat.gte2000.dt$n.tot, 
#                                                                                                     data = ndvi.by.warming.cat.gte2000.dt, cex=0.85, col='black'))
# 
# plot.pcnt.trend.cat.by.warming <- c(plot.pcnt.trend.cat.by.warming.gte1985, plot.pcnt.trend.cat.by.warming.gte2000, layout=c(2,1))
# 
# plot.pcnt.trend.cat.by.warming
# 
# jpeg('figures/lsat_ndvi_trnd_cat_by_warming.jpg', 10, 6, units = 'in', res=400)
# plot.pcnt.trend.cat.by.warming
# dev.off()
# 
# 
# # COMPUTE FREQUENCY OF WARMING TRENDS BY NDVI TREND  =============================================================================
# # by zone
# swi.by.ndvi.cats.zonal.dt <- site.dt %>% group_by(period, lat.zone, trend.cat, swi.trend.cat) %>% summarize(n.swi.cat = n()) %>% 
#   group_by(period, lat.zone, trend.cat) %>% mutate(n.tot = sum(n.swi.cat), pcnt.sites = n.swi.cat/n.tot*100)
# 
# # by biome
# swi.by.ndvi.cats.biome.dt <- site.dt %>% group_by(period, trend.cat, swi.trend.cat) %>% summarize(n.swi.cat = n()) %>% 
#   group_by(period, trend.cat) %>% mutate(n.tot = sum(n.swi.cat), pcnt.sites = n.swi.cat/n.tot*100) %>%
#   mutate(lat.zone = 'Arctic')
# 
# 
# swi.by.ndvi.cat.dt <- rbind(swi.by.ndvi.cats.biome.dt, swi.by.ndvi.cats.zonal.dt)
# swi.by.ndvi.cat.dt
# swi.by.ndvi.cats.biome.dt
# 
# 
# # COMPUTE FREQUENCY OF NDVI - SWI CORRELATION TYPES  =============================================================================
# cor.freq.dt <- site.dt %>% group_by(period, cor.cat) %>% summarise(n.sites = n()) %>% group_by(period) %>%
#   mutate(n.sites.tot = sum(n.sites)) %>% group_by(period) %>% mutate(pcnt.sites = n.sites / n.sites.tot * 100)
# 
# cor.freq.dt
# 
# cor.freq.zone.dt <- site.dt %>% group_by(lat.zone, period, cor.cat) %>% summarise(n.sites = n()) %>% group_by(lat.zone, period) %>%
#   mutate(n.sites.tot = sum(n.sites)) %>% group_by(lat.zone, period) %>% mutate(pcnt.sites = n.sites / n.sites.tot * 100)
# 
# cor.freq.zone.dt
# 
# # COMPUTE FREQUENC OF WARMING AND GREENING CO-OCCURRING  ============================================================================= 
# site.smry.dt <- site.dt %>% group_by(period, trend.cat, swi.trend.cat) %>% summarise(n.sites = n()) %>% 
#   group_by(period, trend.cat) %>% mutate(n.sites.ndvi.trend.cat = sum(n.sites)) %>% 
#   group_by(period, swi.trend.cat) %>% mutate(n.sites.swi.trend.cat = sum(n.sites)) %>% 
#   group_by() %>% mutate(pcnt.sites.ndvi.trend.cat = n.sites / n.sites.ndvi.trend.cat * 100, pcnt.sites.swi.trend.cat = n.sites / n.sites.swi.trend.cat * 100)
# 
# site.smry.dt
# 

