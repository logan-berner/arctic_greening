rm(list=ls())
require(data.table)
require(dplyr)
require(lattice)
require(latticeExtra)
require(grid)

setwd('C:/Users/Logan/Google Drive/research/nau/nsf_arctic/arctic_greening/')

ndvi.trnd.site.cat.pcnt.by.zone <- fread('output/lsat_zonal_freq_trend_summary.csv', fill=T)
ndvi.swi.cor.site.cat.pcnt.by.zone <- fread('output/lsat_zonal_freq_swi_cors_summary.csv', fill=T)

# % TRENDING BY ZONE #--------------------------------------------------------------------------------
ndvi.trnd.site.cat.pcnt.by.zone[, trend.cat := factor(trend.cat, levels = c('browning (p<0.05)','browning (p<0.10)','no trend','greening (p<0.10)','greening (p<0.05)'))]
ndvi.trnd.site.cat.pcnt.by.zone[, lat.zone := factor(lat.zone, levels = c('High Arctic','Low Arctic','Oro Arctic','Arctic'),
                                                     labels = c('High','Low','Oro','Arctic'))]

ndvi.trnd.site.cat.pcnt.by.zone[, pcnt.sites.label := paste(round(pcnt.sites), '%', sep='')]
ndvi.trnd.site.cat.pcnt.by.zone[pcnt.sites.label == '0%' | pcnt.sites.label == '1%' | pcnt.sites.label == '2%', pcnt.sites.label:= '']

ndvi.trnd.site.cat.pcnt.by.zone[, pcnt.position := cumsum(pcnt.sites), by = c('period','lat.zone')]
ndvi.trnd.site.cat.pcnt.by.zone[, pcnt.position2 := lag(pcnt.position,1), by = c('period','lat.zone')]
ndvi.trnd.site.cat.pcnt.by.zone[, pcnt.position := (pcnt.position + pcnt.position2)/2]
ndvi.trnd.site.cat.pcnt.by.zone[is.na(pcnt.position), pcnt.position := 2]

ndvi.trnd.site.cat.pcnt.by.zone[, n.sites.zone := round(sum(n.sites)), by = c('period','lat.zone')]

ndvi.trnds.gte1985.site.by.zone <- ndvi.trnd.site.cat.pcnt.by.zone[period == '1985-2016']
ndvi.trnds.gte2000.site.by.zone <- ndvi.trnd.site.cat.pcnt.by.zone[period == '2000-2016']

# % COR CATEGORY BY ZONE #--------------------------------------------------------------------------------
ndvi.swi.cor.site.cat.pcnt.by.zone <- ndvi.swi.cor.site.cat.pcnt.by.zone[swi.var == 'swi']
ndvi.swi.cor.site.cat.pcnt.by.zone <- ndvi.swi.cor.site.cat.pcnt.by.zone[, cor.cat := factor(cor.cat, levels = c('negative (p<0.05)','negative (p<0.10)','no cor.','positive (p<0.10)','positive (p<0.05)'))]
ndvi.swi.cor.site.cat.pcnt.by.zone <- ndvi.swi.cor.site.cat.pcnt.by.zone[, lat.zone := factor(lat.zone, levels = c('High Arctic','Low Arctic','Oro Arctic','Arctic'),
                                                                                              labels = c('High','Low','Oro','Arctic'))]

ndvi.swi.cor.site.cat.pcnt.by.zone[, pcnt.sites.label := paste(round(pcnt.sites), '%', sep='')]
ndvi.swi.cor.site.cat.pcnt.by.zone[pcnt.sites.label == '0%' | pcnt.sites.label == '1%' | pcnt.sites.label == '2%', pcnt.sites.label := '']

# determine label positition 
ndvi.swi.cor.site.cat.pcnt.by.zone[, pcnt.position := cumsum(pcnt.sites), by = c('period','lat.zone')]
ndvi.swi.cor.site.cat.pcnt.by.zone[, pcnt.position2 := lag(pcnt.position,1), by = c('period','lat.zone')]
ndvi.swi.cor.site.cat.pcnt.by.zone[, pcnt.position := (pcnt.position + pcnt.position2)/2]
ndvi.swi.cor.site.cat.pcnt.by.zone[is.na(pcnt.position), pcnt.position := 2]

# detremine number of sites per zone 
ndvi.swi.cor.site.cat.pcnt.by.zone[, n.sites.zone := round(sum(n.sites)), by = c('period','lat.zone')]

# break into periods for easier plotting
ndvi.swi.cor.gte1985.site.by.zone <- ndvi.swi.cor.site.cat.pcnt.by.zone[period == '1985-2016']
ndvi.swi.cor.gte2000.site.by.zone <- ndvi.swi.cor.site.cat.pcnt.by.zone[period == '2000-2016']


# PLOT DETAILS #--------------------------------------------------------------------------------

xlab.zone <- 'Arctic zone'
# ylab.pcnt.zone <- expression('(% of sites)')
ylab.pcnt.zone <- '% of sites'

my.cex=1.5
my.cex.adj=1

trend.cols <- c('lightsalmon4','lightsalmon3','gray50','springgreen3','springgreen4')
trend.key <- list(title='NDVI trend', corner=c(0.9,0.9),text=list(rev(levels(ndvi.trnds.gte2000.site.by.zone$trend.cat))), rect=list(col=rev(trend.cols)), cex=my.cex, ncols =5)

## define a lattice "axis function"
axis.L <-function(side, ..., line.col){
    if (side %in% c("bottom", "left")) {
      col <- trellis.par.get("axis.text")$col
      axis.default(side, ..., line.col = col)
      if (side == "bottom")
        grid::grid.lines(y = 0)
      if (side == "left")
        grid::grid.lines(x = 0)
    }
  }

## hide panel and strip borders by using col = NA
sty <- list()
sty$axis.line$col <- NA
sty$strip.border$col <- NA
sty$strip.background$col <- NA


# CREATE LEGEND FIGURE --------------------------------------------------------------------------
legend.dt <- data.table(cat=c('negative(p < 0.05)','negative (p < 0.10)', 'none','positive (p < 0.10)', 'positive (p < 0.05)'),
                        cat.label=c("negative",'negative', 'none','positive', 'positive'), 
                        cat.sig.lab=c(expression(alpha~'= 0.05'), expression(alpha~'= 0.10'), expression(alpha~'> 0.10'), expression(alpha~'= 0.10'), expression(alpha~'= 0.05')),
                        blnk = 1, 
                        pcnt = rep(20,5), pcnt.position = seq(10,90,20))

plot.legend <- barchart(pcnt~blnk, legend.dt, groups = cat, stack = T, horizontal = F, col = trend.cols, scales=list(draw=F),
                        ylab = '', xlab='', tck = c(1,0), par.settings=list(axis.line=list(col="transparent")))

plot.legend <- plot.legend + layer(panel.text(x-0.2, legend.dt$pcnt.position, legend.dt$cat.label, srt = 90,  
                                              data = legend.dt, cex=0.85, col='white'))

plot.legend <- plot.legend + layer(panel.text(x+0.2, legend.dt$pcnt.position, legend.dt$cat.sig.lab, srt = 90,  
                                              data = legend.dt, cex=0.85, col='white'))

plot.legend
#-------------------------------------------------------------------------------------------------------------------
# PLOT TRENDS
plot.pcnt.trnds.gte1985 <- barchart(pcnt.sites ~ lat.zone, ndvi.trnds.gte1985.site.by.zone, groups = trend.cat, stack = T, col = trend.cols,
                                    ylab = '', xlab='', scales=list(x=list(draw=F), cex=my.cex, tck = c(1,0)), par.settings = sty, axis = axis.L,
                                    panel=function(x,y,groups,subscripts,...){
                                      panel.barchart(x,y,groups=groups,subscripts=subscripts,...)
                                      panel.abline(v = 3.5, lty=3, lwd=2)
                                    })

plot.pcnt.trnds.gte1985 <- plot.pcnt.trnds.gte1985 + layer(panel.text(x, ndvi.trnds.gte1985.site.by.zone$pcnt.position, ndvi.trnds.gte1985.site.by.zone$pcnt.sites.label, 
                                                                      data = ndvi.trnds.gte1985.site.by.zone, cex=0.85, col='white'))

plot.pcnt.trnds.gte1985 <- plot.pcnt.trnds.gte1985 + layer(panel.text(x, 104, ndvi.trnds.gte1985.site.by.zone$n.sites.zone, 
                                                                      data = ndvi.trnds.gte1985.site.by.zone, cex=0.85, col='black'))

plot.pcnt.trnds.gte2000 <- barchart(pcnt.sites ~ lat.zone, ndvi.trnds.gte2000.site.by.zone, groups = trend.cat, stack =T, col = trend.cols, 
                                    ylab='', xlab='', scales=list(cex=my.cex, tck = c(1,0)), par.settings = sty, axis = axis.L,
                                    panel=function(x,y,groups,subscripts,...){
                                      panel.barchart(x,y,groups=groups,subscripts=subscripts,...)
                                      panel.abline(v = 3.5, lty=3, lwd=2)
                                    })
plot.pcnt.trnds.gte2000 <- plot.pcnt.trnds.gte2000 + layer(panel.text(x, ndvi.trnds.gte2000.site.by.zone$pcnt.position, ndvi.trnds.gte2000.site.by.zone$pcnt.sites.label, 
                                                                      data = ndvi.trnds.gte2000.site.by.zone, cex=0.85, col='white'))

plot.pcnt.trnds.gte2000 <- plot.pcnt.trnds.gte2000 + layer(panel.text(x, 104, ndvi.trnds.gte2000.site.by.zone$n.sites.zone, 
                                                                      data = ndvi.trnds.gte2000.site.by.zone, cex=0.85, col='black'))



# PERCENT OF SITES WITH LANDSAT NDVI - SWI CORELATION BY ZONE #--------------------------------------------------------------------------------

plot.pcnt.clim.cor.gte1985 <- barchart(pcnt.sites ~ lat.zone, ndvi.swi.cor.gte1985.site.by.zone, groups = cor.cat, stack =T, col = trend.cols, 
                                       ylab='', xlab='', scales=list(x = list(draw=F), y = list(labels=NULL), tck = c(1,0)), par.settings = sty, axis = axis.L,
                                       panel=function(x,y,groups,subscripts,...){
                                         panel.barchart(x,y,groups=groups,subscripts=subscripts,...)
                                         panel.abline(v = 3.5, lty=3, lwd=2)
                                       })

plot.pcnt.clim.cor.gte1985 <- plot.pcnt.clim.cor.gte1985 + layer(panel.text(x, ndvi.swi.cor.gte1985.site.by.zone$pcnt.position, ndvi.swi.cor.gte1985.site.by.zone$pcnt.sites.label, 
                                                                            data = ndvi.swi.cor.gte1985.site.by.zone, cex=0.85, col='white'))

plot.pcnt.clim.cor.gte1985 <- plot.pcnt.clim.cor.gte1985 + layer(panel.text(x, 104, ndvi.swi.cor.gte1985.site.by.zone$n.sites.zone, 
                                                                            data = ndvi.swi.cor.gte1985.site.by.zone, cex=0.85, col='black'))

plot.pcnt.clim.cor.gte2000 <- barchart(pcnt.sites ~ lat.zone, ndvi.swi.cor.gte2000.site.by.zone, groups = cor.cat, stack =T, col = trend.cols, 
                                       ylab='', xlab='', scales=list(y=list(labels = NULL), cex=my.cex, tck = c(1,0)), par.settings = sty, axis = axis.L,
                                       panel=function(x,y,groups,subscripts,...){
                                         panel.barchart(x,y,groups=groups,subscripts=subscripts,...)
                                         panel.abline(v = 3.5, lty=3, lwd=2)
                                       })
plot.pcnt.clim.cor.gte2000 <- plot.pcnt.clim.cor.gte2000 + layer(panel.text(x, ndvi.swi.cor.gte2000.site.by.zone$pcnt.position, ndvi.swi.cor.gte2000.site.by.zone$pcnt.sites.label, 
                                                                            data = ndvi.swi.cor.gte2000.site.by.zone, cex=0.85, col='white'))

plot.pcnt.clim.cor.gte2000 <- plot.pcnt.clim.cor.gte2000 + layer(panel.text(x, 104, ndvi.swi.cor.gte2000.site.by.zone$n.sites.zone, 
                                                                            data = ndvi.swi.cor.gte2000.site.by.zone, cex=0.85, col='black'))



pdf('figures/revised_figures/Landsat_NDVI_trend_NDVIvsSWI_cor_zonal_freq.pdf', 9,6)
# jpeg('figures/revised_figures/Landsat_NDVI_trend_NDVIvsSWI_cor_zonal_freq.jpg', 10,6, units = 'in', res = 400)
# ndvi trends
print(plot.pcnt.trnds.gte1985, position=c(0.05,0.50,0.51,0.95), more=T)
print(plot.pcnt.trnds.gte2000, position=c(0.05,0.03,0.51,0.56), more=T)

# ndvi-clim cors
print(plot.pcnt.clim.cor.gte1985, position=c(0.45,0.50,0.85,0.95), more=T)
print(plot.pcnt.clim.cor.gte2000, position=c(0.45,0.03,0.85,0.56), more=T)

# add legend
print(plot.legend, position = c(0.89, 0.0, 1, 1))

# labels
grid.text(expression(bold("(a)")), .14, 0.93, gp=gpar(fontsize=18))
grid.text(expression(bold("(b)")), .49, 0.93, gp=gpar(fontsize=18))
grid.text(expression(bold("(c)")), .14, .51, gp=gpar(fontsize=18))
grid.text(expression(bold("(d)")), .49, .51, gp=gpar(fontsize=18))

grid.text(expression('NDVI'['max']~'trend'), 0.30, 0.97, gp=gpar(fontsize=19), rot=0)
grid.text(expression('NDVI'['max']~'- SWI correlation'), 0.64, 0.97, gp=gpar(fontsize=19), rot=0)
grid.text('1985 - 2016', 0.83, 0.73, gp=gpar(fontsize=19), rot=90)
grid.text('2000 - 2016', 0.83, 0.30, gp=gpar(fontsize=19), rot=90)
grid.text('Percent of sites', 0.05, 0.55, gp=gpar(fontsize=19), rot=90)
grid.text('Arctic zone', 0.48, 0.02, gp=gpar(fontsize=20))

dev.off()


# END SCRIPT #--------------------------------------------------------------------------------