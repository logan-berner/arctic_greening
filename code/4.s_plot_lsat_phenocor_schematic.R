# This R script create a 2 panel figure that shows the (a) Landsat seasonal NDVI phenology and (b) single-year example of the phenological correction 
# Date: 2020-07-02
rm(list=ls())
require(tidyr)
require(dplyr)
require(data.table)
require(lattice)
require(ggplot2)
require(ggpubr)
setwd('/projects/arctic/users/lberner/arctic_greening/')
source('/home/lb968/code/arctic_greening/')

# READ IN NDVI TIME SERIES FROM GEE AND CLEAN UP / REORDER COLUMNS ==============================================
lsat.dt <- fread('output/pheno_timeseries/tundra_landsat_ndvi_pheno_corrected_timeseries_rep_0001.csv')
lsat.max.dt <- fread('output/ndvi_max_timeseries/tundra_site_lsat_ndvi_max_timeseries_rep_0001.csv')

pheno.curves.dt <- fread('output/pheno_curves/pheno_curves_rep_0001.csv')

# rename some columns
lsat.dt <- setnames(lsat.dt, 'ndvi.xcal','ndvi')
lsat.max.dt <- setnames(lsat.max.dt, 'pos.doy.avg', 'doy')
lsat.max.dt <- setnames(lsat.max.dt, 'ndvi.max', 'ndvi')

# find example site 
lsat.ss.dt <- lsat.dt[2000:3000,]
ggplot(lsat.ss.dt, aes(doy, ndvi)) + geom_point() + facet_grid(~site)
unique(lsat.ss.dt$site)
# '0000cb391533f39123f7'
# '0000d93deb6ee1acb7cd'
exmpl.site <- '0000001115ae6b8f7f82'
# exmpl.site <- '00000013b0014a873f63'
# exmpl.site <-  '000000248a602081a1cc'
site.obs.dt <- lsat.dt[site == exmpl.site]
site.max.dt <- lsat.max.dt[site == exmpl.site]
site.curves.dt <- pheno.curves.dt[site == exmpl.site]

exmpl.yr <- 2000
site.obs.exmpl.dt <- site.obs.dt[year == exmpl.yr]
site.max.exmpl.dt <- site.max.dt[year == exmpl.yr]
site.curve.exmpl.dt <- site.curves.dt[focal.yr == exmpl.yr]

# --------------------------------------------------------------------------------
fig1 <- ggplot(site.obs.dt, aes(doy, ndvi)) + ylim(0.25,0.75) + labs(y='Landsat NDVI', x='Day of Year') + ggtitle('All years')
fig1 <- fig1 + geom_point(aes(fill = year), pch=21, color = 'black', size = 2) # plot obs
fig1 <- fig1 + scale_fill_gradient2(name = 'Observation', low='violetred4',mid='brown', high = 'green4', midpoint = 2000) # add color to points
fig1 <- fig1 + geom_line(data = site.curves.dt, mapping = aes(doy, spl.fit, group = focal.yr, color = focal.yr)) # add pheno curves
fig1 <- fig1 + scale_color_gradient2(name = 'Curve', low='violetred4',mid = 'brown', high = 'green4', midpoint = 2000) # add color to curves
fig1 <- fig1 + theme_bw() + theme(legend.position = c(0.6, 0.25), legend.box = 'horizontal', legend.text=element_text(size=8), legend.title=element_text(size=10), 
                                   axis.text=element_text(size=12), axis.title=element_text(size=14),
                                   plot.title=element_text(hjust = 0.5))

fig2 <- ggplot(data = site.obs.exmpl.dt, mapping = aes(doy, ndvi)) + labs(y='Landsat NDVI', x='Day of Year') + ylim(0.25,0.75) + ggtitle('Example: 2000') # set up plot
fig2 <- fig2 + geom_line(data = site.curve.exmpl.dt, mapping = aes(doy, spl.fit), color = 'brown') # add pheno curves
fig2 <- fig2 + geom_line(data = site.curve.exmpl.dt, mapping = aes(doy, spl.fit.max), lty=2) # add dashed horizontal line for NDVImax
fig2 <- fig2 + annotate(geom = 'text', x = 165, y = max(site.obs.exmpl.dt$spl.fit.max)+0.02, label=expression('typical NDVI'[max])) # add NDVI max label
fig2 <- fig2 + geom_segment(data = site.obs.exmpl.dt, mapping = aes(x=doy, xend=doy, y=spl.fit, yend=spl.fit.max), color='lightblue') # add pheno adjust lines
fig2 <- fig2 + geom_point(pch=21, color='black', fill='brown', size = 2) # add obs points
fig2 <- fig2 + geom_point(data = site.max.exmpl.dt, mapping = aes(doy, ndvi), pch=21, color='black', fill='black', size = 4) # add max point
fig2 <- fig2 + geom_errorbar(data = site.max.exmpl.dt, mapping = aes(ymin=ndvi.max.lower, ymax=ndvi.max.upper)) # add max error bars
fig2 <- fig2 + theme_bw() + theme(legend.position = c(0.8, 0.2), axis.text=element_text(size=12), axis.title=element_text(size=14), plot.title=element_text(hjust = 0.5))

combo.fig <- ggarrange(fig1, fig2, labels = c("(a)", "(b)"), label.x = 0.15, label.y = 0.90, ncol = 2, nrow = 1)

jpeg('figures/Landsat_pheno_correction_schematic_multicurve.jpg', 9, 4, units = 'in', res = 400)
combo.fig
dev.off()

# END SCRIPT #--------------------------------------------------------------------------------------------------------
# + theme(legend.position = c(0.8,0.2), axis.text=element_text(size=12), axis.title=element_text(size=14),
#         legend.text=element_text(size=10), legend.title=element_text(size=12))