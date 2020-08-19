# ABOUT THIS SCRIPT  ==============================================================================================================
# This R script summarizes Landsat cross-sensor calibration random forest models
# AUTHOR: LOGAN BERNER, NAU
# DATE: 2020-02-25

# SET UP WORKSPACE ==============================================================================================================
rm(list=ls())
require(data.table)
require(reshape2)
require(ggplot2)
require(ggpubr)

setwd('/projects/arctic/users/lberner/arctic_greening/')

# READ IN FILES
biome.trnd.dt <- do.call("rbind", lapply(list.files('output/lsat_sample_size_test/trend_mag_mcreps/', full.names = T), fread))
site.trnd.dt <- do.call("rbind", lapply(list.files('output/lsat_sample_size_test/trend_freq_mcreps/', full.names = T), fread))

select.pcnts <- c(100, 1000,10000,40000)

# TRENDS IN ZONAL AVERAGE NDVI -------------------------------------------------------------------
biome.trnd.dt[, int.2000 := int + slope * 2000]
biome.trnd.dt[, total.change.pcnt := (slope * 16)/int.2000*100]

biome.trnd.smry.dt <- biome.trnd.dt[, .(slope=median(slope), slope.q025=quantile(slope,0.025), slope.q975=quantile(slope,0.975),
                                        total.change.pcnt=median(total.change.pcnt), total.change.pcnt.q025=quantile(total.change.pcnt,0.025), 
                                        total.change.pcnt.q975=quantile(total.change.pcnt,0.975)), by = sample.size]

biome.trnd.smry.dt[, CI95 := total.change.pcnt.q975 - total.change.pcnt.q025]

# FREQ OF NDVI TRENDS BY ZONE ------------------------------------------------------------------------
site.trnd.smry.dt <- site.trnd.dt[, .(pcnt.sites=round(median(pcnt.sites),1), pcnt.sites.q025=round(quantile(pcnt.sites,0.025),1), pcnt.sites.q975=round(quantile(pcnt.sites,0.975),1)),
                                  by = c('sample.size','trend.cat')]

site.trnd.smry.dt[, CI95 := pcnt.sites.q975-pcnt.sites.q025]


# SUMMARY TABLE =========================================================================

biome.trnd.smry.fancy.dt <- biome.trnd.smry.dt[sample.size %in% select.pcnts]
biome.trnd.smry.fancy.dt<- biome.trnd.smry.fancy.dt[, c('slope','slope.q025','slope.q975') := NULL]
biome.trnd.smry.fancy.dt[, total.change.pcnt := paste0(sprintf('%.2f', total.change.pcnt), ' [', sprintf('%.2f', total.change.pcnt.q025), ', ', sprintf('%.2f', total.change.pcnt.q975), ']')]
biome.trnd.smry.fancy.dt <- biome.trnd.smry.fancy.dt[, c('total.change.pcnt.q025','total.change.pcnt.q975') := NULL]

site.trnd.smry.fancy.dt <- site.trnd.smry.dt[sample.size %in% select.pcnts]
site.trnd.smry.fancy.dt[, pcnt.sites := paste0(sprintf('%.2f', pcnt.sites), ' [', sprintf('%.2f', pcnt.sites.q025), ', ', sprintf('%.2f', pcnt.sites.q975), ']')]
site.trnd.smry.fancy.dt <- site.trnd.smry.fancy.dt[, c('pcnt.sites.q025','pcnt.sites.q975') := NULL]
site.trnd.smry.fancy.dt <- dcast(site.trnd.smry.fancy.dt, sample.size ~ trend.cat)
site.trnd.smry.fancy.dt <- data.table(site.trnd.smry.fancy.dt)

smry.table <- biome.trnd.smry.fancy.dt[site.trnd.smry.fancy.dt, on= 'sample.size']

fwrite(smry.table, 'output/lsat_sample_size_test/sample_size_analysis_smry.csv')


# FIGURES =========================================================================

# mean trend figure 
biome.trnd.fig <- ggplot(biome.trnd.smry.dt, aes(sample.size, total.change.pcnt))
biome.trnd.fig <- biome.trnd.fig + geom_ribbon(aes(x = sample.size, ymin=total.change.pcnt.q025, ymax=total.change.pcnt.q975), alpha = 0.5)
biome.trnd.fig <- biome.trnd.fig + geom_line() + labs(y="Change in mean Arctic NDVI (%)", x="Sample size") 
biome.trnd.fig <- biome.trnd.fig + theme_bw() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14))

# trend freq figure
site.trnd.smry.dt <- site.trnd.smry.dt[trend.cat != 'insig']
site.trnd.smry.dt[, trend.cat := factor(trend.cat, levels = c('greening','browning'))]
trend.cat.cols <- c('green','brown')
site.trnd.fig <- ggplot(site.trnd.smry.dt, aes(sample.size, pcnt.sites, group = trend.cat))
site.trnd.fig <- site.trnd.fig + geom_ribbon(aes(x=sample.size, ymin=pcnt.sites.q025, ymax=pcnt.sites.q975, fill=trend.cat), alpha = 0.5)
site.trnd.fig <- site.trnd.fig + geom_line(aes(colour=trend.cat)) + scale_fill_manual(guide = F, values = trend.cat.cols) + scale_color_manual(name = 'NDVI trend', values = trend.cat.cols)
site.trnd.fig <- site.trnd.fig + labs(y="Percent of Arctic sampling sites", x="Sample size")
site.trnd.fig <- site.trnd.fig + theme_bw() + theme(legend.position = c(0.8,0.45), axis.text=element_text(size=12), axis.title=element_text(size=14), 
                     legend.text=element_text(size=12), legend.title=element_text(size=14))

# plot width of 95% CI as a function of sample size 

biome.CI.fig <- ggplot(biome.trnd.smry.dt, aes(sample.size, CI95))
biome.CI.fig <- biome.CI.fig + geom_line() + labs(y="Width of 95% CI (%)", x="Sample size") 
biome.CI.fig <- biome.CI.fig + theme_bw() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14))

site.CI.fig <- ggplot(site.trnd.smry.dt, aes(sample.size, CI95, group = trend.cat))
site.CI.fig <- site.CI.fig + geom_line(aes(colour=trend.cat)) + scale_fill_manual(guide = F, values = trend.cat.cols) + scale_color_manual(name = 'NDVI trend', values = trend.cat.cols)
site.CI.fig <- site.CI.fig + labs(y="Width of 95% CI (%)", x="Sample size")
site.CI.fig <- site.CI.fig + theme_bw() + theme(legend.position='none', axis.text=element_text(size=12), axis.title=element_text(size=14))

# combine figures
combo.fig <- ggarrange(biome.trnd.fig, site.trnd.fig, labels = c("(a)", "(b)"), label.x = 0.20, label.y = 0.98, ncol = 2, nrow = 1)

jpeg('figures/Landsat_NDVI_trend_sample_size.jpg', width = 8, height = 3.75, units = 'in', res = 400)
combo.fig
dev.off()


# combine figures
combo.fig <- ggarrange(biome.trnd.fig, site.trnd.fig, biome.CI.fig, site.CI.fig, labels = c("(a)", "(b)","(c)","(d)"), label.x = 0.15, label.y = 0.98, ncol = 2, nrow = 2)

jpeg('figures/Landsat_NDVI_trend_sample_size_4panel.jpg', width = 9, height = 7, units = 'in', res = 400)
combo.fig
dev.off()
