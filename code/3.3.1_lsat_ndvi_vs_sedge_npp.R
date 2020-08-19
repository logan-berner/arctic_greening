# ABOUT THIS SCRIPT =========================================================================================================
# This R script compares annual Landsat NDVI with graminoid productivity for sites on Bylot Island in northern Canada. 
# Author: Logan Berner, NAU
# Date: 2019-08-01

# SET UP WORKSPACE =========================================================================================================
rm(list=ls())
require(dplyr)
require(data.table)
require(tidyr)
require(lattice)
require(plotrix)
require(zyp)

setwd('/projects/arctic/users/lberner/arctic_greening/')
source('/home/lb968/code/arctic_greening/0.3_fun_shade_CIs.r')

# LOAD DATA SETS =========================================================================================================
npp <- fread('data/field_data/bylot_sedge_npp/bylot_bcv_biomass.csv')
npp <- npp[treatment == 'grazed']
npp <- na.omit(npp)
head(npp)

lsat.max <- fread('output/site_comparisons/lsat_ndvi_max_for_sedge_npp_subsites_mcreps.csv')
lsat.max

# check autocorrelation in the ANPP time series
npp.site.smry <- npp[, .(npp.gm2yr.med = median(npp.gm2yr, na.rm=T)), by='year']
npp.acf <- acf(npp.site.smry$npp.gm2yr.med)
plot(npp.acf)

# compute average date of harvests
npp %>% group_by() %>% summarise(date.med = mean(date), date.sd = sd(date), n.quadrats = n())

# START MONTE CARLO ITERATIONS =========================================================================================================
bio.ts.mc.list <- list()
bio.cor.mc.list <- list()
bio.trnds.mc.list <- list()
mc.reps <- 1000
for (i in 1:mc.reps){ 

  # subsample quadrats in sample
  npp.mc <- npp[, .SD[sample(.N,.N*0.9)], by = 'year']

  # summarize ANPP across plots
  npp.mc.smry <- npp.mc[, .(npp.gm2yr.med = median(npp.gm2yr, na.rm=T), 
                            npp.gm2yr.log.med = median(log(npp.gm2yr), na.rm=T), n.plots = .N),
                        by='year']
  
  # detrend ANPP time series
  npp.mc.smry[, npp.gm2yr.med.dt := lm(npp.gm2yr.med ~ year)$resid]
  npp.mc.smry[, npp.gm2yr.log.med.dt := lm(npp.gm2yr.log.med ~ year)$resid]
  
  # lsat NDVI permutation of xcal and pheno-correction 
  lsat.max.mc <- lsat.max[rep == i]

  # summarize NDVI across sites
  lsat.max.mc.smry <- lsat.max.mc[, .(ndvi.med = median(ndvi, na.rm=T), 
                                      ndvi.2yr.med = median(ndvi.2yr.mean, na.rm=T),
                                      ndvi.3yr.med = median(ndvi.3yr.mean, na.rm=T),
                                      ndvi.gs.med = median(ndvi.xcal.gs.med, na.rm=T),
                                      ndvi.q90.med = median(ndvi.xcal.gs.q90, na.rm=T)), 
                                  by=year]
  
  lsat.max.mc.smry <- lsat.max.mc.smry[, ':='(ndvi.lag1.med = lag(ndvi.med), ndvi.lag2.med = lag(ndvi.med, 2), ndvi.2yr.med.lag1 = lag(ndvi.2yr.med,1))]
  
  lsat.max.mc.smry <- lsat.max.mc.smry[year >= 1990]
  
  # detrend select NDVI variables
  lsat.max.mc.smry <- lsat.max.mc.smry[, ndvi.med.dt := lm(ndvi.med ~ year)$resid]
  lsat.max.mc.smry <- lsat.max.mc.smry[, ndvi.3yr.med.dt := lm(ndvi.3yr.med ~ year)$resid]
  lsat.max.mc.smry <- lsat.max.mc.smry[, ndvi.2yr.med.lag1.dt := lm(ndvi.2yr.med.lag1 ~ year)$resid]
  
  # combine npp and lsat data into one datatable
  bio.dt <- lsat.max.mc.smry[npp.mc.smry, on = 'year']
  bio.dt <- bio.dt[, rep := i]
  bio.ts.mc.list[[i]] <- bio.dt
  
  # subsample years used in comparson
  bio.dt <- bio.dt[, .SD[sample(.N,.N*0.9)]]
  
  # recast data table into long format
  bio.long.dt <- melt(bio.dt, id.vars = c('year','npp.gm2yr.med', 'npp.gm2yr.log.med','npp.gm2yr.med.dt','npp.gm2yr.log.med.dt'), value.name = 'ndvi', variable.name = 'ndvi.var')
  bio.long.dt <- melt(bio.long.dt, id.vars = c('year','ndvi.var','ndvi'), value.name = 'npp', variable.name = 'npp.var')
  bio.long.dt <- data.table(bio.long.dt)
  bio.long.dt <- bio.long.dt[ndvi.var != 'rep']
  
  # correlate NDVI variables with ANPP
  bio.cors.dt <- bio.long.dt[, .(rho = cor.test(npp, ndvi, alternative = 'greater', method = 'spearman')$estimate,
                                 pval = cor.test(npp, ndvi, alternative = 'greater', method = 'spearman')$p.val), by = c('ndvi.var','npp.var')]
  bio.cors.dt <- bio.cors.dt[, rep := i]
  bio.cor.mc.list[[i]] <- bio.cors.dt
  
  print(paste0('Finished: ', i/mc.reps))
}

# SUMMARIZE MONTE CARLO ITERATIONS =====================================================================================

# time series for each site
bio.ts.mc.reps.dt <- data.table(rbindlist(bio.ts.mc.list))

bio.ts.mc.smry.dt <- bio.ts.mc.reps.dt[, .(ndvi.med=median(ndvi.med), ndvi.q025=quantile(ndvi.med,0.025), ndvi.q975=quantile(ndvi.med,0.975),
                                           ndvi.lag2.med=median(ndvi.lag2.med), ndvi.lag2.med.q025=quantile(ndvi.lag2.med,0.025), ndvi.lag2.med.q975=quantile(ndvi.lag2.med,0.975),
                                           ndvi.3yr.med=median(ndvi.med), ndvi.3yr.med.q025=quantile(ndvi.med,0.025), ndvi.3yr.med.q975=quantile(ndvi.med,0.975),
                                           ndvi.2yr.med.lag1=median(ndvi.2yr.med.lag1), ndvi.2yr.med.lag1.q025=quantile(ndvi.2yr.med.lag1,0.025), ndvi.2yr.med.lag1.q975=quantile(ndvi.2yr.med.lag1,0.975),
                                           npp.gm2yr.med=median(npp.gm2yr.med), npp.gm2yr.med.q025=quantile(npp.gm2yr.med,0.025), npp.gm2yr.med.q975=quantile(npp.gm2yr.med,0.975)),
                                         by = c('year')]

bio.ts.mc.smry.dt


# correlations- ndvi vs npp
bio.cor.mc.reps.dt <- data.table(rbindlist(bio.cor.mc.list))
bio.cor.mc.reps.dt

bio.cor.mc.smry.dt <- bio.cor.mc.reps.dt[, .(rho=round(median(rho),2), rho.q025=round(quantile(rho,0.025),2), rho.q975=round(quantile(rho,0.975),2),
                                             pval=median(pval), pval.q95=quantile(pval,0.95)),
                                             by = c('ndvi.var','npp.var')]

bio.cor.mc.smry.dt

fwrite(bio.cor.mc.smry.dt, 'output/site_comparisons/lsat_ndvi_cor_with_graminoid_anpp_summary.csv')


# TIME SERIES PLOT and SCATTER PLOT (ANPP VS 3-year MEAN NDVI) =============================================================================
npp.col = 'orange4'
npp.ci.col = adjustcolor('orange', alpha.f = 0.5)
ndvi.col = 'darkgreen'
ndvi.ci.col = adjustcolor('green', alpha.f = 0.5)

my.cex = 1.3
axis.cex=1.7
pch.cex=1.2

npp.lab  <- expression('Graminoid ANPP [g m'^-2*' yr'^-1*']')
ndvi.lab <- expression("Prior 2-yr mean Landsat NDVI"[max]~'[unitless]')

r.coef.dt <- bio.cor.mc.smry.dt[ndvi.var == 'ndvi.2yr.med.lag1' & npp.var == 'npp.gm2yr.med']
r.coef <- paste0(sprintf('%.2f', r.coef.dt$rho), ' [',
                 sprintf('%.2f', r.coef.dt$rho.q025), ', ',
                 sprintf('%.2f', r.coef.dt$rho.q975), ']')
r.lab <- bquote('r'[s]~' = '~.(r.coef))


jpeg('figures/Landsat_NDVI_and_sedgeANPP_1990to2017.jpg', 12, 5, units = 'in', res = 400)
par.op <- par(mfrow=c(1,2), mar=c(5,5.5,1,5.5))

# TIME SERIES PLOTS 
# plot mean sedge anpp time series
plot(npp.gm2yr.med ~ year, bio.ts.mc.smry.dt, type='l', xlim=c(1990,2017), ylim=c(0,70), yaxt='n', xaxt='n', ylab='', xlab='', col=npp.col, lwd=2, pch=16)
axis(1, seq(1990,2017,5), cex.axis = my.cex)
mtext('Year', 1, line = 3, cex = my.cex)
axis(2, seq(0,80,20), las=2, cex.axis = my.cex)
mtext(npp.lab, 2, line = 3.5, cex = my.cex, col=npp.col)
text(1992, 68, '(a)', cex = 2, font=2)
shade.confidence.interval(bio.ts.mc.smry.dt$year, lowerCI = bio.ts.mc.smry.dt$npp.gm2yr.med.q025, upperCI = bio.ts.mc.smry.dt$npp.gm2yr.med.q975, band.col = npp.ci.col)

# plot mean ndvi time series
par(new=T)
plot(ndvi.2yr.med.lag1 ~ year, bio.ts.mc.smry.dt, type='l', xlim=c(1990,2017), ylim=c(0.40,0.65), yaxt='n', xaxt='n', ylab='', xlab='', col=ndvi.col, lwd=2, pch=16)
axis(4, seq(0.40,0.65,0.05), las=2, cex.axis = my.cex)
mtext(ndvi.lab, 4, line = 5, cex = my.cex, col=ndvi.col)
shade.confidence.interval(bio.ts.mc.smry.dt$year, lowerCI = bio.ts.mc.smry.dt$ndvi.2yr.med.lag1.q025, upperCI = bio.ts.mc.smry.dt$ndvi.2yr.med.lag1.q975, band.col = ndvi.ci.col)

# SCATTER PLOT
plotCI(bio.ts.mc.smry.dt$npp.gm2yr.med, bio.ts.mc.smry.dt$ndvi.2yr.med.lag1, li = bio.ts.mc.smry.dt$ndvi.2yr.med.lag1.q025, ui = bio.ts.mc.smry.dt$ndvi.2yr.med.lag1.q975, , err='y', xaxt='n', yaxt='n', 
       xlim=c(0, 75), ylim=c(0.40,0.65),xlab = '', ylab='',cex = axis.cex, pch=16, scol=adjustcolor('grey50', alpha.f = 0.5))
plotCI(bio.ts.mc.smry.dt$npp.gm2yr.med, bio.ts.mc.smry.dt$ndvi.2yr.med.lag1, li = bio.ts.mc.smry.dt$npp.gm2yr.med.q025, ui = bio.ts.mc.smry.dt$npp.gm2yr.med.q975, err='x', xaxt='n', yaxt='n', 
       xlim=c(0,75), ylim=c(0.40,0.65),xlab = '', ylab='',cex = axis.cex, pch=16, scol=adjustcolor('grey50', alpha.f = 0.5), add=T)
axis(1, seq(0,80,20), cex.axis=my.cex)
axis(2, at = seq(0.40, 0.65,0.05), las=2, cex.axis=my.cex, labels = T)
mtext(side = 1, line = 3.0, cex = pch.cex, npp.lab, col=npp.col)
# mtext(side = 2, line = 4.0, cex = pch.cex, ndvi.3yr.lab)
text(55, 0.44, r.lab, cex = my.cex)
text(8, 0.565, '(b)', cex = 2, font=2)

par(par.op)
dev.off()


# number of subsites with landsat data each year 
lsat.max[rep==1][, .N, by = c('site','year')][, .N, by = year]

# END SCRIPT ===========================================================================================================================