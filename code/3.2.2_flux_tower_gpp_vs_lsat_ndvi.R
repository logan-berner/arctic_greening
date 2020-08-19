# ABOUT THIS SCRIPT =========================================================================================================
# This R script compares mean annual Landsat NDVI with gross primary productivity estimates at flux towers across the Arctic. 
# Author: Logan Berner, NAU
# Date: 2020-03-25

# SET UP WORKSPACE =========================================================================================================
rm(list=ls())
require(lattice)
require(data.table)
require(dplyr)
require(plotrix)

setwd('/projects/arctic/users/lberner/arctic_greening/')

# LOAD DATA SETS =========================================================================================================
lsat.dt <- fread('output/site_comparisons/lsat_ndvi_max_for_flux_tower_sites_mcreps.csv')
flux.dt <- fread('data/field_data/tundra_flux_towers/tundra_flux_tower_gpp_annual_pcntiles.csv')
flux.sites <- fread('data/field_data/tundra_flux_towers/flux_tower_tundra_sites.csv')

# PREP FLUX DATA =====================================================================  
flux.dt <- flux.dt[flux.sites, on = 'site']

# drop some of the flux data
flux.dt <- flux.dt[fluxnet.tier != 2] # only take open-access data 
flux.dt <- flux.dt[gpp.gCm2yr > 0] # GPP can't be negative 
flux.dt <- flux.dt[site != 'FI-Lom'] # forested wetland outside of tundra
flux.dt <- flux.dt[site != 'SJ-Adv'] # no landsat data at this site during flux years
flux.dt <- flux.dt[site != 'SJ-Blv'] # no landsat data at this site during flux years

# compute log GPP 
flux.dt$gpp.gCm2yr.log <- log(flux.dt$gpp.gCm2yr)

# flux site summary table 
site.yrs <- flux.dt[, .(yr.min = min(year), yr.max = max(year), n.yrs = length(unique(year))), by = 'site']
flux.sites <- flux.sites[, c('class','fluxnet.tier') := NULL]
site.meta.dt <- flux.sites[site.yrs, on = 'site']

mean(site.meta.dt$n.yrs)
sd(site.meta.dt$n.yrs)


# START MONTE CARLO ITERATIONS =========================================================================================================
bio.avg.mc.list <- list()
bio.cor.mc.list <- list()

mc.reps <- 1000
for (i in 1:mc.reps){ 
  # get lsat NDVI permutation
  lsat.mc.dt <- lsat.dt[rep == i]
  
  # get flux GPP permutation
  flux.mc.dt <- flux.dt[, .SD[sample(length(5:95),1, replace = T)], by = c('site','year')]

  # combine lsat and flux tower data -----------------
  bio.dt <- lsat.mc.dt[flux.mc.dt, on = c('site','year')]
  bio.dt <- bio.dt[is.na(ndvi) == F]
  
  # number of years of overlap at each site
  bio.dt[, n.yrs := .N, by = 'site']
  
  # subsample years used from each site ----------------------------
  bio.dt <- bio.dt[, .SD[sample(.N,.N*0.9)], by = 'site']

  # compute average ndvi and fluxs across years for each site -------------------------
  bio.avg.dt <- bio.dt[, .(ndvi = median(ndvi, na.rm=T),
                           gpp.gCm2yr.med = median(gpp.gCm2yr, na.rm=T),
                           gpp.gCm2yr.log.med = median(gpp.gCm2yr.log, na.rm=T),
                           n.yrs = median(n.yrs)),
                       by = 'site']

  bio.avg.dt <- bio.avg.dt[, rep := i]
  bio.avg.mc.list[[i]] <- bio.avg.dt
  
  # recast data table into long format
  bio.long.dt <- melt(bio.avg.dt, id.vars = c('site','ndvi','n.yrs'), value.name = 'gpp', variable.name = 'gpp.var')
  bio.long.dt <- data.table(bio.long.dt)
  bio.long.dt <- bio.long.dt[gpp.var != 'rep']
  
  # correlate NDVI variables with ANPP
  bio.cors.dt <- bio.long.dt[, .(rho = cor.test(gpp, ndvi, method = 'spearman')$estimate), by = c('gpp.var')]
  bio.cors.dt <- bio.cors.dt[, rep := i]
  bio.cor.mc.list[[i]] <- bio.cors.dt
  
  print(i/mc.reps)
}

# SUMMARIZE MONTE CARLO ITERATIONS ===================================================================
# correlations- ndvi vs npp
bio.cor.mc.reps.dt <- data.table(rbindlist(bio.cor.mc.list))
bio.cor.mc.smry.dt <- bio.cor.mc.reps.dt[, .(rho=median(rho), rho.q025=quantile(rho,0.025), rho.q975=quantile(rho,0.975)), by = 'gpp.var']

bio.cor.mc.smry.dt
fwrite(bio.cor.mc.smry.dt, 'output/site_comparisons/lsat_ndvi_cor_with_flux_tower_gpp_summary.csv')

# site averages
bio.avg.mc.dt <- data.table(rbindlist(bio.avg.mc.list))

bio.avg.mc.smry.dt <- bio.avg.mc.dt[, .(ndvi=median(ndvi), 
                                        ndvi.q025=quantile(ndvi,0.025),
                                        ndvi.q975=quantile(ndvi,0.975),
                                        gpp.gCm2yr.med=median(gpp.gCm2yr.med),
                                        gpp.gCm2yr.med.q025=quantile(gpp.gCm2yr.med,0.025), 
                                        gpp.gCm2yr.med.q975=quantile(gpp.gCm2yr.med,0.975),
                                        n.yrs = median(n.yrs)),
                                       by = c('site')]

bio.avg.mc.smry.dt


# 'fancy' output table summarizing site
site.meta.dt <- site.meta.dt[bio.avg.mc.smry.dt, on = 'site'] # add sites years and coords
site.meta.dt <- site.meta.dt[order(source,site)] 

site.meta.fancy.dt <- site.meta.dt
site.meta.fancy.dt <- site.meta.fancy.dt[, ':='(lat = round(lat, 3), lon = round(lon, 3),
                                                years = paste0(yr.min,'-',yr.max),
                                                ndvi = paste0(sprintf('%.2f', ndvi),' [', sprintf('%.2f', ndvi.q025),',', sprintf('%.2f', ndvi.q975),']'),
                                                gpp = paste0(sprintf('%.0f',gpp.gCm2yr.med),' [',sprintf('%.0f',gpp.gCm2yr.med.q025),',',sprintf('%.0f',gpp.gCm2yr.med.q975),']'))]

site.meta.fancy.dt[, c('yr.min','yr.max','n.yrs','ndvi.q025','ndvi.q975','gpp.gCm2yr.med','gpp.gCm2yr.med.q025','gpp.gCm2yr.med.q975'):=NULL]
site.meta.fancy.dt <- site.meta.fancy.dt %>% dplyr::select('source', 'site','lat','lon','years','ndvi','gpp') %>% data.table()

fwrite(site.meta.fancy.dt, 'output/site_comparisons/flux_tower_meta_smry.csv')

# check correlation using site averages across all MC iterations
cor.test(bio.avg.mc.smry.dt$ndvi, bio.avg.mc.smry.dt$gpp.gCm2yr.med, alternative = 'greater', method = 'spearman')

# SCATTER PLOT OF LANDSAT NDVI VS GPP ============================================================================================
ylab.lsat <- expression('Landsat NDVI'[max]~'[unitless]')
xlab.gpp <- expression('Ecosystem GPP [g C m'^-2~' yr'^-1~']')

my.pch.cex=1.5
my.cex.axis=1.5

r.coef <- paste0(sprintf('%.2f', bio.cor.mc.smry.dt[gpp.var == 'gpp.gCm2yr.med']$rho), ' [',
                 sprintf('%.2f', bio.cor.mc.smry.dt[gpp.var == 'gpp.gCm2yr.med']$rho.q025), ', ',
                 sprintf('%.2f', bio.cor.mc.smry.dt[gpp.var == 'gpp.gCm2yr.med']$rho.q975), ']')
                
r.lab <- bquote('r'[s]~' = '~.(r.coef))

jpeg('figures/Landsat_NDVI_vs_EcoGPP_avgs.jpg', 5, 5, units = 'in', res = 400)
par.op <- par(mar=c(5,5.5,1,1))
plotCI(bio.avg.mc.smry.dt$gpp.gCm2yr.med, bio.avg.mc.smry.dt$ndvi, li = bio.avg.mc.smry.dt$gpp.gCm2yr.med.q025, ui = bio.avg.mc.smry.dt$gpp.gCm2yr.med.q975, err='x', xaxt='n', yaxt='n', 
       xlim=c(0, 450), ylim=c(0.39,0.75),xlab = '', ylab='', cex = my.cex.axis, pch=16, scol=adjustcolor('grey50', alpha.f = 0.5))
plotCI(bio.avg.mc.smry.dt$gpp.gCm2yr.med, bio.avg.mc.smry.dt$ndvi, li = bio.avg.mc.smry.dt$ndvi.q025, ui = bio.avg.mc.smry.dt$ndvi.q975, err='y', xaxt='n', yaxt='n', 
       xlim=c(0,450), ylim=c(0.39,0.75),xlab = '', ylab='', cex = my.cex.axis, pch=16, scol=adjustcolor('grey50', alpha.f = 0.5), add=T)
axis(1, seq(0,400,200), cex.axis=my.cex.axis)
axis(2, at = seq(0.4, 0.7,0.1), las=2, cex.axis=my.cex.axis, labels = T)
mtext(xlab.gpp, 1, line = 4, cex = my.pch.cex)
mtext(ylab.lsat, 2, line = 4, cex = my.pch.cex)
text(300, 0.40, r.lab, cex = 1.5)
text(bio.avg.mc.smry.dt$gpp.gCm2yr.med, bio.avg.mc.smry.dt$ndvi, bio.avg.mc.smry.dt$n.yrs, col = 'lightblue', cex = 0.5) # sample size 
par(par.op)
dev.off()

# END SCRIPT =========================================================================================================