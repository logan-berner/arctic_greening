# ABOUT THIS SCRIPT =====================================================================================================
# This R scirpt compiles together flux tower data from sites in the Arctic and outputs monthly and annual GPP.
# Author: Logan Berner, NAU
# Date: 2019-12-06

# SET UP WORKSPACE =====================================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(tidyr)
require(lattice)
require(data.table)
require(rgeos)
require(maptools)
require(raster)
# require(rgdal)
require(lubridate)
require(chron)
require(dplyr)

setwd('/projects/above_gedi/lberner/arctic_greening/data/')

# FLUXNET: COMBINE ALL ANNUAL SITE FILES INTO ONE DATA TABLE =========================================================================================
files.annual <- list.files('field_data/tundra_flux_towers/unzip/', pattern = 'FULLSET_YY_', full.names = T)
n.files.annual <- length(files.annual)

sites <- list.files('field_data/tundra_flux_towers/unzip/', pattern = 'FULLSET_MM_', full.names = F)
sites <- substr(sites, 5, 10)

for (i in 1:n.files.annual){
  if (i == 1){
    fluxnet.gpp.yrly.dt <- fread(files.annual[i])
    fluxnet.gpp.yrly.dt$site <- sites[i]
    # cnames <- colnames(fluxnet.gpp.yrly.dt)
  } else {
    tmp.ts <- fread(files.annual[i])
    # colnames(tmp.ts) <- cnames # some of the sample output files.annual are missing the column names
    tmp.ts$site <- sites[i]
    fluxnet.gpp.yrly.dt <- rbind(fluxnet.gpp.yrly.dt, tmp.ts, fill = T)
  }
  print(i/n.files.annual)
}

colnames(fluxnet.gpp.yrly.dt) <- tolower(colnames(fluxnet.gpp.yrly.dt))

fluxnet.gpp.yrly.dt$year <- fluxnet.gpp.yrly.dt$timestamp

# select GPP variables
fluxnet.gpp.yrly.dt[, ':='(gpp.gCm2yr.q05 = gpp_nt_vut_05,
                           gpp.gCm2yr.q15 = gpp_nt_vut_16,
                           gpp.gCm2yr.q25 = gpp_nt_vut_25,
                           gpp.gCm2yr.q50 = gpp_nt_vut_50,
                           gpp.gCm2yr.q75 = gpp_nt_vut_75,
                           gpp.gCm2yr.q84 = gpp_nt_vut_84,
                           gpp.gCm2yr.q95 = gpp_nt_vut_95)] 

cols.to.keep <- c('site','year','gpp.gCm2yr.q05', 'gpp.gCm2yr.q15', 'gpp.gCm2yr.q25', 
                  'gpp.gCm2yr.q50', 'gpp.gCm2yr.q75', 'gpp.gCm2yr.q84', 'gpp.gCm2yr.q95')

fluxnet.gpp.yrly.dt <- fluxnet.gpp.yrly.dt[, ..cols.to.keep]

fluxnet.gpp.yrly.dt <- fluxnet.gpp.yrly.dt[gpp.gCm2yr.q50 != -9999]

# interpolate yrly GPP between uncertainty percentiles
fluxnet.gpp.yrly.dt <- melt(fluxnet.gpp.yrly.dt, id.vars = c('site','year'), variable.name = 'percentile', value.name = 'gpp.gCm2yr')
fluxnet.gpp.yrly.dt[, percentile := as.numeric(gsub('gpp.gCm2yr.q','', percentile))]
fluxnet.gpp.yrly.dt <- fluxnet.gpp.yrly.dt[, .(percentile = approx(percentile, gpp.gCm2yr, xout = 5:95)$x, 
                                               gpp.gCm2yr = approx(percentile, gpp.gCm2yr, xout = 5:95)$y), by = c('site','year')]

# compute ratio between 50th percentile and all other percentiles
fluxnet.gpp.yrly.p50.dt <- fluxnet.gpp.yrly.dt[percentile == 50][, percentile := NULL]
setnames(fluxnet.gpp.yrly.p50.dt, 'gpp.gCm2yr','gpp.gCm2yr.p50')
fluxnet.gpp.yrly.dt <- fluxnet.gpp.yrly.p50.dt[fluxnet.gpp.yrly.dt, on = c('site','year')]
fluxnet.gpp.yrly.dt[, uncertainty.frac := gpp.gCm2yr/gpp.gCm2yr.p50][, gpp.gCm2yr.p50 := NULL]

fluxnet.gpp.med.uncert.dt <- fluxnet.gpp.yrly.dt[, .(uncertainty.frac = median(uncertainty.frac)), by = c('site','percentile')][, .(uncertainty.frac = median(uncertainty.frac)), by = 'percentile']


# AON: COMBINE ALL 30 MIN MEASUREMENTS INTO ONE FILE =====================================================================================================
aon.30m.files <- list.files('field_data/tundra_flux_towers/aon/', pattern = '_gapfilled', full.names = T)
n.aon.30m.files <- length(aon.30m.files)
aon.sites <- substr(aon.30m.files, 44,47)

aon.gpp.30m.dt <- data.table(Year=NA,DoY=NA,Hour=NA,GPP=NA,GPP.f=NA, GPP.gC.m2.30min=NA, GPP.f.gC.m2.30min=NA)
for (i in 1:n.aon.30m.files){
  tmp.ts <- fread(aon.30m.files[i])
  tmp.ts <- tmp.ts[-1,] # drop row of units
  tmp.ts$GPP <- as.numeric(tmp.ts$GPP)
  tmp.ts$GPP.f <- as.numeric(tmp.ts$GPP_f)
  tmp.ts <- tmp.ts[GPP != -9999][GPP.f != -9999]
  tmp.ts <- tmp.ts[, c('Year','DoY','Hour','GPP','GPP.f')]
  tmp.ts$site <- aon.sites[i]
  tmp.ts$GPP.gC.m2.30min <- tmp.ts$GPP * (10^-6) * (44.0095) * (12.07 / 44.0095) * (1800)
  tmp.ts$GPP.f.gC.m2.30min <- tmp.ts$GPP.f * (10^-6) * (44.0095) * (12.07 / 44.0095) * (1800)
  aon.gpp.30m.dt <- rbind(aon.gpp.30m.dt, tmp.ts, fill = T)
  print(i/n.aon.30m.files)
  print(dim(aon.gpp.30m.dt))
}
names(aon.gpp.30m.dt) <- tolower(names(aon.gpp.30m.dt))
aon.gpp.30m.dt <- aon.gpp.30m.dt[-1,] # drop first row of NAs 
aon.gpp.30m.dt[, year := as.numeric(as.character(aon.gpp.30m.dt$year))]
aon.gpp.30m.dt[, doy := as.numeric(as.character(aon.gpp.30m.dt$doy))]
aon.gpp.30m.dt[, month := month.day.year(aon.gpp.30m.dt$doy)$month] # determine month from julian day

# aon.gpp.30m.dt <- aon.gpp.30m.dt[year != 2018]
aon.gpp.30m.dt[site == 1523, site := 'IVO-FEN']
aon.gpp.30m.dt[site == 1991, site := 'IVO-RID']
aon.gpp.30m.dt[site == 1993, site := 'IVO-TUS']
aon.gpp.30m.dt[site == 2044, site := 'CHE']

# compute monthly GPP
aon.gpp.monthly.dt <- aon.gpp.30m.dt[, .(gpp.gCm2mon = sum(gpp.f.gc.m2.30min, na.rm = T)), by=c('site', 'year', 'month')]

# compute annual GPP
aon.gpp.yrly.dt <- aon.gpp.30m.dt[, .(gpp.gCm2yr = sum(gpp.f.gc.m2.30min, na.rm = T)), by=c('site', 'year')]
aon.gpp.yrly.dt <- na.omit(aon.gpp.yrly.dt)
aon.gpp.yrly.dt <- aon.gpp.yrly.dt[order(site,year)]

# add uncertainty estimates derived from fluxnet
full.fac <- aon.gpp.yrly.dt %>% expand(site, year, percentile = 5:95)
full.fac <- data.table(full.fac)

aon.gpp.yrly.dt <- aon.gpp.yrly.dt[full.fac, on = c('site','year')]
aon.gpp.yrly.dt <- na.omit(aon.gpp.yrly.dt)

aon.gpp.yrly.dt <- aon.gpp.yrly.dt[fluxnet.gpp.med.uncert.dt, on = 'percentile']
aon.gpp.yrly.dt <- aon.gpp.yrly.dt[order(site,year,percentile)]

aon.gpp.yrly.dt <- aon.gpp.yrly.dt[, gpp.gCm2yr := gpp.gCm2yr * uncertainty.frac]

# COMBINE FLUXNET AND AON  =====================================================================================================

# annual fluxes
tundra.gpp.yrly <- rbind(fluxnet.gpp.yrly.dt, aon.gpp.yrly.dt)
tundra.gpp.yrly <- tundra.gpp.yrly[gpp.gCm2yr != -9999]

# xyplot(gpp.gCm2Yr ~ year | site, tundra.gpp.yrly)

write.table(tundra.gpp.yrly, 'field_data/tundra_flux_towers/tundra_flux_tower_gpp_annual_pcntiles.csv', sep = ',', row.names = F, col.names = T)

# END SCRIPT  =====================================================================================================