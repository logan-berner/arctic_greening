# This R script generates computes trend in Arctic summer temperatures (SWI) at multiple spatial scales in a Monte Carlo framework.
# Date: 2020-04-06
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(data.table)
require(reshape2)
require(dplyr)
require(tidyr)
require(raster)
require(rgdal)
require(R.utils)
require(maptools)
require(zyp)

args <- commandArgs(TRUE)
i = as.numeric(args[1])
# i = 1

setwd('/projects/arctic/users/lberner/arctic_greening/')
source('/home/lb968/code/arctic_greening/0.2_fun_stack_trend.R')

tmp.dir <- paste('/scratch/lb968/tmp/R_',i, sep='')
tempfile(tmpdir=tmp.dir)
rasterOptions(tmpdir=tmp.dir)

calc.trends <- function(x,y){
  xx <- zyp.yuepilon(y,x) ## note the order of x and y are switched in this call!!!
  return(data.table(slope=round(xx['trend'],5), int=round(xx['intercept'],5), tau=round(xx['tau'],3), pval=round(xx['sig'],4)))
}

# reformat i....
if (i < 10){
  i <- paste0('000',i)
} else if (i < 100){
  i <- paste0('00',i)
} else if (i < 1000){
  i <- paste0('0',i)
}

# LOAD DATA SETS ==============================================================================================

# specify time periods
yoi <- 1985:2016
yoi2 <- 2000:2016
n.yoi <- length(yoi)

# specify projections
wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
laea <- CRS("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# load arctic zones and site coorditanes
arctic.lat.zones <- raster('data/gis_data/arctic_zones/arctic_oroarctic_lat_zones_laea_50km.tif')
site.cords.dt <- fread('output/tundra_site_conditions.csv', fill=T)

# load raster stack
swi.files <- list.files('output/swi_ensemble_grids/mc_reps/', full.names = T)
swi.files <- grep(paste0('rep_',i), swi.files, swi.files, value=T) # files for specific MC rep
swi.stk <- stack(swi.files) / 10 # rescale to deg C months
names(swi.stk) <- 1982:2016


# DUMP DATA INTO DATA TABLE  =================================================================================
swi.dt <- data.table(values(swi.stk))
swi.dt$pxl.id <- 1:nrow(swi.dt)
swi.dt$lat.zone <- values(arctic.lat.zones)
swi.dt <- melt.data.table(swi.dt, id.vars = c('pxl.id','lat.zone'), variable.name = 'year', value.name = 'swi')
swi.dt <- na.omit(swi.dt)
swi.dt$year <- as.numeric(as.character(gsub("X","",swi.dt$year)))

swi.dt[, lat.zone := as.character(lat.zone)]
swi.dt[lat.zone == 1 | lat.zone == 2, lat.zone := "High Arctic"]
swi.dt[lat.zone == 3 | lat.zone == 4, lat.zone := "Low Arctic"]
swi.dt[lat.zone == 5, lat.zone := "Oro Arctic"]


# COMPUTE SWI TRENDS FOR EACH GRID CELL ======================================================================================

# break into two periods and recombine
swi.pxl.yrly.get1985.dt <- swi.dt[year >= 1985][, ':='(period = '1985-2016', year.rsc = year - 1985)]
swi.pxl.yrly.get2000.dt <- swi.dt[year >= 2000][, ':='(period = '2000-2016', year.rsc = year - 2000)]
swi.pxl.yrly.periods.dt <- rbind(swi.pxl.yrly.get1985.dt, swi.pxl.yrly.get2000.dt)

# fit and summarize trend models
swi.pxl.trnd.periods.dt <- swi.pxl.yrly.periods.dt %>% group_by(pxl.id, period) %>% do(out=calc.trends(x=.$year.rsc, y=.$swi)) %>% unnest(cols=c(out)) %>% data.table()
swi.pxl.trnd.periods.dt <- swi.pxl.trnd.periods.dt[order(pxl.id,period)]

# code rasters with trend slopes
swi.trnd.gte1985.r <- arctic.lat.zones
swi.trnd.gte1985.r[] <- NA
swi.trnd.gte1985.r[swi.pxl.trnd.periods.dt[period == '1985-2016']$pxl.id] <- swi.pxl.trnd.periods.dt[period == '1985-2016']$slope

swi.trnd.gte2000.r <- arctic.lat.zones
swi.trnd.gte2000.r[] <- NA
swi.trnd.gte2000.r[swi.pxl.trnd.periods.dt[period == '2000-2016']$pxl.id] <- swi.pxl.trnd.periods.dt[period == '2000-2016']$slope

# make output directories and write out rasters
mkdirs('output/swi_gridded_trends/mc_reps/1985to2016')
mkdirs('output/swi_gridded_trends/mc_reps/2000to2016')
writeRaster(swi.trnd.gte1985.r, overwrite=T, paste0('output/swi_gridded_trends/mc_reps/1985to2016/arctic_tundra_ensemble_swi_degC_trend_1985to2016_rep_',i,'.tif'))
writeRaster(swi.trnd.gte2000.r, overwrite=T, paste0('output/swi_gridded_trends/mc_reps/2000to2016/arctic_tundra_ensemble_swi_degC_trend_2000to2016_rep_',i,'.tif'))


# COMPUTE AVERAGE ZONAL / BIOME SWI TIME SERIES =================================================================================
swi.zonal.yrly.dt <- swi.dt[, .(swi.avg = mean(swi)), by = c('lat.zone','year')]
swi.biome.yrly.dt <- swi.dt[, .(swi.avg = mean(swi)), by = c('year')][, lat.zone := 'Arctic']
swi.zonal.yrly.dt <- rbind(swi.zonal.yrly.dt, swi.biome.yrly.dt)

# write out
swi.zonal.yrly.dt$rep <- i
mkdirs('output/swi_zonal_avg_timeseries/mc_reps')
fwrite(swi.zonal.yrly.dt, paste0('output/swi_zonal_avg_timeseries/mc_reps/swi_ensemble_zonal_avg_timeseries_1985to2016_rep_',i,'.csv'))


# COMPUTE ZONAL TRENDS IN SWI FROM 1985/2000 TO 2016 ===============================================================================
# break data into two periods
swi.zonal.yrly.gte1985.dt <- swi.zonal.yrly.dt[year >= 1985][, ':='(period = '1985-2016', year.rsc = year - 1985)]
swi.zonal.yrly.gte2000.dt <- swi.zonal.yrly.dt[year >= 2000][, ':='(period = '2000-2016', year.rsc = year - 2000)]
swi.zonal.yrly.periods.dt <- rbind(swi.zonal.yrly.gte1985.dt, swi.zonal.yrly.gte2000.dt)

# fit and summarize trend models
swi.trnd.zone.avg.smry.dt <- swi.zonal.yrly.periods.dt %>% group_by(period, lat.zone) %>% do(out=calc.trends(x=.$year.rsc, y=.$swi.avg)) %>% unnest(cols=c(out)) %>% data.table()
swi.trnd.zone.avg.smry.dt <- swi.trnd.zone.avg.smry.dt[, n.yrs :=  c(rep(32,4), rep(17,4))]
swi.trnd.zone.avg.smry.dt <- swi.trnd.zone.avg.smry.dt[, total.change := round(slope * n.yrs, 1)]
swi.trnd.zone.avg.smry.dt <- swi.trnd.zone.avg.smry.dt[, total.change.pcnt := round(total.change / int * 100, 1)]

# write out zonal average trends 
swi.trnd.zone.avg.smry.dt$rep <- i
mkdirs('output/swi_zonal_trends/mc_reps')
fwrite(swi.trnd.zone.avg.smry.dt, paste0('output/swi_zonal_trends/mc_reps/swi_ensemble_zonal_avg_trends_rep_',i,'.csv'))


# COMPUTE SITE LEVELS TRENDS IN SWI FROM 1985/2000 TO 2016 ==============================================================================================

# spatialize sites
pts.wgs84 <- SpatialPoints(coords = site.cords.dt[,c(3,2)], proj4string = wgs84)
pts.laea <- spTransform(pts.wgs84, CRSobj = laea)

# extract climate data for each site
swi.site.dt <- raster::extract(swi.stk, pts.laea)
swi.site.dt <- data.table(swi.site.dt)
swi.site.dt$site <- site.cords.dt$site
swi.site.dt$lat.zone <- site.cords.dt$lat.zone
swi.site.dt[, lat.zone := as.character(lat.zone)]
swi.site.dt[lat.zone == 'AB' | lat.zone == 'C', lat.zone := "High Arctic"]
swi.site.dt[lat.zone == 'D' | lat.zone == 'E', lat.zone := "Low Arctic"]
swi.site.dt[lat.zone == 'O', lat.zone := "Oro Arctic"]

swi.site.dt <- melt.data.table(swi.site.dt, id.vars = c('site','lat.zone'), variable.name = 'year', value.name = 'swi')
swi.site.dt$year <- as.numeric(as.character(gsub("X","",swi.site.dt$year)))

# write out site data 
mkdirs('output/swi_site_timeseries/mc_reps/')
swi.site.dt$rep <- i
fwrite(swi.site.dt, paste0('output/swi_site_timeseries/mc_reps/tundra_site_swi_timeseries_rep_',i,'.csv'))

# split into periods and recombine
swi.site.gte1985.dt <- swi.site.dt[year >= 1985][, ':='(period = '1985-2016', year.rsc = year - 1985)]
swi.site.gte2000.dt <- swi.site.dt[year >= 2000][, ':='(period = '2000-2016', year.rsc = year - 2000)]
swi.site.periods.dt <- rbind(swi.site.gte1985.dt, swi.site.gte2000.dt)

# compute site trends
swi.sites.trnds.dt <- swi.site.periods.dt %>% group_by(period, site) %>% do(out=calc.trends(x=.$year.rsc, y=.$swi)) %>% unnest(cols=c(out)) %>% data.table()
swi.sites.trnds.dt <- swi.sites.trnds.dt[, n.yrs :=  0][period == '1985-2016', n.yrs := 32][period == '2000-2016', n.yrs := 17]
swi.sites.trnds.dt <- swi.sites.trnds.dt[, total.change := round(slope * n.yrs, 1)]
swi.sites.trnds.dt <- swi.sites.trnds.dt[, total.change.pcnt := round(total.change / int * 100, 1)]
swi.sites.trnds.dt <- swi.sites.trnds.dt[, sig := cut(pval, c(-Inf, 0.050, 0.10, Inf), c('sig.p5','sig.p10','insig'))]
swi.sites.trnds.dt <- swi.sites.trnds.dt[, slope.cat := cut(slope, c(-Inf, 0, Inf), c('cooling','warming'))]
swi.sites.trnds.dt <- swi.sites.trnds.dt[, swi.trend.cat := paste(slope.cat, sig, sep = '.')]
swi.sites.trnds.dt <- swi.sites.trnds.dt[sig == 'insig', swi.trend.cat := 'insig']

# write out swi trends for each sampling site
mkdirs('output/swi_site_trends/mc_reps/')
swi.sites.trnds.dt$rep <- i
fwrite(swi.sites.trnds.dt, paste0('output/swi_site_trends/mc_reps/tundra_site_swi_trends_rep_',i,'.csv'))


# CLEAN UP ==============================================================================================
gc()
removeTmpFiles()
unlink(tmp.dir, recursive = T)
print('All done!!')

# END SCRIPT ======================================================================================