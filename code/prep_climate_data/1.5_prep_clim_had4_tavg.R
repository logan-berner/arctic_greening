# COMPUTE ANNUAL SUMMER WARTH INDEX USING THE "HadCRUT4 hybrid with UAH" MONTHLY TEMPERATURE ANOMALIES
rm(list=ls())
require(ncdf4)
require(raster)
require(maptools)
require(R.utils)
require(foreach)
require(doParallel)

setwd('C:/Users/Logan/Google Drive/research/nau/nsf_arctic/arctic_greening/data/')

wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
laea <- CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

aoi.wgs84 <- readShapePoly('gis_data/arctic_oroarctic_wgs84_buf50km.shp', proj4string = wgs84)
aoi.laea <- readShapePoly('gis_data/arctic_oroarctic_zones_laea.shp', proj4string = laea)
tmplt.laea.50km <- raster('gis_data/arctic_oroarctic_laea_buf50km_50km.tif')

anom.stk <- stack('gis_data/climate/hadcru4/tavg/had4_short_uah_v2_0_0.nc', varname = 'temperature_anomaly') # anomalies from 1981-2010 climatology
plot(anom.stk[[50]])

months <- c('01','02','03','04','05','06','07','08','09','10','11','12')

n.ts <- nlayers(anom.stk)
ts.year <- as.numeric(substr(names(anom.stk), 2,5))
ts.months <- rep(1:12, n.ts / 12)

yrs <- unique(ts.year)
n.yrs <- length(yrs)


# output directories
mkdirs('A:/research/data/climate/hadcru4/tavg/monthly')
mkdirs('gis_data/climate/hadcru4/tavg/monthly')
mkdirs('gis_data/climate/hadcru4/swi/monthly/')
mkdirs('gis_data/climate/hadcru4/swi/annual/')


# CREATE CLIMATOLOGY (1981-2010) USING BERKELEY EARTH, CRU, and UDEL DATA SETS ----------------------------------------------
climotology.yrs <- 1981:2010
berk.files <- list.files('gis_data/climate/berkearth/tavg/monthly/', full.names = T)
cru.files <- list.files('gis_data/climate/cru_ts4.01/tavg/monthly/', full.names = T)
udel.files <- list.files('gis_data/climate/udel/tavg/monthly/', full.names = T)

for (i in 1:12){ # for each month
  # get files for the month of interest
  yr.mons <- paste(climotology.yrs, months[i], sep="_")
  berk.stk <- stack(unique(grep(paste(yr.mons,collapse="|"), berk.files, value=TRUE)))
  cru.stk <- stack(unique(grep(paste(yr.mons,collapse="|"), cru.files, value=TRUE)))
  udel.stk <- stack(unique(grep(paste(yr.mons,collapse="|"), udel.files, value=TRUE)))
  
  # compute average across years for each data sets
  berk.stk.avg <- mean(berk.stk, na.rm=T)
  cru.stk.avg <- mean(cru.stk, na.rm=T)
  udel.stk.avg <- mean(udel.stk, na.rm=T)
  
  # combine data sets
  ens.stk <- stack(berk.stk.avg, cru.stk.avg, udel.stk.avg)
  ens.avg <- stackApply(ens.stk, rep(nlayers(ens.stk),1), fun = mean, na.rm=T)
  
  # store ensemble avg climatology
  if (i == 1){
    abs.stk <- ens.avg 
  } else {
    abs.stk <- addLayer(abs.stk, ens.avg)
  }
  print(i)
}


# MONTHLY MEAN TEMPERATURE FROM 1979 - 2017 -------------------------------------------
clust <- makeCluster(8)
registerDoParallel(clust)
foreach(i=1:n.ts) %dopar% {
    # tmp files
    require(raster)
    tmp.dir <- paste('C:/tmp/R_',i, sep='')
    tempfile(tmpdir=tmp.dir)
    rasterOptions(tmpdir=tmp.dir)
    
    # get baseline from ensemble (already multiplied by 10) 
    tavg.aoi.laea <- subset(abs.stk, ts.months[i])
        
    # get anom for time period and scale
    tanom <- subset(anom.stk, i) * 10 
    
    # crop / mask to Arctic tundra biome
    tanom.aoi <- raster::crop(tanom, aoi.wgs84)
    tanom.aoi.laea <- projectRaster(tanom.aoi, tmplt.laea.50km, method = 'bilinear')
    tanom.aoi.laea <- raster::mask(tanom.aoi.laea, tmplt.laea.50km)
    
    # compute absolute temperature from anomaly + monthly climatology
    r.aoi.laea <- tavg.aoi.laea + tanom.aoi.laea 
    
    # write out
    r.aoi.laea <- round(r.aoi.laea)
    arctic.outname <- paste('gis_data/climate/hadcru4/tavg/monthly/hadcru4_tavg_', ts.year[i], '_', months[ts.months[i]],'_degCx10_50km_laea.tif', sep='')
    writeRaster(r.aoi.laea, arctic.outname, overwrite=T, dataType='INT2S')
    
    # clean up
    gc()
    removeTmpFiles()
    unlink(tmp.dir, recursive = T)
  }
stopCluster(clust)


# COMPUTE MONTHLY SUMMER WARMTH INDEX ----------------------------------------------------------
swi.outdir <- paste('gis_data/climate/hadcru4/swi/monthly', sep='')
mkdirs(swi.outdir)

tmp.files <- list.files('gis_data/climate/hadcru4/tavg/monthly/', full.names = T)
swi.outnames <- gsub('tavg','swi', tmp.files)

clust <- makeCluster(8)
registerDoParallel(clust)
foreach(i=1:length(tmp.files)) %dopar% {
  # tmp files
  require(raster)
  tmp.dir <- paste('C:/tmp/R_',i, sep='')
  tempfile(tmpdir=tmp.dir)
  rasterOptions(tmpdir=tmp.dir)
  
  # swi computation
  r <- raster(tmp.files[i])
  r[r<0] <- 0
  writeRaster(r, swi.outnames[i], overwrite=T, dataType='INT2U')

  # clean up
  gc()
  removeTmpFiles()
  unlink(tmp.dir, recursive = T)
}
stopCluster(clust)


# COMPUTE ANNUAL SUMMER WARMTH INDEX ------------------------------------------------------------------------------
swi.yrly.outdir <- 'gis_data/climate/hadcru4/swi/annual/'
mkdirs(swi.yrly.outdir)

n.yrs <- length(yrs)
swi.monthly.files <- list.files('gis_data/climate/hadcru4/swi/monthly/', full.names = T)
swi.annual.outnames <- paste(swi.yrly.outdir, 'hadcru4_swi_', yrs, '_degCx10_50km_laea.tif', sep='') 

clust <- makeCluster(8)
registerDoParallel(clust)
foreach(i=1:n.yrs) %dopar% {
  # tmp files
  require(raster)
  tmp.dir <- paste('C:/tmp/R_',i, sep='')
  tempfile(tmpdir=tmp.dir)
  rasterOptions(tmpdir=tmp.dir)
  
  yr.files <- swi.monthly.files[grep(yrs[i], swi.monthly.files)]
  stk <- stack(yr.files)
  stk.sum <- sum(stk)
  writeRaster(stk.sum, swi.annual.outnames[i], overwrite=T, dataType='INT2U')
  
  # clean up
  gc()
  removeTmpFiles()
  unlink(tmp.dir, recursive = T)
}
stopCluster(clust)