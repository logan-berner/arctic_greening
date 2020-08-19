# Prepare Berkely Earth air temperature data set for the Arctic domain
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

abs.stk <- stack('gis_data/climate/berkearth/tavg/Land_and_Ocean_LatLong1.nc', varname = 'climatology') # absolute temperature climatology
anom.stk <- stack('gis_data/climate/berkearth/tavg/Land_and_Ocean_LatLong1.nc', varname = 'temperature') # anomalies from climatology

months <- c('01','02','03','04','05','06','07','08','09','10','11','12')

n.ts <- nlayers(anom.stk)
ts.year <- as.numeric(substr(names(anom.stk), 2,5))
ts.months <- rep(1:12, n.ts / 12)

yrs <- unique(ts.year)
n.yrs <- length(yrs)

plot(abs.stk[[6]])
# output directories
mkdirs('A:/research/data/climate/berkearth/tavg/monthly')
mkdirs('gis_data/climate/berkearth/tavg/monthly')
mkdirs('gis_data/climate/berkearth/swi/monthly/')
mkdirs('gis_data/climate/berkearth/swi/annual/')


# MONTHLY MEAN TEMPERATURE --------------------------------------------------------------------------------------
clust <- makeCluster(8)
registerDoParallel(clust)
foreach(i=1:n.ts) %dopar% {
    # tmp files
    require(raster)
    tmp.dir <- paste('C:/tmp/R_',i, sep='')
    tempfile(tmpdir=tmp.dir)
    rasterOptions(tmpdir=tmp.dir)
    
    # subset data for period of interest
    tavg <- subset(abs.stk, ts.months[i]) 
    tanom <- subset(anom.stk, i) 
    
    # add anomaly and baseline, then scale
    r.global <- (tavg + tanom) * 10 # scale by 10 
    
    # crop / mask to Arctic tundra biome
    r.aoi <- raster::crop(r.global, aoi.wgs84)
    r.aoi.laea <- projectRaster(r.aoi, tmplt.laea.50km, method = 'bilinear')
    r.aoi.laea <- raster::mask(r.aoi.laea, tmplt.laea.50km)
    
    # write arctic rasters 
    r.aoi.laea <- round(r.aoi.laea)
    arctic.outname <- paste('gis_data/climate/berkearth/tavg/monthly/berkearth_tavg_', ts.year[i], '_', months[ts.months[i]],'_degCx10_50km_laea.tif', sep='')
    writeRaster(r.aoi.laea, arctic.outname, overwrite=T, dataType='INT2S')
    
    # clean up
    gc()
    removeTmpFiles()
    unlink(tmp.dir, recursive = T)
  }
stopCluster(clust)


# COMPUTE MONTHLY SUMMER WARMTH INDEX --------------------------------------------------------------------------
swi.outdir <- paste('gis_data/climate/berkearth/swi/monthly', sep='')
mkdirs(swi.outdir)

tmp.files <- list.files('gis_data/climate/berkearth/tavg/monthly/', full.names = T)
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

# COMPUTE ANNUAL SUMMER WARMTH INDEX -----------------------------------------------------------
swi.yrly.outdir <- 'gis_data/climate/berkearth/swi/annual/'
mkdirs(swi.yrly.outdir)

n.yrs <- length(yrs)
swi.monthly.files <- list.files('gis_data/climate/berkearth/swi/monthly/', full.names = T)
swi.annual.outnames <- paste(swi.yrly.outdir, 'berkearth_swi_', yrs, '_degCx10_50km_laea.tif', sep='') 

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
