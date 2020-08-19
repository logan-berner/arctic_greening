# COMPUTE ANNUAL SUMMER WARTH INDEX USING UNIV DELEWARE MONTHLY TEMPERATURE
rm(list=ls())
require(R.matlab)
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
tmplt.wgs84 <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, resolution = 0.5, crs=wgs84)

yrs <- 1900:2017
n.yrs <- length(yrs)
months <- c('01','02','03','04','05','06','07','08','09','10','11','12')


# MONTHLY MEAN TEMPERATURE FROM 1900 - 2014 ------------------------------------------------------------
tavg.files <- list.files('A:/research/data/climate/udel/tavg/ascii/', full.names = T)
n.files <- length(tavg.files)

# output directories
mkdirs('gis_data/climate/udel/tavg/monthly')

for (i in 1:n.yrs){
  clim <- read.table(tavg.files[i])

  clust <- makeCluster(6)
  registerDoParallel(clust)
  foreach(j=1:12) %dopar% {
    # tmp files
    require(raster)
    tmp.dir <- paste('C:/tmp/R_',j, sep='')
    tempfile(tmpdir=tmp.dir)
    rasterOptions(tmpdir=tmp.dir)
    
    # rasterize ascii data and scale
    r.global <- rasterize(x = clim[,1:2], tmplt.wgs84, field = clim[,j+2])
    r.global <- r.global * 10
    
    # crop / mask to Arctic tundra biome
    r.aoi <- raster::crop(r.global, aoi.wgs84)
    r.aoi.laea <- projectRaster(r.aoi, tmplt.laea.50km, method = 'bilinear')
    r.aoi.laea <- raster::mask(r.aoi.laea, tmplt.laea.50km)
    
    # write out arctic rasters 
    r.aoi.laea <- round(r.aoi.laea)
    arctic.outname <- paste('gis_data/climate/udel/tavg/monthly/udel_tavg_', yrs[i], '_', months[j],'_degCx10_50km_laea.tif', sep='')
    writeRaster(r.aoi.laea, arctic.outname, overwrite=T, dataType='INT2S')
    
    # clean up
    gc()
    removeTmpFiles()
    unlink(tmp.dir, recursive = T)
  }
  stopCluster(clust)
  print(i/length(n.yrs))
}


# COMPUTE MONTHLY SUMMER WARMTH INDEX ---------------------------------------------------------
swi.outdir <- paste('gis_data/climate/udel/swi/monthly', sep='')
mkdirs(swi.outdir)

tmp.files <- list.files('gis_data/climate/udel/tavg/monthly/', full.names = T)
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


# COMPUTE ANNUAL SUMMER WARMTH INDEX -------------------------------------------------------------
swi.yrly.outdir <- 'gis_data/climate/udel/swi/annual/'
mkdirs(swi.yrly.outdir)

n.yrs <- length(yrs)
swi.monthly.files <- list.files('gis_data/climate/udel/swi/monthly/', full.names = T)
swi.annual.outnames <- paste(swi.yrly.outdir, 'udel_swi_', yrs, '_degCx10_50km_laea.tif', sep='') 

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