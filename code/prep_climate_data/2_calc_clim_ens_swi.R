#-------------------------------------------------------------------------------------------------------------------
# CREATE ANNUAL SUMMER WARTH INDEX FOR ARCTIC USING ENSEMBLE OF CLIMATE DATA SETS FOR 1984 - 2016
# AUTHOR: LOGAN BERNER, NAU
# DATE: 2018-06-05
#-------------------------------------------------------------------------------------------------------------------
rm(list=ls())
require(raster)
require(R.utils)
require(maptools)
require(foreach)
require(doParallel)
# setwd('C:/Users/lb968/Google Drive/research/')
setwd('C:/Users/Logan/Google Drive/research/nau/nsf_arctic/arctic_greening/data/gis_data/')
mkdirs('climate/ensemble/swi/annual/')

laea <- CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
aoi.laea <- readShapePoly('arctic_oroarctic_laea', proj4string = laea)

yrs <- 1980:2016
n.yrs <- length(yrs)

berk.files <- list.files('climate/berkearth/swi/annual/', full.names = T, pattern = glob2rx('*.tif'))
cru.files <- list.files('climate/cru_ts4.01/swi/annual/', full.names = T, pattern = glob2rx('*.tif'))
giss.files <- list.files('climate/giss/swi/annual/', full.names = T, pattern = glob2rx('*.tif'))
hadcru.files <- list.files('climate/hadcru4/swi/annual/', full.names = T, pattern = glob2rx('*.tif'))
udel.files <- list.files('climate/udel/swi/annual/', full.names = T, pattern = glob2rx('*.tif'))

#--------------------------------------------------------------------------------------------------------------------
# DERIVE COMMON MASK FOR ALL OF THE CLIMATE DATA SETS
#--------------------------------------------------------------------------------------------------------------------
berk.r <- raster(berk.files[1], value = T)
cru.r <- raster(grep(yrs[1], cru.files, value = T))
giss.r <- raster(grep(yrs[1], giss.files, value = T))
hadcru.r <- raster(grep(yrs[1], hadcru.files, value = T))
udel.r <- raster(grep(yrs[1], udel.files, value = T))

clim.msk <- berk.r * cru.r * giss.r * hadcru.r * udel.r
clim.msk[clim.msk >= 0] <- 1

plot(clim.msk)
plot(aoi.laea, add=T)

writeRaster(clim.msk, 'climate/ensemble/ensemble_common_mask.tif', overwrite=T)

#--------------------------------------------------------------------------------------------------------------------
# COMPUTE ENSEMBLE MEDIAN FOR EACH YEAR (SPATIAL DATA SETS) 
#--------------------------------------------------------------------------------------------------------------------
clust <- makeCluster(6)
registerDoParallel(clust)
foreach(i=1:n.yrs) %dopar% {
  # tmp files
  require(raster)
  tmp.dir <- paste('C:/tmp/R_',i, sep='')
  tempfile(tmpdir=tmp.dir)
  rasterOptions(tmpdir=tmp.dir)
  
  # load in climate data for given year
  berk.r <- raster(grep(yrs[i], berk.files, value = T))
  cru.r <- raster(grep(yrs[i], cru.files, value = T))
  giss.r <- raster(grep(yrs[i], giss.files, value = T))
  hadcru.r <- raster(grep(yrs[i], hadcru.files, value = T))
  udel.r <- raster(grep(yrs[i], udel.files, value = T))
  
  # create ensemble stack
  ens.stk <- stack(berk.r, cru.r, giss.r, hadcru.r, udel.r)
  ens.stk <- mask(ens.stk, clim.msk)
  
  # compute ensemble median
  ens.med <- stackApply(ens.stk, rep(1, nlayers(ens.stk)), fun = mean, na.rm = F)
  
  outname <- paste('climate/ensemble/swi/annual/ensemble_avg_swi_', yrs[i], '_degCx10_50km_laea.tif', sep='')
  writeRaster(ens.med, outname, overwrite = T, dataType = 'INT2U')

  # clean up
  gc()
  removeTmpFiles()
  unlink(tmp.dir, recursive = T)
}
stopCluster(clust)

