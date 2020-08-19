# This R script preps several geospatial data sets use in analysis of Landsat NDVI trends across the circumboreal domain.
# Author: Logan Berner, NAU
# Date: 2019-07-15
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
rm(list=ls())
require(raster)
require(data.table)
require(dplyr)
require(tidyr)
require(zyp)
require(rgdal)
require(R.utils)
require(rgeos)
require(gdalUtils)

tmp.dir <- '/scratch/lb968/R_tmp_alt'
tempfile(tmpdir=tmp.dir)
rasterOptions(tmpdir=tmp.dir)

epsg.laea <-  'EPSG:3571'
epsg.polarstereo <- 'EPSG:3995'
setwd('/projects/arctic/users/lberner/arctic_greening/')

source('/home/lb968/code/arctic_greening/0.2_fun_stack_trend.R')

##### LOAD FILES  ============================================================================

# domain files
arctic.laea.shp <- readOGR(dsn = 'data/gis_data/arctic_zones/arctic_oroarctic_laea_buf50km.shp')
extnt.laea <- round(extent(arctic.laea.shp))

highlat.shp.file <- '/projects/arctic/users/lberner/arctic_greening/data/gis_data/arctic_zones/45N.shp'

# active layer thickness
alt.files <- list.files('/projects/arctic/share/geol/esa_cci_active_layer_thickness/0_orig/', full.names = T)
alt.names <- list.files('/projects/arctic/share/geol/esa_cci_active_layer_thickness/0_orig/', full.names = F)

#### CHANGE FORMAT AND CLIP/REPROJECT RASTERS ============================================================================

for (i in 1:length(alt.names)){
  # nc --> geotiff
  tmpfile1 <- paste0('/scratch/lb968/pf_alt_format_',i,'.tif')
  gdal_translate(alt.files[i], tmpfile1)
  
  # project from polar stereo to LAEA
  tmpfile2 <- paste0('/scratch/lb968/pf_alt_proj_',i,'.tif')
  gdalwarp(srcfile = tmpfile1, dstfile = tmpfile2, s_srs = epsg.polarstereo, t_srs = epsg.laea, r = 'bilinear', tr = c(1000,1000), te = extnt.laea[c(1,3,2,4)])
  print('finished projecting')
  
  # load raster into R
  alt.r <- raster(tmpfile2)
  
  # mask to arctic
  alt.r  <- crop(alt.r, arctic.laea.shp)
  alt.r <- mask(alt.r, arctic.laea.shp)
  print('finished masking...')
  
  # write out
  mkdirs('data/gis_data/geol/esa_cci_pf_alt/annual/')
  
  alt.names <- gsub('fv01.0.nc','1km-laea.tif', alt.names)
  outname <- paste0('data/gis_data/geol/esa_cci_pf_alt/annual/', alt.names[i])
  
  writeRaster(alt.r, outname, overwrite = T)
  print(i/length(alt.names))
}

#### COMPUTE MEAN ANNUAL ALT AND TREND IN ALT ============================================================================
alt.files <- list.files('data/gis_data/geol/esa_cci_pf_alt/annual/', full.names = T)
alt.files <- alt.files[-length(alt.files )] # drop 2017
alt.stk <- stack(alt.files)

# mean annual alt
alt.avg.r <- mean(alt.stk, na.rm=T)
writeRaster(alt.avg.r, 'data/gis_data/geol/esa_cci_pf_alt/esa_cci_pf_alt_cm_avg_2003to2016_1km_laea.tif', overwrite=T)

# alt trend
alt.trnd.stk <- stack.trend(alt.stk)
writeRaster(alt.trnd.r, 'data/gis_data/geol/esa_cci_pf_alt/esa_cci_pf_alt_cm_trend_2003to2016_1km_laea.tif', overwrite=T)


#### CLEAN UP ============================================================================
gc()
removeTmpFiles()
unlink(tmp.dir, recursive = T)

#### END SCRIPT  ============================================================================