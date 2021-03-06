# This R script preps several geospatial data sets use in analysis of Landsat NDVI trends across the circumboreal domain.
# Author: Logan Berner, NAU
# Date: 2019-07-15
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
rm(list=ls())
require(raster)
require(data.table)
require(rgdal)
require(R.utils)
require(rgeos)
require(gdalUtils)

# args <- commandArgs(TRUE)
# i = as.numeric(args[1])
# i = 1

tmp.dir <- paste0('/scratch/lb968/R_tmp_pfe')
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
pfe.files <- list.files('/projects/arctic/share/geol/esa_cci_permafrost_extent/0_orig/', full.names = T)
pfe.names <- list.files('/projects/arctic/share/geol/esa_cci_permafrost_extent/0_orig/', full.names = F)

#### CHANGE FORMAT AND CLIP/REPROJECT RASTERS ============================================================================

for (i in 1:length(pfe.names)){
  # nc --> geotiff
  tmpfile1 <- paste0('/scratch/lb968/pf_pfe_format_',i,'.tif')
  gdal_translate(pfe.files[i], tmpfile1)
  
  # project from polar stereo to LAEA
  tmpfile2 <- paste0('/scratch/lb968/pf_pfe_proj_',i,'.tif')
  gdalwarp(srcfile = tmpfile1, dstfile = tmpfile2, s_srs = epsg.polarstereo, t_srs = epsg.laea, r = 'bilinear', tr = c(1000,1000), te = extnt.laea[c(1,3,2,4)])
  print('finished projecting')
  
  # load raster into R
  pfe.r <- raster(tmpfile2)
  
  # mask to arctic
  pfe.r  <- crop(pfe.r, arctic.laea.shp)
  pfe.r <- mask(pfe.r, arctic.laea.shp)
  print('finished masking...')
  
  # write out
  mkdirs('data/gis_data/geol/esa_cci_pf_pfe/annual/')
  
  pfe.names <- gsub('fv01.0.nc','1km-laea.tif', pfe.names)
  outname <- paste0('data/gis_data/geol/esa_cci_pf_pfe/annual/', pfe.names[i])
  
  writeRaster(pfe.r, outname, overwrite = T)
  print(i/length(pfe.names))
}

#### COMPUTE MEAN ANNUAL PFE AND TREND IN PFE ============================================================================
pfe.files <- list.files('data/gis_data/geol/esa_cci_pf_pfe/annual/', full.names = T)
pfe.files <- pfe.files[-length(pfe.files )] # drop 2017
pfe.stk <- stack(pfe.files)

# mean annual pfe
pfe.avg.r <- mean(pfe.stk, na.rm=T)
writeRaster(pfe.avg.r, 'data/gis_data/geol/esa_cci_pf_pfe/esa_cci_pf_pfe_pcnt_avg_2003to2016_1km_laea.tif', overwrite=T)

# pfe trend
pfe.trnd.stk <- stack.trend(pfe.stk)
writeRaster(pfe.trnd.r, 'data/gis_data/geol/esa_cci_pf_pfe/esa_cci_pf_pfe_pcnt_trend_2003to2016_1km_laea.tif', overwrite=T)

#### CLEAN UP ============================================================================
gc()
removeTmpFiles()
unlink(tmp.dir, recursive = T)

#### END SCRIPT  ============================================================================ 