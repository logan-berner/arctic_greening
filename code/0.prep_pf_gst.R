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

tmp.dir <- paste0('/scratch/lb968/R_tmp_pft_gst')
tempfile(tmpdir=tmp.dir)
rasterOptions(tmpdir=tmp.dir)

epsg.laea <-  'EPSG:3571'
epsg.polarstereo <- 'EPSG:3995'
# setwd('/projects/arctic/users/lberner/arctic_greening/')

source('/home/lb968/code/arctic_greening/0.2_fun_stack_trend.R')

# ##### LOAD FILES  ============================================================================
# 
# # domain files 
# arctic.laea.shp <- readOGR(dsn = 'data/gis_data/arctic_zones/arctic_oroarctic_laea_buf50km.shp')
# extnt.laea <- round(extent(arctic.laea.shp))
# 
# highlat.shp.file <- '/projects/arctic/users/lberner/arctic_greening/data/gis_data/arctic_zones/45N.shp'
# 
# # active layer thickness
# gst.files <- list.files('/projects/arctic/share/geol/esa_cci_ground_temperature/0_orig/', full.names = T)
# gst.names <- list.files('/projects/arctic/share/geol/esa_cci_ground_temperature/0_orig/', full.names = F)
# gst.names <- gsub('fv01.0.nc','1mDepth-degC-1km-laea.tif', gst.names)
# 
# # check bands
# gdalinfo(gst.files[1])
# 
# 
# #### CHANGE FORMAT AND CLIP/REPROJECT RASTERS ============================================================================
# 
# for (i in 1:length(gst.names)){
#   # nc --> geotiff
#   tmpfile1 <- paste0('/scratch/lb968/pf_gst_format_',i,'.tif')
#   gdal_translate(gst.files[i], tmpfile1, sd_index = 1) # ground surface temperature
#   # gdal_translate(gst.files[i], tmpfile1, sd_index = 2) # ground temperature at 1 m
#   
#   # project from polar stereo to LAEA
#   tmpfile2 <- paste0('/scratch/lb968/pf_gst_proj_',i,'.tif')
#   gdalwarp(srcfile = tmpfile1, dstfile = tmpfile2, s_srs = epsg.polarstereo, t_srs = epsg.laea, r = 'bilinear', tr = c(1000,1000), te = extnt.laea[c(1,3,2,4)])
#   print('finished projecting')
#   
#   # load raster into R
#   gst.r <- raster(tmpfile2)
#   
#   # mask to arctic
#   gst.r  <- crop(gst.r, arctic.laea.shp)
#   gst.r <- mask(gst.r, arctic.laea.shp)
#   print('finished masking...')
#   
#   # rescale and convert deg K to
#   gst.r <- (gst.r/100)-272.15
#   
#   # write out
#   mkdirs('data/gis_data/geol/esa_cci_pf_gst/annual/')
#   outname <- paste0('data/gis_data/geol/esa_cci_pf_gst/annual/', gst.names[i])
#   writeRaster(gst.r, outname, overwrite = T)
# 
#   print(i/length(gst.names))
# }
# 
# #### COMPUTE MEAN ANNUAL GST AND TREND IN ALT ============================================================================
gst.files <- list.files('/projects/arctic/users/lberner/arctic_greening/data/gis_data/geol/esa_cci_pf_gst/annual/', full.names = T)
gst.files <- gst.files[-length(gst.files )] # drop 2017
gst.stk <- stack(gst.files)
# 
# # mean annual gst
# gst.avg.r <- mean(gst.stk, na.rm=T)
# writeRaster(gst.avg.r, 'data/gis_data/geol/esa_cci_pf_gst/esa_cci_pf_gst_degC_avg_2003to2016_1km_laea.tif', overwrite=T)

# gst trend
gst.trnd.stk <- stack.trend(gst.stk)
writeRaster(gst.trnd.stk, '/projects/arctic/users/lberner/arctic_greening/data/gis_data/geol/esa_cci_pf_gst/esa_cci_pf_gst_1m_degC_trend_2003to2016_1km_laea.tif', overwrite=T)

#### CLEAN UP ============================================================================
gc()
removeTmpFiles()
unlink(tmp.dir, recursive = T)

#### END SCRIPT  ============================================================================