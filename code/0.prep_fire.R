# This R script preps several geospatial data sets use in analysis of Landsat NDVI trends across the circumboreal domain.
# Author: Logan Berner, NAU
# Date: 2019-07-15
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
rm(list=ls())
require(raster)
require(data.table)
require(rgdal)
require(R.utils)
require(gdalUtils)

tmp.dir <- '/scratch/lb968/R_tmp_fire'
tempfile(tmpdir=tmp.dir)
rasterOptions(tmpdir=tmp.dir)

setwd('/projects/arctic/users/lberner/arctic_greening/data/gis_data/')

r <- raster('fires/arctic_modis_v6_burn_year_2000to2016_500m_laea.tif')


r <- raster('climate/worldclim/arctic_worldclim_bio_13_1km_laea.tif')