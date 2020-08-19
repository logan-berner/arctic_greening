# This R script takes takes global Terra Climate data and clips / reprojects to high latitudes
# Author: Logan Berner
# Date: 2020-05-14
#---------------------------------------------------------------------------------------------------------------------
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.6.2/'))
require(R.utils)
require(raster)
require(ncdf4)
require(rgdal)
wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
laea <- CRS("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
setwd('/projects/arctic/')

# SET UP TEMPORARY DIRECTORY FOR RASTER  #---------------------------------------------------------------------------
tmp.dir <- paste('/scratch/lb968/Rtmp/tclim')
mkdirs(tmp.dir)
rasterOptions(tmpdir = tmp.dir)

# METADATA ABOUT VARIABLE #---------------------------------------------------------------------------
var.meta.df <- data.frame(var=c('soil'), units = c('mm'), bitdepth = c('INT2U'), fun = c('min'))

# IDENTIFY INPUT CLIMATE FILE AND CREATE MAIN OUTPUT DIR #-----------------------------

# get input climate file and parse variable and year
yoi <- 1999:2016
in.files <- list.files('share/clim/terra_climate/0_orig/', full.names = T)
in.files <- unique(grep(paste(yoi,collapse="|"), in.files, value=TRUE)) # subset time period
in.files <- grep('soil',in.files, value=T)

file.names <- list.files('share/clim/terra_climate/0_orig/', full.names = F)
file.names <- unique(grep(paste(yoi,collapse="|"), file.names, value=TRUE)) # subset time period
file.names <- grep('soil', file.names, value=T)

# load highlat and domain shapefiles
land.45n.wgs84.shp <- readOGR('users/lberner/arctic_greening/data/gis_data/arctic_zones/45N.shp')
arctic.laea.shp <- readOGR('users/lberner/arctic_greening/data/gis_data/arctic_zones/arctic_oroarctic_laea_buf50km.shp')

# create main output directory
mkdirs('users/lberner/arctic_greening/data/gis_data/climate/terra_clim/monthly')

# LOOP THROUGHT EACH ANNUAL CLIMATE FILE, OUTPUT MONTHLY LAYER AS A GEOTIF, THEN PROJECT TO LAEA, CLIPPING TO HIGH LATITUDES --------------
month.nums <- c('01','02','03','04','05','06','07','08','09','10','11','12')    

print('starting loop...')

# loop through the June - August bands in each file
for (i in yoi){}
  for(j in 1:12){
    in.file <- grep(i, in.files, value=T)
    r <- suppressWarnings(raster::raster(in.file, band = j))
    crs(r) <- wgs84
    r.45n <- raster::crop(r, land.45n.wgs84.shp)
    r.45n.laea <- raster::projectRaster(r.45n, res = 4000, crs = laea, method = 'bilinear')
    r.45n.laea.arctic <- raster::crop(r.45n.laea, arctic.laea.shp) 
    r.45n.laea.arctic <- raster::mask(r.45n.laea.arctic, arctic.laea.shp) 
    r.45n.laea.arctic <- round(r.45n.laea.arctic) # ensure everything is an integer
    out.file <- paste('users/lberner/arctic_greening/data/gis_data/climate/terra_clim/monthly/terraclim_arctic_soil_',i,'_',month.nums[j],'_', var.meta.df$units, '_laea_4km.tif', sep='')
    raster::writeRaster(r.45n.laea.arctic, out.file, datatype = as.character(var.meta.df$bitdepth), overwrite = T)
    print(j)
  }
}
# delete tmp directory
unlink(tmp.dir, recursive = T)
print('done')


# LOOP THROUGHT MONTHLY FILES AND GENERATE SEASONAL COMPOSITES --------------------------------------------------------------
yoi <- 2000:2016
mkdirs('seasonal')

# get files for specific variable
var.files <- list.files('users/lberner/arctic_greening/data/gis_data/climate/terra_clim/monthly/', pattern = 'soil', full.names = T)

season.files <- var.files[-c(1:11, length(var.files)-1)]
season.df <- data.frame(year = sort(rep(yoi,4)), 
                        season = c('djf','mam','jja','son'), 
                        file.index = seq(1, length(yoi)*12, 3))

n.yr.seas <- length(yoi) * 4
print('starting loop...')

# loop through the monthly bands in each file
for(j in 1:n.yr.seas){
  yr <- season.df$year[j]
  season <- season.df$season[j]
  in.files <- season.files[c(season.df$file.index[j]:(season.df$file.index[j]+2))]
  stk <- stack(in.files)
  seas.r <- stackApply(stk, rep(1,3), fun = as.character(var.meta.df$fun))
  seas.r <- mask(seas.r, arctic.laea.shp)
  out.file <- paste('seasonal/terraclim_arctic_', var,'_',yr,'_',season,'_', var.meta.df$fun, '_', var.meta.df$units, '_laea_4km.tif', sep='')
  writeRaster(seas.r, out.file, datatype = as.character(var.meta.df$bitdepth), overwrite = T)
  print(j)
}

# delete tmp directory
unlink(tmp.dir, recursive = T)
print('done')

# END SCRIPT #----------------------------------------------------------------------------------------------------------------------