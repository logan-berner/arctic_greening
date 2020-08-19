# This R script takes Landsat sampling sites across the arctic tunbra biome and extracts geospatial data for the sites
# Date: 2020-04-23
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(data.table)
require(reshape2)
require(maptools)
require(raster)
require(rgdal)
require(tidyr)
wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
laea <- CRS("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

setwd('/projects/arctic/users/lberner/arctic_greening/')

# LOAD DATA SETS =====================================================================================
site.dt <- fread('data/lsat_samples/lsat_tundra_samples_n50000.csv')

# arctic zones and keys
arctic.lat.zones.laea <- readOGR('data/gis_data/arctic_zones/arctic_three_lat_zones_laea.shp')
arctic.lon.zones.laea <- readOGR('data/gis_data/arctic_zones/arctic_oroarctic_lon_zones_laea.shp')
arctic.countries.laea <- readOGR('data/gis_data/arctic_zones/arctic_oroarctic_countries_laea.shp')

# create raster with ID from each pixel in climate data 
clim.pxl.id.r <- raster('data/gis_data/climate/ensemble/ensemble_common_mask.tif')
values(clim.pxl.id.r) <- 1:ncell(clim.pxl.id.r)

# swi annual files from all models 
clim.yrs <- 1982:2016
swi.models <- c('berkearth','cru_ts4.01','giss','hadcru4','udel')
swi.file.lst <- list()

for (i in swi.models){
  mdl <- data.table(file = list.files(paste0('data/gis_data/climate/',i,'/swi/annual'), full.names = T, pattern = glob2rx('*.tif')),
                    name = list.files(paste0('data/gis_data/climate/',i,'/swi/annual'), full.names = F, pattern = glob2rx('*.tif')))
  mdl$col.name <- paste('swi',i, sep='.')
  mdl <- mdl[grep(paste(clim.yrs, collapse="|"), mdl$name, value=F)] # get specific years
  mdl$year <- clim.yrs
  swi.file.lst[[i]] <- mdl
}
swi.file.dt <- rbindlist(swi.file.lst)



# BUILD OUTPUT TABLES ====================================================================================
site.dt <- site.dt[,.(latitude=mean(latitude),longitude=mean(longitude)), by='site'] # store static data in this data table
full.fac <- site.dt %>% expand(site, year = 1982:2016)
full.fac <- data.table(full.fac)
site.yr.dt <- site.dt[full.fac, on = c('site')]
site.yr.dt <- site.yr.dt[order(site,year)] # store time series data in this data table

# SPATIALIZE SITES ====================================================================================
pts.wgs84 <- SpatialPoints(coords = site.dt[,c(3,2)], proj4string = wgs84)
pts.laea <- spTransform(pts.wgs84, CRSobj = laea)

# double check alignment against zones and climate data
# plot(arctic.lat.zones.laea, col = rainbow(5))
# points(pts.laea, col = 'white', pch = '*', cex = 0.5)


# EXTRACT STATIC GEOSPATIAL DATA FOR SITES  ============================================================================================================== 

# lat / zone subzones
site.dt$country <- over(pts.laea, arctic.countries.laea)[,1]
site.dt$lat.zone <- over(pts.laea, arctic.lat.zones.laea)[,2]
site.dt$lon.zone <- over(pts.laea, arctic.lon.zones.laea)[,2]

# climate grid cell id 
site.dt$clim.pxl <- as.numeric(raster::extract(clim.pxl.id.r, pts.laea))

# topo
site.dt$elev.m <- raster::extract(elev.r, pts.laea)
site.dt$slope.deg <- raster::extract(slope.r, pts.laea)
site.dt$southness <- raster::extract(southness.r, pts.laea)
site.dt$westness <- raster::extract(westness.r, pts.laea)


# EXTRACT CLIMATE TIME SERIES GEOSPATIAL DATA FOR SITES  ============================================================================================================== 

## extract annual SWI
for (j in 1:nrow(swi.file.dt)){
  r <- raster(as.character(swi.file.dt$file[j]))
  site.yr.dt[year==swi.file.dt$year[j], swi.file.dt$col.name[j] := raster::extract(r, pts.laea) / 10]
  print(j/nrow(swi.file.dt))
}
site.yr.dt$swi.ensemble <- rowMeans(site.yr.dt[, 5:9], na.rm = T)


# WRITE OUT CSV WITH ALL THE DATA ==============================================================================================================
fwrite(site.dt, 'output/tundra_site_conditions.csv')
fwrite(site.yr.dt, 'output/tundra_site_climate_timeseries.csv')
