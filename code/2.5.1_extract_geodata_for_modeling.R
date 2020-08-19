# This R script takes Landsat sampling sites across the arctic tunbra biome and extracts geospatial data for the sites
# Date: 2020-04-23
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(data.table)
require(missForest)
require(reshape2)
require(maptools)
require(raster)
require(rgdal)
require(tidyr)
require(dplyr)
require(zyp)
require(ggplot2)
require(R.utils)
wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
laea <- CRS("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

setwd('/projects/arctic/users/lberner/arctic_greening/')

calc.trends <- function(x,y){
  xx <- zyp.sen(y~x) 
  return(data.table(slope=xx$coefficients[2], int=xx$coefficients[1]))
}


# LOAD AND SPATIALIZE SAMPLING SITES  =====================================================================================
site.dt <- fread('data/lsat_samples/lsat_tundra_samples_n50000.csv')
site.dt <- site.dt[,.(latitude=mean(latitude),longitude=mean(longitude)), by='site']
pts.wgs84 <- SpatialPoints(coords = site.dt[,c(3,2)], proj4string = wgs84)
pts.laea <- spTransform(pts.wgs84, CRSobj = laea)


# EXTRACT ARCTIC ZONE ====================================================================================
zone.shp <- readOGR('data/gis_data/arctic_zones/arctic_three_lat_zones_laea.shp')
site.dt[, lat.zone := sp::over(pts.laea, zone.shp)[,1]]
site.dt[lat.zone == 1, lat.zone := "High Arctic"]
site.dt[lat.zone == 2, lat.zone := "Low Arctic"]
site.dt[lat.zone == 3, lat.zone := "Oro Arctic"]


# EXTRACT SWI TREND  ====================================================================================
swi.trnds.dt <- fread('output/swi_site_trends/lsat_swi_site_trend_summary.csv')
swi.trnds.dt <- swi.trnds.dt[period == '2000-2016']
cols <- c('site','swi.int','swi.change','swi.change.pcnt')
swi.trnds.dt <- swi.trnds.dt[, ..cols]
setnames(swi.trnds.dt, 'swi.int', 'swi.2000')
site.dt <- swi.trnds.dt[site.dt, on = 'site']


# EXTRACT TERRA CLIMATE SOIL MOISTURE  =====================================================

tclim.meta.df <- data.frame(var=c('soil'), scaler=c(1), units = c('mm'), fun = c('min'))

tclim.files <- list.files('data/gis_data/climate/terra_clim/seasonal/', full.names = T)
tclim.names <- list.files('data/gis_data/climate/terra_clim/seasonal/', full.names = F)

seasons <- c('djf','mam','jja','son')

# extract data
for (i in 1:nrow(tclim.meta.df)){
  for (j in seasons){
    in.names <- unique(grep(paste(tclim.meta.df$var[i],collapse="|"), tclim.names, value=TRUE)) # subset time period
    in.files <- unique(grep(paste(tclim.meta.df$var[i],collapse="|"), tclim.files, value=TRUE)) # subset time period
    in.files <- in.files[grep(j, in.names)]
    stk <- stack(in.files)
    clim.dt <- data.table(raster::extract(stk, pts.laea))
    names(clim.dt) <- paste0('x', 2000:2016)
    clim.dt <- clim.dt[, site := site.dt$site]
    clim.dt <- melt.data.table(clim.dt, id.vars = c('site'), variable.name = 'year', value.name = 'clim')
    clim.dt <- clim.dt[is.na(clim)==F]
    clim.dt <- clim.dt[, n.yrs := .N, by = site]
    clim.dt[, year := as.numeric(gsub('x','', year))]
    clim.dt[, year.rsc := year - min(year)]
    clim.dt[, clim := clim * as.numeric(tclim.meta.df$scaler[i])]
    
    # compute trend
    clim.trnd.dt <- clim.dt %>% group_by(site) %>% do(out=calc.trends(x=.$year.rsc, y=.$clim)) %>% unnest(cols=c(out)) %>% data.table()
    clim.trnd.dt[, abs.chng := slope * 17, by = 'site']
    clim.trnd.dt <- clim.trnd.dt[, c('slope') := NULL]
    var.name <- paste0(tclim.meta.df$var[i],'.',j, '.', tclim.meta.df$fun[i])
    clim.trnd.dt <- setnames(clim.trnd.dt, 
                             c('int','abs.chng'),
                             c(paste(var.name,'2000',tclim.meta.df$units[i], sep='.'), paste(var.name, 'change', tclim.meta.df$units[i], sep='.')))
    
    # store 
    site.dt <- clim.trnd.dt[site.dt, on = 'site']
    print(paste(i/nrow(tclim.meta.df), j))
  }
}


# EXTRACT MODIS BURN YEAR  ===================================================================================
burnyr.r <- raster('data/gis_data/fires/arctic_modis_v6_burn_year_2000to2016_500m_laea_5x5majfilt.tif')
site.dt[, burnyr := raster::extract(burnyr.r, pts.laea)] # extract raster values straight into the main data table
rm(burnyr.r)
   

# EXTRACT THERMOKARST VULNERABILITY ==========================================================================
thermokarst.shp <- readOGR('data/gis_data/geol/thermokarst/arctic_thermokarst_landscapes_laea.shp')
site.dt[, tk_wetland := sp::over(pts.laea, thermokarst.shp)[,1]]
site.dt[, tk_lake := sp::over(pts.laea, thermokarst.shp)[,2]]
site.dt[, tk_hill := sp::over(pts.laea, thermokarst.shp)[,3]]

# set some blank missings values to NA
site.dt[tk_hill == '', tk_hill := NA]
site.dt[tk_wetland == '', tk_wetland := NA]
site.dt[tk_lake == '', tk_lake := NA]
rm(thermokarst.shp)


# EXTRACT PERMAFROST ACTIVE LAYER THICKNESS ==========================================================================
alt.files <- list.files('data/gis_data/geol/esa_cci_pf_alt/annual/', full.names = T)[-15]
alt.stk <- stack(alt.files)
names(alt.stk) <- 2003:2016

# extract and reformat
alt.dt <- data.table(raster::extract(alt.stk, pts.laea))
alt.dt <- alt.dt[, site := site.dt$site]
alt.dt <- melt.data.table(alt.dt, id.vars = c('site'), variable.name = 'year', value.name = 'alt')
alt.dt <- alt.dt[is.na(alt)==F]
alt.dt <- alt.dt[, n.yrs := .N, by = site][n.yrs == 14]
alt.dt[, year := as.numeric(as.character(gsub("X","",alt.dt$year)))][, year.rsc := year - 2003]

# compute trend in ALT
alt.trnd.dt <- alt.dt %>% group_by(site) %>% do(out=calc.trends(x=.$year.rsc, y=.$alt)) %>% unnest(cols=c(out)) %>% data.table()
alt.trnd.dt[, alt.change.cm := slope * 14, by = 'site']
alt.trnd.dt <- setnames(alt.trnd.dt, 'int','alt.2003.cm')
alt.trnd.dt <- alt.trnd.dt[, c('slope') := NULL]
fivenum(alt.trnd.dt$alt.change.pcnt)

# store 
site.dt <- alt.trnd.dt[site.dt, on = 'site']
rm(alt.stk); rm(alt.dt)


# EXTRACT PERMAFROST GROUND SURFACE TEMPERATURE ==========================================================================
gst.files <- list.files('data/gis_data/geol/esa_cci_pf_gst/annual/', full.names = T)[-15]
gst.stk <- stack(gst.files)
names(gst.stk) <- 2003:2016

# extract and reformat
gst.dt <- data.table(raster::extract(gst.stk, pts.laea))
gst.dt <- gst.dt[, site := site.dt$site]
gst.dt <- melt.data.table(gst.dt, id.vars = c('site'), variable.name = 'year', value.name = 'gst')
gst.dt <- gst.dt[is.na(gst)==F]
gst.dt <- gst.dt[, n.yrs := .N, by = site][n.yrs == 14]
gst.dt[, year := as.numeric(as.character(gsub("X","",gst.dt$year)))][, year.rsc := year - 2003]

# compute trend in gst
gst.trnd.dt <- gst.dt %>% group_by(site) %>% do(out=calc.trends(x=.$year.rsc, y=.$gst)) %>% unnest(cols=c(out)) %>% data.table()
gst.trnd.dt[, gst.change.degC := slope * 14, by = 'site']
gst.trnd.dt <- setnames(gst.trnd.dt, 'int','gst.2003.degC')
gst.trnd.dt <- gst.trnd.dt[, c('slope') := NULL]
fivenum(gst.trnd.dt$gst.change.degC)

# store 
site.dt <- gst.trnd.dt[site.dt, on = 'site']
rm(gst.stk); rm(gst.dt)


# EXTRACT PERMAFROST PROBABILITY ==========================================================================
pfe.files <- list.files('data/gis_data/geol/esa_cci_pf_pfe/annual/', full.names = T)[-15]
pfe.stk <- stack(pfe.files)
names(pfe.stk) <- 2003:2016

# extract and reformat
pfe.dt <- data.table(raster::extract(pfe.stk, pts.laea))
pfe.dt <- pfe.dt[, site := site.dt$site]
pfe.dt <- melt.data.table(pfe.dt, id.vars = c('site'), variable.name = 'year', value.name = 'pfe')
pfe.dt <- pfe.dt[is.na(pfe)==F]
pfe.dt <- pfe.dt[, n.yrs := .N, by = site][n.yrs == 14]
pfe.dt[, year := as.numeric(as.character(gsub("X","",pfe.dt$year)))][, year.rsc := year - 2003]

# compute trend in pfe
pfe.trnd.dt <- pfe.dt %>% group_by(site) %>% do(out=calc.trends(x=.$year.rsc, y=.$pfe)) %>% unnest(cols=c(out)) %>% data.table()
pfe.trnd.dt[, pfe.change.pcnt := slope * 14, by = 'site']
pfe.trnd.dt <- setnames(pfe.trnd.dt, 'int','pfe.2003.pcnt')
pfe.trnd.dt <- pfe.trnd.dt[, c('slope') := NULL]

# store 
site.dt <- pfe.trnd.dt[site.dt, on = 'site']
rm(pfe.stk); rm(pfe.dt)


# EXTRACT TANDEM X TOPOGRAPHY VARIABLES =====================================================
topo.files <- data.frame(file = list.files('data/gis_data/topo/', full.names = T),
                         name = list.files('data/gis_data/topo/', full.names = F))
topo.files$name <- tolower(gsub('TDM1_DEM__30_', '', gsub('_clip.tif','', topo.files$name)))
topo.files$name <- gsub('dem','elev', topo.files$name)

# extract topo variables
for (i in 1:nrow(topo.files)){
  r <- raster(as.character(topo.files$file[i]))
  site.dt[, topo.files$name[i] := raster::extract(r, pts.laea)]
  print(i/nrow(topo.files))
}
site.dt[elev < 0, elev := 0]

# cap tpi at 5% - 95% of range and tri at 0 - 95% of range
tpi.hilo <- quantile(na.omit(site.dt$tpi), probs = c(0.05,0.95))
tri.hilo <- quantile(na.omit(site.dt$tri), probs = 0.95)
site.dt[tpi <= tpi.hilo[1], tpi := tpi.hilo[1]]
site.dt[tpi >= tpi.hilo[2], tpi := tpi.hilo[2]]
site.dt[tri >= tri.hilo[1], tri := tri.hilo[1]]


# EXTRACT LAND COVER VFC VARIABLES =====================================================
# land cover and forest type keys
landcov.key <- fread('data/gis_data/landcover/probaV_vcf/probav_landcov_key.csv')
fortype.key <- data.table(forest.type=c(0:5), forest.class=c('unknown','enf','ebf','dnf','dfb','mix'))

# get files 
landcov.files <- data.frame(file = list.files('data/gis_data/landcover/probaV_vcf/', full.names = T),
                            name = list.files('data/gis_data/landcover/probaV_vcf/', full.names = F))
landcov.files <- landcov.files[7, ]

# set up layer names
landcov.files$name <- gsub('ProbaV_LC100_epoch2015_global_v2.0.2_','', landcov.files$name)
landcov.files$name <- gsub('-layer_EPSG-4326_reproj_LAEA_clip_arctic.tif','', landcov.files$name)
landcov.files$name <- gsub('discrete-classification_EPSG-4326_reproj_LAEA_clip_arctic.tif','landcov', landcov.files$name)

# extract land cover
r <- raster(as.character(landcov.files$file))
site.dt[, landcov.files$name := raster::extract(r, pts.laea)]
  
# replace land cover code with class name 
landcov.key <- setnames(landcov.key, 'code','landcov')
site.dt <- landcov.key[site.dt, on = 'landcov']
site.dt <- site.dt[, landcov := NULL]
site.dt <- setnames(site.dt, 'class','landcov')

# WRITE OUT CSV WITH ALL THE DATA ==============================================================================================================
mkdirs('output/lsat_site_trend_modeling')
fwrite(site.dt, 'output/lsat_site_trend_modeling/tundra_site_biogeoclim.csv')

# END SCRIPT ======================================================================================================