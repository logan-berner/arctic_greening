# This R script generates 
# Date: 2020-04-06
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(data.table)
require(reshape2)
require(raster)
require(rgdal)
require(R.utils)

args <- commandArgs(TRUE)
i = as.numeric(args[1])
# i = 1

setwd('/projects/above_gedi/lberner/arctic_greening/')

# LOAD DATA SETS =====================================================================================
mc.reps <- 1000

# common mask
arctic.tmplt.r <- raster('data/gis_data/climate/ensemble/ensemble_common_mask.tif')

# create raster with ID from each pixel in climate data 
clim.pxl.id.r <- arctic.tmplt.r
values(clim.pxl.id.r) <- 1:ncell(clim.pxl.id.r)

# generate list of swi files from all models 
clim.yrs <- 1982:2016
clim.yr.mc <- clim.yrs[i]
swi.models <- c('berkearth','cru_ts4.01','giss','hadcru4','udel')
swi.file.lst <- list()

for (j in swi.models){
  mdl <- data.table(file = list.files(paste0('data/gis_data/climate/',j,'/swi/annual'), full.names = T, pattern = glob2rx('*.tif')),
                    name = list.files(paste0('data/gis_data/climate/',j,'/swi/annual'), full.names = F, pattern = glob2rx('*.tif')))
  mdl$col.name <- paste('swi',i, sep='.')
  mdl <- mdl[grep(paste(clim.yrs, collapse="|"), mdl$name, value=F)] # get specific years
  mdl$year <- clim.yrs
  swi.file.lst[[j]] <- mdl
}
swi.file.dt <- rbindlist(swi.file.lst)

# subset files associated with specific MC iteration
swi.file.dt <- swi.file.dt[year == clim.yr.mc]


# GENERATE RANDOMIZED RASTERS =====================================================================================

# load raster stack, dump values to data table, cast in a long-format
swi.stk <- stack(arctic.tmplt.r, clim.pxl.id.r, stack(swi.file.dt$file))
swi.dt <- data.table(values(swi.stk))
names(swi.dt) <- c('common.mask','pxl.id', swi.models)
swi.dt <- melt(swi.dt, id.vars = c('common.mask','pxl.id'), variable.name = 'dataset', value.name = 'swi')
swi.dt <- swi.dt[common.mask == 1]

# repeatedly generate annual SWI rasters by randomly selecting data from a data set for each pixel
mkdirs('output/swi_ensemble_grids/mc_reps')

for (k in 1:mc.reps){
  swi.mc.r <- arctic.tmplt.r
  swi.mc.dt <- swi.dt[, .SD[sample(.N, 1)], by = 'pxl.id']
  swi.mc.r[swi.mc.dt$pxl.id] <- swi.mc.dt$swi
  
  # reformat i for naming output
  if (k < 10){
    k <- paste0('000',k)
  } else if (k < 100){
    k <- paste0('00',k)
  } else if (k < 1000){
    k <- paste0('0',k)
  }
  
  outname <- paste0('output/swi_ensemble_grids/mc_reps/ensemble_swi_', clim.yr.mc, '_degCx10_rep_', k, '.tif')
  writeRaster(swi.mc.r, outname, overwrite=T)
}

print('All done!!')