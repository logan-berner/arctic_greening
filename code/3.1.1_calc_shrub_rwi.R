# ABOUT THIS SCRIPT ================================================================================================
# This R script computes shrub growth indices from ring-width measurements
# AUTHOR: LOGAN BERNER, NAU
# DATE: 2020-04-13
# SET UP WORKSPACE ================================================================================================
rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.6.2'))
require(dplyr)
require(data.table)
require(reshape2)
require(dplR)
require(R.utils)

args <- commandArgs(TRUE)
i = as.numeric(args[1])
i = 1

# reformat i for sorting
if (i < 10){
  i <- paste0('000',i)
} else if (i < 100){
  i <- paste0('00',i)
} else if (i < 1000){
  i <- paste0('0',i)
}

setwd('/projects/arctic/users/lberner/arctic_greening/')

# LOAD SHRUB RING WIDTH DATA ====================================================================
shrub.rw.dt <- fread('data/field_data/tundra_shrub_growth/tundra_shrub_ring_width.csv')
shrub.rw.dt

# DERIVE RING-WIDTH INDICES  ====================================================================

# compute mean ring width (for plant where two radii were measured)
shrub.rw.dt <- shrub.rw.dt[, .(ring.width.mm.avg = mean(ring.width.mm, na.rm = T)), by=c('cluster', 'site', 'plant', 'genus', 'species', 'year', 'contributor', 'shrub.hub')]

# set missing rings to have rw of 0.01 mm
shrub.rw.dt <- shrub.rw.dt[ring.width.mm.avg == 0, ring.width.mm.avg:= 0.01]

# determine plant age at each site step, dealing with the first 5 years being trimmed off the shrub hub data set 
shrub.rw.dt[shrub.hub == 'yes', age := seq_along(ring.width.mm.avg)+5, by = 'plant']
shrub.rw.dt[shrub.hub == 'no', age := seq_along(ring.width.mm.avg), by = 'plant']
shrub.rw.dt <- shrub.rw.dt[age > 5]

# compute rwi using spline to detrend and standardize series
shrub.rw.dt[, spl.fit := suppressWarnings(ffcsaps(ring.width.mm.avg, year, nyrs = runif(1, 15, 25), f = runif(1, 0.45, 0.55))), by = 'plant']
shrub.rw.dt[, rwi := ring.width.mm.avg / spl.fit]


# COMPUTE MEAN INTER-SERIES CORRELATIONS FOR EACH SHRUB  ====================================================================
clst.rbar.dt <- data.table(cluster=NA, genus=NA, plant=NA, res.cor=NA, p.val=NA)

clusters <- unique(shrub.rw.dt$cluster)
genera <- unique(shrub.rw.dt$genus)
for (j in clusters){
  for (k in genera){
    clst <- shrub.rw.dt %>% filter(cluster == j & genus == k)
    if (nrow(clst)==0){
      next()
    } else {
      clst.wide <- dcast(clst, year ~ plant, value.var = 'rwi')
      rownames(clst.wide) <- clst.wide$year
      clst.wide <- clst.wide[,-1]
      
      # shrub-wise interseries correlations
      clst.rbar <- interseries.cor(clst.wide, method='spearman')
      clst.rbar$plant <- rownames(clst.rbar)
      clst.rbar$cluster <- j
      clst.rbar$genus <- k
      clst.rbar.dt <- rbind(clst.rbar.dt, clst.rbar)
    }
  }
}

clst.rbar.dt <- data.table(na.omit(clst.rbar.dt))
clst.rbar.smry.dt <- clst.rbar.dt[, .(rbar.avg = mean(res.cor)), by = c('cluster', 'genus')]


# BUILD CHRONOLOGY =====================================================================================
# subsample shrubs used in each chronology
clst.plants.dt <- shrub.rw.dt[, .(plant = unique(plant)), by = 'cluster']
clst.plants.smpl.dt <- clst.plants.dt[, .SD[sample(.N, .N*0.9)], by = 'cluster']
shrub.rw.dt <- shrub.rw.dt[plant %in% clst.plants.smpl.dt$plant]
  
# compute median chron for each cluster
crn.dt <- shrub.rw.dt[, .(rwi.med = median(rwi, na.rm=T), n.shrubs = .N), by = c('cluster','genus','shrub.hub','year')]
  
# add rbar statistic 
crn.dt <- clst.rbar.smry.dt[crn.dt, on = c('cluster','genus')]

# store MC rep
crn.dt$rep <- i
  
# write out 
mkdirs('output/site_comparisons/shrub_chrons_mcreps/')
fwrite(crn.dt, paste0('output/site_comparisons/shrub_chrons_mcreps/tundra_shrub_cluster_chron_timeseries_rep_',i,'.csv'))

# END SCRIPT =====================================================================================
print('All done!!')