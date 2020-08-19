# ABOUT THIS SCRIPT =========================================================================
# This R script takes Landsat data downloaded in batch from GEE, applies QAQC steps, and combines
# data into one large spreadsheet. 
# Date: 2018-10-17
# Author: Logan Berner, Northern Arizona University

# SET UP WORKSPACE ============================================================================================
rm(list=ls())
require(tidyr)
require(lattice)
require(data.table)
require(R.utils)
setwd('/projects/arctic/users/lberner/arctic_greening/')
source('/home/lb968/code/arctic_greening/0.1_fun_lsat_tools.R')


# CLEAN BATCHES OF LANDSAT DATA DOWNLOADED FROM GOOGLE EARTH ENGINE ============================================================================================
mkdirs('data/lsat_samples/cleaned/')

files <- list.files('data/lsat_samples/raw/lsat_tundra_samples/', full.names = T)
n.files <- length(files)

# parse file names to get reginon and sample
samples <- list.files('C:/tmp/lsat_tundra_samples/', full.names = F)
samples <- gsub('.csv','',gsub('lsat_tundra_samples_','',samples))
samples <- data.frame(matrix(unlist(strsplit(samples, '_')), ncol = 2, byrow = T))
colnames(samples) <- c('region','sample')
samples$sample <- gsub('srs','',samples$sample)

# build output table in which to store summary stats of samples
sample.smry <- data.frame(region=NA, sample=NA, grab=NA, n.sites.all=NA, n.obs.all=NA, 
                          n.sites.land=NA, n.obs.land=NA, n.sites.land.clear=NA,
                          n.obs.land.clear=NA, n.obs.land.clear.ngb=NA)

col.names <- colnames(fread(files[1], nrows=0))

for (i in 1:n.files){
  n.rows <- nrow(fread(files[i], select = 1L, fill=T))
  rows.per.grab <- 3*10^6 # grab 3 million rows at a time
  row.seq <- seq(0, n.rows, rows.per.grab)
  n.row.seq <- length(row.seq)
  
  for (j in 1:n.row.seq){
    lsat.samples <- fread(files[i],nrows=rows.per.grab, skip=row.seq[j], fill=T)
    colnames(lsat.samples) <- col.names
    lsat.samples[,c("CID"):=NULL] # remove several unnecessary colunms
    lsat.samples <- lsat_general_prep(lsat.samples)
    
    # select from specified years and seasonal window
    lsat.samples <- lsat.samples[year >= 1985][year <= 2016]
    lsat.samples <- lsat.samples[doy >= 152][doy <= 243]
    n.sites.all <- length(unique(lsat.samples$site))
    n.obs.all <- nrow(lsat.samples)
    
    # filter out water
    lsat.samples$water <- apply(X = data.table(pixel.qa=lsat.samples$pixel.qa), MARGIN = 1, FUN = water_flag)
    lsat.samples <- lsat.samples[water == 0]
    lsat.samples <- lsat.samples[jrc.water == 0]
    
    n.sites.land <- length(unique(lsat.samples$site))
    n.obs.land <- nrow(lsat.samples)
    
    # QAQC flags
    lsat.samples <- lsat_qaqc_flags(lsat.samples, filter.water = F) # already filtered water
    n.sites.land.clear <- length(unique(lsat.samples$site))
    n.obs.land.clear <- nrow(lsat.samples)
    
    # neighborhood average
    lsat.samples <- lsat_ngb_mean(lsat.samples)
    n.obs.land.clear.ngb <- nrow(lsat.samples)
    
    # add ndvi and other spectral inidices
    lsat.samples <- lsat_spec_index(lsat.samples, 'ndvi')
    
    # log grab info to data frame
    sample.smry <- rbind(sample.smry,
                         data.frame(region = samples$region[i], sample=samples$sample[i], grab=j, n.sites.all=n.sites.all, n.obs.all=n.obs.all,
                                    n.sites.land=n.sites.land, n.obs.land=n.obs.land, n.sites.land.clear=n.sites.land.clear,
                                    n.obs.land.clear=n.obs.land.clear, n.obs.land.clear.ngb=n.obs.land.clear.ngb))
    
    
    # write out grab of samples and print progress report
    outname <- paste('C:/tmp/lsat_tundra_samples_cleaned/', i, '_', j, '.csv', sep='')
    fwrite(lsat.samples, outname)
    print(paste(j, ' grab of ', n.row.seq, '; ', i, ' of ', n.files, ' files', sep=''))
  }
}

# comibe intermediate files and write to main directory on Google Drive
lsat <- do.call("rbind", lapply(list.files('C:/tmp/lsat_tundra_samples_cleaned/', full.names = T), fread))
fwrite(lsat, 'output/lsat_tundra_samples.csv')

# summarize and write out information on biomwe-wide samples
sample.smry <- na.omit(sample.smry)
sample.smry <- data.table(sample.smry)
biome.sample.smry <- sample.smry[, .(n.sites.all=sum(n.sites.all), n.obs.all=sum(n.obs.all),
                                     n.sites.land=sum(n.sites.land), n.obs.land=sum(n.obs.land),
                                     n.sites.land.clear=sum(n.sites.land.clear), n.obs.land.clear=sum(n.obs.land.clear),
                                     n.obs.land.clear.ngb=sum(n.obs.land.clear.ngb))]
biome.sample.smry

fwrite(biome.sample.smry, 'output/lsat_tundra_biome_sample_summary.csv')
biome.sample.smry <- fread('output/lsat_tundra_biome_sample_summary.csv')
