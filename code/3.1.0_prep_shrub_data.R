# This R script complise shrub ring-width data from several sources into one data base
rm(list=ls())
require(dplyr)
require(data.table)
require(reshape2)
require(dplR)
setwd('/projects/above_gedi/lberner/arctic_greening/data/')

# READ IN, CLEAN UP, AND COMPILE RING-WIDTH MEASUREMENTS ..... 

# BEN GAGLIOTI'S SHRUB DATA  ====================================================================
ben.burn.wide <- read.rwl('field_data/tundra_shrub_growth/gaglioti_ak_shrubs/Alder_B_Final')
ben.burn.wide$year <- rownames(ben.burn.wide)
ben.burn.wide <- as.data.frame(ben.burn.wide)
ben.burn <- melt(ben.burn.wide, id.vars=c('year'), variable.name='series',value.name='ring.width.mm')
ben.burn <- na.omit(ben.burn)
ben.burn$series <- substr(ben.burn$series, 1,6)
ben.burn$series <- gsub(pattern = " ", replacement = 'B', ben.burn$series)
ben.burn$series <- gsub(pattern = "_", replacement = '', ben.burn$series)
ben.burn$site <- "Noatak_RB"
ben.burn$plant <- paste('RB', substr(ben.burn$series, 3,4), sep='_')
ben.burn$radii <- substr(ben.burn$series, 5,5)
ben.burn$radii[ben.burn$radii != "B"] <- 'A'
ben.burn$species <- 'Alnus fruticosa'
ben.burn$genus <- 'Alnus'
head(ben.burn)

ben.oldburn.wide <- read.rwl('field_data/tundra_shrub_growth/gaglioti_ak_shrubs/Alder_OB_Final')
ben.oldburn.wide$year <- rownames(ben.oldburn.wide)
ben.oldburn.wide <- as.data.frame(ben.oldburn.wide)
ben.oldburn <- melt(ben.oldburn.wide, id.vars=c('year'), variable.name='series',value.name='ring.width.mm')
ben.oldburn <- na.omit(ben.oldburn)
ben.oldburn$series <- substr(ben.oldburn$series, 1,6)
ben.oldburn$series <- gsub(pattern = " ", replacement = 'B', ben.oldburn$series)
ben.oldburn$site <- 'Noatak_OB'
ben.oldburn$plant <- paste('OB', substr(ben.oldburn$series, 4,5), sep='_')
ben.oldburn$radii <- substr(ben.oldburn$series, 6,7)
ben.oldburn$radii[ben.oldburn$radii != "B"] <- 'A'
ben.oldburn$species <- 'Alnus fruticosa'
ben.oldburn$genus <- 'Alnus'
head(ben.oldburn)

ben.noburn.wide <- read.rwl('field_data/tundra_shrub_growth/gaglioti_ak_shrubs/Alder_NB_Final')
ben.noburn.wide$year <- rownames(ben.noburn.wide)
ben.noburn.wide <- as.data.frame(ben.noburn.wide)
ben.noburn <- melt(ben.noburn.wide, id.vars=c('year'), variable.name='series',value.name='ring.width.mm')
ben.noburn <- na.omit(ben.noburn)
ben.noburn$series <- substr(ben.noburn$series, 1,6)
ben.noburn$series <- gsub(pattern = " ", replacement = 'B', ben.noburn$series)
ben.noburn$site <- 'Noatak_NB'
ben.noburn$plant <- paste('NB', substr(ben.noburn$series, 4,5), sep='_')
ben.noburn$radii <- substr(ben.noburn$series, 6,7)
ben.noburn$radii[ben.noburn$radii != "B"] <- 'A'
ben.noburn$species <- 'Alnus fruticosa'
ben.noburn$genus <- 'Alnus'

ben.nslope <- rbind(ben.noburn, ben.oldburn, ben.burn)
ben.nslope$contributor <- "B. Gaglioti"
ben.nslope$cluster <- 'Noatak'
ben.nslope$shrub.hub <- 'no'
ben.nslope <- ben.nslope %>% select(site,cluster,plant,radii,year,ring.width.mm,genus, species, contributor, shrub.hub)
head(ben.nslope)


# LAIA ANDREU-HAYLES SHRUB DATA FROM NORTH SLOPE OF AK  ====================================================================
laia.nslope.wide <- read.rwl('field_data/tundra_shrub_growth/hayles_ak_shrubs/All_Shrubs.txt')
laia.nslope.wide$year <- rownames(laia.nslope.wide)
laia.nslope.wide <- as.data.frame(laia.nslope.wide)
laia.nslope <- melt(laia.nslope.wide, id.vars=c('year'), variable.name='series',value.name='ring.width.mm')
laia.nslope <- na.omit(laia.nslope)
laia.nslope$site <- substr(laia.nslope$series, 1,2)
laia.nslope <- subset(laia.nslope, site == 'AL' | site == 'SG' | site == 'IC')
laia.nslope$plant <- paste(laia.nslope$site, substr(laia.nslope$series, 3,5), sep='_')
laia.nslope$radii <-substr(laia.nslope$series, 6,6)
laia.nslope$species <- 'Salix glauca' # all are willow except for AL
laia.nslope$species[laia.nslope$site == "AL"] <- 'Alnus fruticosa'
laia.nslope$genus <- 'Salix' # all are willow except for AL
laia.nslope$genus[laia.nslope$site == "AL"] <- 'Alnus'
laia.nslope$cluster[laia.nslope$site == "AL"] <- "Sagwon_N"
laia.nslope$cluster[laia.nslope$site == "IC"] <- 'Sagwon_N'
laia.nslope$cluster[laia.nslope$site == "SG"] <- 'Sagwon_S'
laia.nslope$site[laia.nslope$site == "AL"] <- "Sagwon_AL"
laia.nslope$site[laia.nslope$site == "IC"] <- 'Sagwon_IC'
laia.nslope$site[laia.nslope$site == "SG"] <- 'Sagwon_SG'
laia.nslope$contributor <- "L. Andrue-Hayles"
laia.nslope$shrub.hub <- 'no'
laia.nslope <- laia.nslope %>% select(site,cluster,plant,radii,year,ring.width.mm,genus, species, contributor, shrub.hub)
head(laia.nslope)


# MARC AND BRUCE'S SHRUB DATA FROM EURASIA  ====================================================================
marc.files <- list.files('field_data/tundra_shrub_growth/forbes_eurasia_shrubs/csv/', full.names = T)
marc.files.short <- list.files('field_data/tundra_shrub_growth/forbes_eurasia_shrubs/csv/', full.names = F)
marc.eurasia <- data.frame(year=NA, series=NA, ring.width.mm=NA, site=NA, plant=NA, radii=NA, species=NA, genus=NA)
for (i in 1:length(marc.files)){
  site.name <- as.character(as.data.frame(strsplit(marc.files.short[i], '_'))[3,1])
  site.sp <- paste(as.data.frame(strsplit(marc.files.short[i], '_'))[1,1], as.data.frame(strsplit(marc.files.short[i], '_'))[2,1], sep=' ')
  site.dat <- read.table(marc.files[i], sep = ',', header = T)  
  colnames(site.dat) <- tolower(colnames(site.dat))
  site.dat <- melt(site.dat, id.vars = 'year', variable.name = 'series', value.name='ring.width.mm')
  site.dat <- na.omit(site.dat)
  site.dat$site <- gsub(" ", "", site.name)
  site.dat$plant <- site.dat$series
  site.dat$radii <- 'A'
  site.dat$species <- site.sp
  site.dat$genus <- if(substr(site.sp,1,1) == "A"){"Alnus"} else{'Salix'}
  if (site.name == 'Cherskii' | site.name == 'Enontekio'){
    site.dat$ring.width.mm <- site.dat$ring.width.mm / 100
  }
  marc.eurasia <- rbind(marc.eurasia, site.dat)
  print(site.name)
  print(unique(site.dat$series))
}
marc.eurasia <- na.omit(marc.eurasia)
marc.eurasia$contributor <- 'B. Forbes and M. Macias-Fauria'
marc.eurasia$cluster <- marc.eurasia$site
marc.eurasia$cluster[marc.eurasia$site == 'MordyYaha'] <- 'Mordy Yaha'
marc.eurasia$shrub.hub <- 'no'
marc.eurasia <- marc.eurasia %>% select(site,cluster,plant,radii,year,ring.width.mm,genus, species, contributor, shrub.hub)
head(marc.eurasia)


# PADDY AND ERIC'S SHRUB DATA FROM GREENLAND  ====================================================================
paddy.betula <- read.csv('field_data/tundra_shrub_growth/sullivan_greenland_shrubs/betula_nana/data/betula.ring.widths.long.csv')
paddy.salix <- read.csv('field_data/tundra_shrub_growth/sullivan_greenland_shrubs/salix_glauca/data/salix.ring.widths.long.csv')

paddy.salix <- rename(paddy.salix, plant=disc, ring.width.mm = ring.width_mm)
paddy.salix$site <- substr(paddy.salix$plant, 0,2)
paddy.salix$site[paddy.salix$site != 'L1' & paddy.salix$site != 'L2'] <- substr(paddy.salix$site[paddy.salix$site != 'L1' & paddy.salix$site != 'L2'], 0,1)
paddy.salix$site <- recode(paddy.salix$site, 'H'='Haredalen', 'L1'='LongLakeI', 'L2'='LongLakeII', 'S'='Sandflugtsdalen', 'T'='West')
paddy.salix$cluster <- "Kanger"
paddy.salix$radii <- 'A'
paddy.salix$genus <- 'Salix'
paddy.salix$species <- 'Salix glauca'
paddy.salix$contributor <- 'P. Sullivan and E. Post'
paddy.salix$shrub.hub <- 'no'
paddy.salix <- paddy.salix %>% select(site,cluster,plant,radii,year,ring.width.mm,genus, species, contributor, shrub.hub)  
head(paddy.salix)

paddy.betula <- rename(paddy.betula, plant=disc, ring.width.mm = ring.width_mm)
paddy.betula$site <- substr(paddy.betula$plant, 0,2)
paddy.betula$site[paddy.betula$site != 'L1' & paddy.betula$site != 'L2'] <- substr(paddy.betula$site[paddy.betula$site != 'L1' & paddy.betula$site != 'L2'], 0,1)
paddy.betula$site <- recode(paddy.betula$site, 'H'='Haredalen', 'K'='Kodalen', 'L1'='LongLakeI', 'S'='Sandflugtsdalen', 'T'='West')
paddy.betula$cluster <- "Kanger"
paddy.betula$radii <- 'A'
paddy.betula$genus <- 'Betula'
paddy.betula$species <- 'Betula nana'
paddy.betula$contributor <- 'P. Sullivan and E. Post'
paddy.betula$shrub.hub <- 'no'
paddy.betula <- paddy.betula %>% select(site,cluster,plant,radii,year,ring.width.mm,genus, species, contributor, shrub.hub)  
head(paddy.betula)

paddy.greenland <- rbind(paddy.salix, paddy.betula)


# ISLA'S SHRUB HUB SYNTHESIS DATA SET  ====================================================================
shrubhub <- fread('field_data/tundra_shrub_growth/shrubhub_synthesis_shrubs/shrubhub_tundra_shrub_growth.csv')

shrubhub <- shrubhub[sampling_method != 'WMS'] # don't use wintermarked septa
shrubhub <- shrubhub[cross_referenced == 'SI' | cross_referenced == 'TC' | cross_referenced == 'Si' | cross_referenced == 'TC, SI'] # only take series that were cross-dated to tree or other shrubs
shrubhub <- shrubhub[pi_full != "Bruce C. Forbes/Marc Macias-Fauria"] # don't take mmnts from Bruce and Marc since we have updated data from them

# simplify several site names
shrubhub[site_full == 'Nowell Lake, NWT', site_full := "Nowell Lake"]
shrubhub[site_full == 'Arsuk Fjord, S. Greenland', site_full := "Arsuk Fjord"]
shrubhub[site_full == 'Sweden,Abisko,Slottatjockasouthslope', site_full := "Abisko"]
shrubhub[site_full == 'Zackenberg Fell-field', site_full := "Zackenberg"]
shrubhub[site == 's24', situ_full := 'Zackenberg']
setnames(shrubhub, old = 'site_full', new = 'cluster')
shrubhub$genus <- recode(shrubhub$genus, 'A'='Alnus', 'P'="Pinus", 'S'='Salix', 'J'="Juniper", 'B'='Betula')
shrubhub$shrub.hub <- 'yes'
shrubhub <- shrubhub %>% dplyr::select(site,cluster,plant=shrub_id,radii,year,ring.width.mm=rw,genus, species, contributor=pi_full, shrub.hub)  
head(shrubhub)


# COMPILE DATA SETS AND WRITE OUT====================================================================
shrub.rw <- rbind(ben.nslope, laia.nslope, marc.eurasia, paddy.greenland, shrubhub)
shrub.rw.check <- shrub.rw %>% group_by(site) %>% summarise(rw.avg = mean(ring.width.mm)) # check that order of magnitude is correct... what order is correct??
fivenum(shrub.rw.check$rw.avg)
length(unique(shrub.rw$plant))

fwrite(shrub.rw, 'field_data/tundra_shrub_growth/tundra_shrub_ring_width.csv')