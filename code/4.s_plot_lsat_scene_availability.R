# This R script fits phenological curves to to site-level NDVI time series and then uses the curves to estimates annual peak summer NDVI 
# Date: 2018-01-09
rm(list=ls())
require(tidyr)
require(dplyr)
require(lattice)
require(zoo)
require(data.table)
setwd('/projects/arctic/users/lberner/arctic_greening/')
source('/home/lb968/code/arctic_greening/0.3_fun_shade_CIs.r')


# READ IN NDVI TIME SERIES FROM GEE AND CLEAN UP / REORDER COLUMNS #------------------------------------
vi.ts <- fread('data/lsat_samples/lsat_tundra_samples_n50000.csv')
head(vi.ts)
dim(vi.ts)

#--------------------------------------------------------------------------------------------------------
# SUMMARIZE NUMBER OF SCENES AVAILABLE EACH YEAR
#--------------------------------------------------------------------------------------------------------
obs.site.yr.smry <- vi.ts %>% group_by(site, year) %>% summarise(n.obs=n()) 
site.yr.combos <- obs.site.yr.smry %>% group_by() %>% expand(site, 1984:2016) %>% select(site, year = `1984:2016`) # compute all site x year combinations to fill in missing data (e.g., years w/ no observations)
obs.site.yr.smry <- obs.site.yr.smry %>% right_join(site.yr.combos)
obs.site.yr.smry$n.obs[is.na(obs.site.yr.smry$n.obs)] <- 0

# Maximum number of scenes for a site x year
max(obs.site.yr.smry$n.obs)


# Number of scenes per YEAR 
obs.yr.smry <- obs.site.yr.smry %>% group_by(year) %>% summarise(n.obs.med=quantile(n.obs, 0.5), n.obs.q25=quantile(n.obs, 0.25), n.obs.q75=quantile(n.obs, 0.75),
                                                        n.obs.q05=quantile(n.obs, 0.05), n.obs.q95=quantile(n.obs, 0.95))

obs.yr.smry %>% filter(year == 1985 | year == 1995 | year == 2005 | year == 2015)

# Number of scenes per SITE
obs.site.smry <- obs.site.yr.smry %>% group_by(site) %>% summarise(n.obs.tot = sum(n.obs), n.obs.med=quantile(n.obs, 0.5), n.obs.q25=quantile(n.obs, 0.25), n.obs.q75=quantile(n.obs, 0.75),
                                                                 n.obs.q05=quantile(n.obs, 0.05), n.obs.q95=quantile(n.obs, 0.95))
fivenum(obs.site.smry$n.obs.tot) # sumamrize the number of obs among sites  

# Number of scenes per DOY
obs.doy.smry <- vi.ts %>% group_by() %>% mutate(n.obs.tot = n()) %>% group_by(doy) %>% summarise(n.obs = n(), n.obs.tot = first(n.obs.tot)) %>%
  group_by() %>% mutate(n.obs.pcnt = n.obs/n.obs.tot*100)

# fraction and cumulative fraction of site x years with given number of scenes 
obs.cumfrac.smry <- obs.site.yr.smry %>% group_by(n.obs) %>% summarise(site.yr.cnt = n()) %>%
  mutate(total.obs = sum(site.yr.cnt), frac = site.yr.cnt / total.obs, cum.frac = cumsum(frac))

sum(obs.cumfrac.smry$frac[2:6]) # frac of site x years with 1-5 scenes
sum(obs.cumfrac.smry$frac[7:length(obs.cumfrac.smry$frac)]) # frac of site x years with > 5 scenes

require(zyp)
zyp.yuepilon(obs.yr.smry$n.obs.med, obs.yr.smry$year)

cor.test(obs.yr.smry$n.obs.med, obs.yr.smry$year, method = 'spearman')

fivenum(obs.smry$n.obs)# compute percentiles


# length of snow free season 
snow.free.season <- vi.ts %>% group_by(site) %>% summarize(doy.first=min(doy), doy.last=max(doy)) %>% mutate(doy.rng = doy.last-doy.first)

hist(snow.free.season$doy.rng, xlab = "Length of snow-free season (days)", ylab='Number of sites')
mean(snow.free.season$doy.rng) # mean number of snow-free days 
sd(snow.free.season$doy.rng)

mean(snow.free.season$doy.first) # mean first snow-free day
mean(snow.free.season$doy.last) # mean last snow-free day

#--------------------------------------------------------------------------------------------------------
# PLOT NUMBER OF SCENES AVAILABLE EACH YEAR
#--------------------------------------------------------------------------------------------------------
jpeg('figures/Landsat_Nscenes.jpg', 7,3.5,units = 'in', res=400)
par.op = par(mfrow = c(1,2))
par.top = par(mar = c(4,4,1,1))
plot(n.obs.med~year, obs.yr.smry, type='l', lwd=2, ylim=c(0,20), xlim=c(1984,2017), yaxt='n',  
     xlab="Year", ylab='Number of Landsat scenes')
axis(2, at = seq(0,20,5), las=1)
shade.confidence.interval(obs.yr.smry$year, obs.yr.smry$n.obs.q25, obs.yr.smry$n.obs.q75, band.col = 'gray30')
shade.confidence.interval(obs.yr.smry$year, obs.yr.smry$n.obs.q05, obs.yr.smry$n.obs.q95, band.col = 'gray60')
lines(obs.yr.smry$year, obs.yr.smry$n.obs.med, lwd=2)
abline(v = c(1984, 1999,2013), lty=2, lwd=0.75, col='gray50')
text(c(1984,1999,2013)+1, 17, c('L5 launched','L7 launced','L8 launched'), srt=90, cex=0.8)
text(2016.5, 20, '(a)', font=2)

par(par.top)

par.bot = par(mar = c(4,4,1,1))
obs.site.yr.smry.lte25 <- subset(obs.site.yr.smry, n.obs <= 25)
hist(obs.site.yr.smry.lte25$n.obs, breaks = 0:26, xlim=c(0,26), ylim=c(0,0.26), probability = T, yaxt='n', xaxt = 'n', col = 'gray80',
     ylab='Proportion of site x year combinations', xlab='Number of Landsat scenes ', main='', right = F)
axis(2, at = seq(0,0.25,0.05), las=1)
axis(1, at = seq(0.5,25.5,5), labels = seq(0,25,5), las=1)
text(25, 0.26, '(b)', font=2)
box()
par(par.bot)
par(par.op)
dev.off()


quantile(obs.site.yr.smry$n.obs, 0.999)


#--------------------------------------------------------------------------------------------------------
# PLOT NUMBER OF SCENES EACH DOY
#--------------------------------------------------------------------------------------------------------

plot(n.obs.pcnt~doy, obs.doy.smry, type='l', xlab = 'Day of Year', ylab='Percent of Landsat scenes')
