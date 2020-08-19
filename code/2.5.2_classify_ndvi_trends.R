rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(data.table)
require(dplyr)
require(randomForest)
require(caret) # confusionMatrix() and train()
require(pdp) # partial dependency plots
require(R.utils)

args <- commandArgs(TRUE)
i = as.numeric(args[1])
# i = 120

setwd('/projects/arctic/users/lberner/arctic_greening/')

# LOAD DATA SETS AND CREATE OUTPUT DIRS ====================================================================================
ndvi.trnds.dt <- fread(list.files('output/lsat_site_trends/mc_reps/',full.names = T)[i])
site.dt <- fread('output/lsat_site_trend_modeling/tundra_site_biogeoclim.csv', fill=T)

# create output dirs
mkdirs('output/lsat_ndvi_trend_rf/')
mkdirs('output/lsat_ndvi_trend_rf/rf_models_mcreps')
mkdirs('output/lsat_ndvi_trend_rf/accuracy_mcreps')
mkdirs('output/lsat_ndvi_trend_rf/conf_mtx_mcreps')
mkdirs('output/lsat_ndvi_trend_rf/var_imp_mcreps')
mkdirs('output/lsat_ndvi_trend_rf/partial_depend_mcreps')

# reformat i for sorting
if (i < 10){
  i <- paste0('000',i)
} else if (i < 100){
  i <- paste0('00',i)
} else if (i < 1000){
  i <- paste0('0',i)
}

# COMBINE LANDSAT TREND DATA WITH SITES CONDITIONS ==================================================

# focus on 2000-2016
ndvi.trnds.dt <- ndvi.trnds.dt[period == '2000-2016']

# reclassify NDVI trend categories
ndvi.trnds.dt[trend.cat == 'greening.sig.p5' | trend.cat == 'greening.sig.p10', trend.cat := 'greening']
ndvi.trnds.dt[trend.cat == 'browning.sig.p5' | trend.cat == 'browning.sig.p10', trend.cat := 'browning']
ndvi.trnds.dt[trend.cat == 'insig', trend.cat := 'none']
ndvi.trnds.dt[, trend.cat := factor(trend.cat, levels = c('browning','none','greening'))]
trend.classes <- c('browning','none','greening')

# add specific columns to site data table
cols <- c('site','trend.cat')
ndvi.trnds.dt <- ndvi.trnds.dt[, ..cols]
site.dt <- ndvi.trnds.dt[site.dt, on = 'site']

# drop any row without a ndvi trend
site.dt <- site.dt[is.na(trend.cat) ==F]

# PREPARE DATA FOR MODELING  ==================================================

# check data
summary(site.dt)

# drop a few cols 
site.dt <- site.dt[, c('latitude','longitude','site') := NULL]

# deal with some black values that should be NA
site.dt[tk_hill == '', tk_hill := NA]
site.dt[tk_wetland == '', tk_wetland := NA]
site.dt[tk_lake == '', tk_lake := NA]

# drop nas
site.dt <- na.omit(site.dt)

# check data again
summary(site.dt)

# balance sample size across trend classes
n.brown <- as.numeric(site.dt[trend.cat == 'browning', .(n.sites = .N)])
site.dt <- site.dt[, .SD[sample(.N, n.brown, replace = F)], by = c('trend.cat')]

# turn burn year into a factor
site.dt[, burnyr := as.factor(burnyr)]

# set characters to factors
site.dt <- site.dt %>% mutate_if(is.character,as.factor) %>% as.data.table()

# screen highly correlated variables if present
print('screening highly correlated values...')
site.fac.dt <- site.dt[, sapply(site.dt, is.numeric)==F, with=F] 
site.num.dt <- site.dt[, sapply(site.dt, is.numeric), with=F]
cors <- cor(site.num.dt, use = 'pair')
high.cors <- findCorrelation(cors, cutoff = 0.75)
site.num.dt <- site.num.dt[, (high.cors) := NULL]
site.dt <- cbind(site.fac.dt, site.num.dt)

# check data again
summary(site.dt)

# # subset data for testing
# site.dt <- site.dt[sample(1:nrow(site.dt), 500)]

# split into training and eval data
train.indx <- sample(1:nrow(site.dt), size = nrow(site.dt)*0.666)
site.train.dt <- site.dt[train.indx]
site.eval.dt <- site.dt[-train.indx]

# PREDICT NDVI TREND CATEGORY  ===========================================================

# fit random forest model
print('fitting models...')
control <- trainControl(method="oob", number=20, search="random")
rf.trained <- train(trend.cat~., data=site.train.dt, method="rf", metric="Accuracy", trControl=control, ntree = 500, importance = T)

# derive confusion matrix
print('starting model evaluation...')
site.eval.dt$predicted <- predict.train(rf.trained, newdata = site.eval.dt[,-1])
conf.mtx <- confusionMatrix(reference = site.eval.dt$trend.cat, data = site.eval.dt$predicted, mode='sens_spec')
conf.mtx.dt <- data.table(conf.mtx$table)
conf.mtx.dt[, rep := i]
accur.dt <- data.table(class = trend.classes, as.data.table(conf.mtx$byClass))
accur.dt[, OverallAccuracy := conf.mtx$overall[1]]
accur.dt[, rep := i]

# variable importance
rf.best.fit <- rf.trained$finalModel
var.imp.dt <- cbind(data.table(var = rownames(rf.best.fit$importance)), data.table(rf.best.fit$importance))[order(MeanDecreaseAccuracy, decreasing = T)]
var.imp.dt[, rep := i]

# partial dependency
print('computing partial dependencies...')

pd.list <- list()
cnt = 1
vars <- names(site.dt)[-1]
for (j in vars){
  for (k in trend.classes){
    pd <- data.table(partial(object = rf.trained, pred.var = j, which.class = k, train = site.dt, smooth = T, type = 'classification', prob = T, trim.outliers = T))
    names(pd) <- c('x.value','prob')
    
    # compute ecdf to help determine where caution needed during interpretation
    fun.ecdf <- ecdf(site.dt[[j]]) # returns NAs for categorical variables
    pd[, ecdf := fun.ecdf(pd$x.value)]
    
    # store output
    pd <- cbind(pd, data.table(var = j, class = k, rep = i))
    pd.list[[cnt]] <- pd
    cnt = cnt + 1
    print(paste(j,k, sep=' '))
  }
}

pd.dt <- rbindlist(pd.list)

# write output to disk
print('writting output...')
save(rf.best.fit, file = paste0('output/lsat_ndvi_trend_rf/rf_models_mcreps/ndvi_trend_rf_classifier_mcrep_',i,'.RDATA'))
fwrite(conf.mtx.dt, paste0('output/lsat_ndvi_trend_rf/conf_mtx_mcreps/ndvi_trend_rf_conf_mtx_',i,'.csv'))
fwrite(accur.dt, paste0('output/lsat_ndvi_trend_rf/accuracy_mcreps/ndvi_trend_rf_accuracy_',i,'.csv'))
fwrite(var.imp.dt, paste0('output/lsat_ndvi_trend_rf/var_imp_mcreps/ndvi_trend_rf_var_imp_',i,'.csv'))
fwrite(pd.dt, paste0('output/lsat_ndvi_trend_rf/partial_depend_mcreps/ndvi_trend_rf_partial_dependency_',i,'.csv'))

print("All done!")

# HOLDING ----------------------------------------------------------------------------------------------------------

# site.dt[burnyr == 0, burned := 'no']
# site.dt[burnyr != 0, burned := 'yes']
# burn.smry <- site.dt[, .(n.sites = .N), by = burned]
# burn.smry <- burn.smry[, pcnt := n.sites / sum(n.sites)]
# 
# burn.smry <- site.dt[, .(n.sites = .N), by = c('trend.cat', 'burned')]
# burn.smry <- burn.smry[, pcnt := n.sites / sum(n.sites), by = trend.cat]
