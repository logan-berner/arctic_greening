# ABOUT THIS SCRIPT  ==============================================================================================================
# This R script summarizes Landsat cross-sensor calibration random forest models
# AUTHOR: LOGAN BERNER, NAU
# DATE: 2020-02-25

# SET UP WORKSPACE ==============================================================================================================
rm(list=ls())
require(data.table)
require(ggplot2)

setwd('/projects/above_gedi/lberner/arctic_greening/')

# READ IN FILES
files <- list.files('output/xcal_ndvi/', full.names = T)
files <- files[grep('eval.csv', files)]

xcal.dt <- do.call("rbind", lapply(files, fread))
xcal.dt$rep <- sort(rep(1:1000, 2))

# SUMMARIZE EVALUATION COEFFICIENTS ACROSS MONTE CARLO ITERATIONS
xcal.smry.dt <- xcal.dt[, .(rf.r2=paste0(round(median(rf.r2),3), ' [', round(quantile(rf.r2,0.025),3),',',round(quantile(rf.r2,0.975),3),']'),
                            rf.rmse=paste0(round(median(rf.rmse),3), ' [', round(quantile(rf.rmse,0.025),3),',',round(quantile(rf.rmse,0.975),3),']'),
                            rf.n=paste0(round(median(rf.n)), ' [', round(quantile(rf.n,0.025)),',',round(quantile(rf.n,0.975)),']'),
                            xval.r2=paste0(round(median(xval.r2),3), ' [', round(quantile(xval.r2,0.025),3),',',round(quantile(xval.r2,0.975),3),']'),
                            xval.rmse=paste0(round(median(xval.rmse),3), ' [', round(quantile(xval.rmse,0.025),3),',',round(quantile(xval.rmse,0.975),3),']'),
                            xval.bias=paste0(round(median(xval.bias),3), ' [', round(quantile(xval.bias,0.025),3),',',round(quantile(xval.bias,0.975),3),']'),
                            xval.n=paste0(round(median(xval.n)), ' [', round(quantile(xval.n,0.025)),',',round(quantile(xval.n,0.975)),']')),
                        by = 'sat']

xcal.smry.dt

# WRITE OUT SUMMARY TABLE
fwrite(xcal.smry.dt, 'output/lsat_ndvi_xcal_smry.csv')
