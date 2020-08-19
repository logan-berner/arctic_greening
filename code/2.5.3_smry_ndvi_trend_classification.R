rm(list=ls())
.libPaths(c(.libPaths(), "~/R/", '/home/lb968/R/3.5'))
require(data.table)
require(dplyr)
require(ggplot2)
require(ggpubr)
require(R.utils)
require(zoo) # rollapply
setwd('/projects/arctic/users/lberner/arctic_greening/')
source('/home/lb968/code/boreal_biome_shift/fun_ggplot_stat_boxplot.R')

# LOAD DATA SETS ====================================================================================
accur.dt <- do.call("rbind", lapply(list.files('output/lsat_ndvi_trend_rf/accuracy_mcreps/', full.names = T), fread))
conf.mtx.dt <- do.call("rbind", lapply(list.files('output/lsat_ndvi_trend_rf/conf_mtx_mcreps/', full.names = T), fread))
var.imp.dt <- do.call("rbind", lapply(list.files('output/lsat_ndvi_trend_rf/var_imp_mcreps/', full.names = T), fread))
pd.dt <- do.call("rbind", lapply(list.files('output/lsat_ndvi_trend_rf/partial_depend_mcreps/', full.names = T), fread))

# adust a few variable names
var.imp.dt[var == 'elev', var := 'elev.m']
var.imp.dt[var == 'gst.2003.degC', var := 'gt.2003.degC']
var.imp.dt[var == 'gst.change.degC', var := 'gt.change.degC']
var.imp.dt[var == 'soil.jja.min.change.mm', var := 'sm.change.mm']
var.imp.dt[var == 'soil.jja.min.2000.mm', var := 'sm.2000.mm']
var.imp.dt$var <- toupper(var.imp.dt$var)

pd.dt[var == 'elev', var := 'elev.m']
pd.dt[var == 'gst.2003.degC', var := 'gt.2003.degC']
pd.dt[var == 'gst.change.degC', var := 'gt.change.degC']
pd.dt[var == 'soil.jja.min.change.mm', var := 'sm.change.mm']
pd.dt[var == 'soil.jja.min.2000.mm', var := 'sm.2000.mm']
pd.dt$var <- toupper(pd.dt$var)
pd.dt[class == 'none', class := 'no trend']

conf.mtx.dt[Prediction == 'none', Prediction := 'no trend']
conf.mtx.dt[Reference == 'none', Reference := 'no trend']

accur.dt[class == 'none', class := 'no trend']

# SUMMARIZE ACROSS MONTE CARLO SIMULATIONS ====================================================================================

# classification accuracy
cols <- names(accur.dt)[c(2,3,12,13)]
accur.med.mtx <- as.matrix(accur.dt[, lapply(.SD, function(x){round(median(x),3)}), by = class, .SDcols = cols])
accur.q025.mtx <- as.matrix(accur.dt[, lapply(.SD, function(x){round(quantile(x,0.025),3)}), by = class, .SDcols = cols])
accur.q975.mtx <- as.matrix(accur.dt[, lapply(.SD, function(x){round(quantile(x,0.975),3)}), by = class, .SDcols = cols])
accur.smry.dt <- data.table(matrix(paste0(accur.med.mtx, ' [', accur.q025.mtx, ', ', accur.q975.mtx,']'), ncol = 5))
colnames(accur.smry.dt) <- capitalize(colnames(accur.med.mtx))
accur.smry.dt$Class <- c('browning','none','greening')
accur.smry.dt
fwrite(accur.smry.dt, 'output/lsat_ndvi_trend_rf/lsat_ndvi_trend_rf_classification_accuracy_summary.csv')

# variable importance
var.imp.smry.dt <- var.imp.dt[, .(vi=median(MeanDecreaseAccuracy), vi.q025=quantile(MeanDecreaseAccuracy,0.025), vi.q975=quantile(MeanDecreaseAccuracy,0.975)), by = c('var')]
var.imp.smry.dt <- var.imp.smry.dt[order(-vi)]
top.vars <- var.imp.smry.dt[1:6]
var.imp.smry.dt
fwrite(var.imp.smry.dt, 'output/lsat_ndvi_trend_rf/lsat_ndvi_trend_rf_variable_importance_summary.csv')

# confusion matrix
conf.mtx.dt[, Prediction := factor(Prediction, levels = c('browning','no trend','greening'))]
conf.mtx.dt[, Reference := factor(Reference, levels = c('browning','no trend','greening'))]
conf.mtx.smry.dt <- conf.mtx.dt[, .(N=median(N), N.q025=round(quantile(N,0.025)), N.q975=round(quantile(N,0.975))), by = c('Prediction','Reference')]
conf.mtx.smry.fncy.dt <- conf.mtx.smry.dt[,1:2]
conf.mtx.smry.fncy.dt$N <- paste0(conf.mtx.smry.dt$N, ' [', conf.mtx.smry.dt$N.q025, ', ', conf.mtx.smry.dt$N.q975, ']')
conf.mtx.smry.fncy.dt <- dcast(conf.mtx.smry.fncy.dt, Reference ~ Prediction, value.var = 'N')
conf.mtx.smry.fncy.dt
fwrite(conf.mtx.smry.fncy.dt, 'output/lsat_ndvi_trend_rf/lsat_ndvi_trend_rf_confusion_matrix_summary.csv')

# partial dependency
vars <- unique(pd.dt$var)
var.fac <- vars[c(1:6)]
var.num <- vars[(vars %in% var.fac) == F]
pd.num.dt <- pd.dt[var %in% var.num]
pd.num.dt <- pd.num.dt[, x.value := as.numeric(x.value)]

pd.num.interp.list <- list()
cnt = 1
for (i in unique(pd.num.dt$var)){
  var.dt <- pd.num.dt[var == i]
  x.min <- min(var.dt$x.value)
  x.max <- max(var.dt$x.value)
  x.rng <- x.max - x.min
  x.seq <- seq(x.min, x.max, by = x.rng / 50)
  for (j in unique(var.dt$class)){
    var.class.dt <- var.dt[class == j]
    for (k in unique(var.class.dt$rep)){
      dt <- var.class.dt[rep == k]
      if(sd(dt$x.value)==0){
        next()
      } else {
        prob.fun <- approxfun(dt$x.value, dt$prob)
        ecdf.fun <- approxfun(dt$x.value, dt$ecdf)
        pd.num.interp.list[[cnt]] <- data.table(var = i, class = j, rep = k, x.value = x.seq, prob = prob.fun(x.seq), ecdf = ecdf.fun(x.seq))
        print(cnt)
        cnt = cnt + 1
      }
    }
  }
}
pd.num.dt <- rbindlist(pd.num.interp.list)

pd.num.smry.dt <- pd.num.dt[, .(prob=median(prob, na.rm=T), prob.q025=quantile(prob,0.025, na.rm=T), prob.q975=quantile(prob,0.975, na.rm=T),
                                ecdf = median(ecdf, na.rm=T)), by = c('var','class','x.value')]

# CREATE FIGURES ====================================================================================

varimp.labs <- c('SWI.CHANGE.DEGC' = expression(Delta~'Summer Warmth Index ('*degree*'C)'),
                 'GT.2003.DEGC' = expression('Annual Soil Temperature ('*degree*'C)'),
                 'SWI.2000.DEGC' = expression('Summer Warmth Index ('*degree*'C)'),
                 'ELEV.M' = 'Elevation (m)',
                 'SM.CHANGE.MM' = expression(Delta~'Summer Soil Moisture (mm)'),
                 'GT.CHANGE.DEGC' = expression(Delta~'Annual Soil Temperature ('*degree*'C)'))

pd.labs <- c(paste0("Delta~", "Summer~Warmth~Index~(", "degree*", "C)"),
             paste0("Annual~Soil~Temperature~(", "degree*", "C)"),
             paste0("Summer~Warmth~Index~(", "degree*", "C)"),
             paste0('Elevation~(m)'), 
             paste0("Delta~", "Summer~Soil~Moisture~(mm)"),
             paste0("Delta~", "Annual~Soil~Temperature~(", "degree*", "C)"))

trend.cols <- c('lightsalmon4','gray50','springgreen4')

# variable importance
vap.imp.top.dt <- var.imp.dt[var %in% top.vars$var]
vap.imp.top.dt[, var := factor(var, levels = rev(top.vars$var))]

var.imp.fig <- ggplot() + coord_flip() + theme_bw()
var.imp.fig <- var.imp.fig + stat_boxplot_custom(data = vap.imp.top.dt, mapping = aes(var, MeanDecreaseAccuracy), 
                                                 qs = c(0.025, 0.25,0.5, 0.75, 0.975), geom ='errorbar', width = 0.2)
var.imp.fig <- var.imp.fig + stat_boxplot_custom(vap.imp.top.dt, mapping = aes(var, MeanDecreaseAccuracy), 
                                                 qs = c(0.025, 0.25,0.5, 0.75, 0.975), geom = 'boxplot', width = 0.7, 
                                                 notch = F, outlier.shape = NA, fill='lightgray')
var.imp.fig <- var.imp.fig + xlab('Predictor variable') + ylab('Mean decrease in accuracy')
var.imp.fig <- var.imp.fig + scale_x_discrete(labels = varimp.labs)
var.imp.fig <- var.imp.fig + theme_bw() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 11),
                                                axis.title = element_text(size = 14))

jpeg('figures/Landsat_NDVI_trend_classification_varimp.jpg', 6, 5, res = 400, units = 'in')
var.imp.fig
dev.off()


# partial dependency of most important vars
pd.smry.top.dt <- pd.num.smry.dt[var %in% top.vars$var]
pd.smry.top.dt[, var := factor(var, levels = top.vars$var, labels = pd.labs)]
pd.smry.top.dt[, class := factor(class, levels = c('browning','no trend','greening'))]
# pd.smry.top.dt <- pd.smry.top.dt[class != 'none']
# pd.smry.top.dt <- pd.smry.top.dt[ecdf > 0.05 & ecdf < 0.95]

pd.fig <- ggplot(pd.smry.top.dt, aes(x.value, prob, group = class)) + geom_line(aes(color=class))
pd.fig <- pd.fig + geom_ribbon(aes(x=x.value, ymin=prob.q025, ymax=prob.q975, fill=class), alpha = 0.25)
pd.fig <- pd.fig + scale_fill_manual(values = trend.cols) + scale_color_manual(values = trend.cols)
pd.fig <- pd.fig + facet_wrap(var ~ ., scales = 'free', labeller = label_parsed)
pd.fig <- pd.fig + xlab('Value of predictor variable') + ylab('Classification probability') 
pd.fig <- pd.fig + labs(color=expression('NDVI'[max]~'trend class: '), fill = expression('NDVI'[max]~'trend class: '))
pd.fig <- pd.fig + theme_bw() + theme(legend.position = 'top', legend.text=element_text(size=14), legend.title=element_text(size=14),
                                      strip.text = element_text(size = 9), axis.text = element_text(size = 12), axis.title = element_text(size = 14))

jpeg('figures/Landsat_NDVI_trend_classification_partial_depend.jpg', 10, 8, res = 400, units = 'in')
pd.fig
dev.off()


# combined var importance and partial dependency
combo.fig <- ggarrange(var.imp.fig, pd.fig, labels = c("(a)", "(b)"), label.x = -0.01, label.y = 1.04, ncol = 1, nrow = 2, heights=c(1,2))

pdf('figures/Landsat_NDVI_trend_classification_varImp_parDepend.pdf',7, 7)
# jpeg('figures/Landsat_NDVI_trend_classification_varImp_parDepend.jpg', 11, 4, res = 500, units = 'in')
combo.fig
dev.off()