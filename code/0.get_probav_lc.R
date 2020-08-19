# THIS R SCRIPT DOWNLOADS PROBA-V LANDCOVER DATA FROM THE INTERNET
# AUTHOR: LOGAN BERNER, NAU
# DATE: 2020-2-17

rm(list=ls())
require(R.utils)
mkdirs('probaV_lc')
setwd('/projects/above_gedi/geodata/landcover/probaV_lc/')

urls <- read.csv('probav_urls.csv')

n.files <- nrow(urls)

for (i in 1:n.files){
  dest <- paste(getwd(), urls$file[i], sep='/')
  url <- as.vector(urls$url[i])
  download.file(url = url, destfile = dest)
  print(i/n.files)
}
