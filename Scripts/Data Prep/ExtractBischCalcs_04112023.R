# Code to extract NDVI to points for Adaptive Fidelity
# Author LW
# Created April 11 2023
# edited April 11 2023

#------------------------------------#
# Set up env ####

library(sf)
library(tidyverse)
library(lubridate)
library(snow)
library(parallel)
library(sp)
library(beepr)

#bring in functions
source("Z:/MODIS_NDVI/info/ndviEXTRACT.R")
source("Z:/MODIS_NDVI/info/ndviEXTRACT2.R")

#------------------------------------#
# load data and prep ####

load("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Raw Data/GPS/Spring Migrations/AllSpring_2014to2022.RData")

#clean points too fast
mig.SpringMigration.sf <- mig.SpringMigration.sf %>% arrange(id_yr, POSIXct) %>% group_by(id_yr) %>% mutate(dist = st_distance(geometry, geometry[row_number() + 1], by_element = T), diff = POSIXct[row_number() + 1] - POSIXct, speed = dist / (as.numeric(diff, units="hours") *3600)) %>% ungroup()
mig.SpringMigration.sf <- mig.SpringMigration.sf %>% filter(as.numeric(speed) < 3)
mig.SpringMigration.sf <- mig.SpringMigration.sf %>% dplyr::select(-c(dist, diff, speed))

#ensure that x and y are there
mig.SpringMigration.sf[,c("x","y")] <- st_coordinates(mig.SpringMigration.sf)


#seperate periods
pre <- mig.SpringMigration.sf %>% mutate(Year = year(POSIXct), POSIXct = as.POSIXct(POSIXct), Date = dmy(POSIXct)) %>% dplyr::filter(Year < 2022) %>% st_drop_geometry() %>% as.data.frame()
post <- mig.SpringMigration.sf %>% mutate(Year = year(POSIXct), POSIXct = as.POSIXct(POSIXct)) %>% dplyr::filter(Year >= 2022) %>% st_drop_geometry() %>% as.data.frame()

str(pre)

#------------------------------------#
# extract function GIS_WyomingPlus ####


metrics <- c("maxIRGdate","fittedVals","IRGVals","sumIRG","springStart", "springEnd", "csumNDVImax") #The NDVImetrics can be maxNDVIdate, maxIRGdate, springStart, springEnd, csumNDVImax, fittedVals, IRGVals, maxBrownDowndate, xmidS_SD,SpringScale, springLength, or sumIRG

#first, years before 2022
f<-as.data.frame(pre); head(f)

#for(w in 1:length(metrics)){
  w=1
  f$temp <- ndviEXTRACT(XYdata=pre[1:500,], NDVImetric=metrics[w], NDVIfolder="Z:/MODIS_NDVI/Bischof_calculations",
                        maxcpus= detectCores()-1, xname="x", yname="y", datesname="POSIXct", scaleIRG=TRUE,
                        xyCRS=CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  names(f)[names(f)=="temp"] <- metrics[w]
}
head(f)
beep(sound = 1)

#throwing error - error in evaluating the argument 'x' in selecting a method for function 'merge': undefined columns selected





#then, years after 2022
ndviEXTRACT2(XYdata=post, NDVImetric="NDVI_scaled", NDVIfolder="Z:/MODIS_NDVI/NDVI2022",
             maxcpus= detectCores()-1, xname="x", yname="y", datesname="POSIXct", scaleIRG=TRUE,
             xyCRS=CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))



#------------------------------------#
# save ####

#------------------------------------#
# ####