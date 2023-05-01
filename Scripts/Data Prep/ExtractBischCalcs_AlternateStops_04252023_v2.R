# Code to extract NDVI to points for Adaptive Fidelity
# Author LW
# Created April 11 2023
# edited April 11 2023

#------------------------------------#
# Set up env ####

gc(); rm(list = ls()); gc()

library(sf)
library(tidyverse)
library(lubridate)
library(snow)
library(parallel)
library(sp)
library(beepr)
library(mapview)
library(snowfall)
#bring in functions
source("Z:/MODIS_NDVI/info/ndviEXTRACT.R")
source("Z:/MODIS_NDVI/info/ndviEXTRACT2.R")

#------------------------------------#
# load data and prep ####
setwd("C:/Users/lwilde2/Documents/Chapter2_StopoverFidelity/")
load("Data_out/Data/Stopover/RDH_AlternateStops_20230428.RData")

data.sf <- st_as_sf(blank)
st_crs(data.sf)
#
#data.sf <- st_as_sf(data.sf, coords = c("Long", "Lat"), crs = "+proj=longlat +zone=12 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")


#create x and y
data.sf[,c("x","y")] <- st_coordinates(data.sf)

#make dataframe
data <- data.sf %>% st_drop_geometry() %>% as.data.frame() 
str(data)

data$date <- data$POSIXct
data$xend <- data$x
data$yend <- data$y

data14t22 <- data %>% filter(Year < 2023)

#---------------------------------#
# Extract NDVI ####

folder2 <- "Z:/MODIS_NDVI/NDVI2022"#"//DESKTOP-ug87ror.uwyo.edu//GIS_WyomingPlus//MODIS_NDVI//NDVI2022"
proj <- "+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs+towgs84=0,0,0"

metrics22 <- c("MaxNDVIDay", "MaxIRGday", "SpringStartDay", "SpringEndDay", 
               "IntegratedNDVI", "NDVI_scaled", "IRG_scaled", "MaxBrownDownDate", "SE_springDLCfit","SpringScale", "SpringLength", "sumIRG")

f<-as.data.frame(data14t22); head(f)

for(w in c(1:12)){
#w=1
f$temp <- ndviEXTRACT2(XYdata=data14t22, NDVImetric=metrics22[w], NDVIfolder=folder2,
                      maxcpus= detectCores()-1, xname="x", yname="y", datesname="date",
                      xyCRS=CRS(proj))
names(f)[names(f)=="temp"] <- metrics22[w]
}
head(f)
#beep(sound = 1)
showConnections(); sfStop(); showConnections()


AlternateStops_Bischof <- f

save(AlternateStops_Bischof, file = "Data_out/Data/Stopover/RDH_AlternateStops_Bischof_20230428.RData")



#------------------------------------#
# save ####

#------------------------------------#
# ####