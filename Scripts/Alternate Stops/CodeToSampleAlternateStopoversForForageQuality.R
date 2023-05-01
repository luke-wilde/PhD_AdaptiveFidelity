# Code to bootstrap the stopovers to compare nearby alternatives

#------------------------------------------------------#
#### Set up ####
library(sf)
library(sp)
library(tidyverse)
library(dplyr)
library(lubridate)
library(mapview)
library(MoveTools)
library(terra)
library(patchwork)
library(nngeo)
library(adehabitatHR)
#library(amt)
library(stringr)
library(rgeos)
library(lubridate)
library(parallel)
library(scales)

#------------------------------------------------------#
#### load data ####
load("C:/Users/lwilde2/Documents/Chapter2_StopoverFidelity/Data_out/Data/Stopover/RDH_StopoverBundle_14t22_20230425.RData")

test <- onstop_yr_fin %>% filter(id_yr == "255_2016")

proj <- "+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs+towgs84=0,0,0"
test <- st_transform(test, proj)

#test <- test 

#------------------------------------------------------#
#### set up sequence ####

set.seed(163)


n = 100
xlim = 7000
ylim = 7000

set <- data.frame(xset = runif(n, -xlim, xlim), yset = runif(n, -ylim, ylim))

set.sf <- set %>% st_as_sf(coords = c("xset", "yset"), crs = proj)

#------------------------------------------------------#
#### create available points  ####

blank <- onstop_yr_fin[0,0]
for(i in unique(onstop_yr_fin$id_yr)){
  #i = "255_2016"
  test <- onstop_yr_fin %>% filter(id_yr == i)
  test$flag <- 1
  n = 30
 set.sf1 <- set.sf %>% sample_n(n) %>% st_as_sf()
 
for(j in 1:n){
  #j = 4
  t <- test %>% dplyr::select(-flag)
  x <- st_as_sf(t$geometry + set.sf1$geometry[j], crs = proj)
  st_geometry(x) <- "geometry"
  x$flag <- 0
  x <- st_as_sf(cbind(test %>% st_drop_geometry() %>% dplyr::select(names(test)[c(1:22,36:38,40:43)]),x), crs = proj)
  blank <- rbind(x,blank)
  rm(x)
  }#n
}#id_yr

length(unique(blank$id_yr))
length(unique(onstop_yr_fin$id_yr))

setwd("C:/Users/lwilde2/Documents/Chapter2_StopoverFidelity/")
save(blank, file = "Data_out/Data/Stopover/RDH_AlternateStops_20230428.RData")



shift <- st_as_sf(test$geometry + set.sf$geometry[1], crs = proj)
orig <- test$geometry
shift2 <-  st_as_sf(test$geometry + set.sf$geometry[2], crs = proj)

#shift <- shift %>% st_cast("POINT")


mapview(shift, col.regions = "red") + mapview(orig) + mapview(shift2, col.regions = "black")


