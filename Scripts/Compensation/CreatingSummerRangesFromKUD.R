#Code for creating 95% KUDs of summer range points
#Written by Jerod Merkle and modified by Anna Ortega, Wyoming Cooperative Fish and Wildlife Research Unit
#December 28, 2021

#Load required libraries
library(adehabitatHR)
library(BBMM)
library(move)
library(snowfall)
library(stringr)

#Clean up working environment
rm(list=ls())

#Import GPS collar data
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange")

load("AllRD2H_2011_thru_2020_LongDistance_SummerRangePoints_CleanForAnalyses.RData")

#Make database spatially aware
proj <- "+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs+towgs84=0,0,0"   #name a new proj
coordinates(data) <- c("x","y")
proj4string(data) <- CRS(proj)

#-------------------------------------------------------#
# calculate kernel home ranges with reference method ####
#-------------------------------------------------------#

kernS <- kernelUD(data[,"id_yr"], h="href", grid=150, same4all=FALSE, kern="bivnorm")
conts <- getverticeshr(kernS, percent=95)   #get the contours as a spatial polygons dataframe

class(kernS)

#Plot the kernels, one individual at a time
image(kernS[[3]])
plot(conts[3,], add=TRUE)

#calculate contour areas from kernel UD
kernShr <- kernel.area(kernS, percent = c(50, 95, 99), unin="m", unout="km2")
kernShr

outname<- paste("RD2H_SummerSeasonalRanges_95KUD")

#write out the file, "E:/R" is the folder that you want, "test" is the output file name
writeOGR(conts,"D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs", outname, driver="ESRI Shapefile",overwrite=TRUE)

#Now export KUDs based on years

#Create year column
head(conts)
conts$year<-str_split_fixed(conts$id, "_", 2)[,2]
head(conts)

d.11<-subset(conts,conts$year=="2011")
d.12<-subset(conts,conts$year=="2012")
d.14<-subset(conts,conts$year=="2014")
d.16<-subset(conts,conts$year=="2016")
d.17<-subset(conts,conts$year=="2017")
d.18<-subset(conts,conts$year=="2018")
d.19<-subset(conts,conts$year=="2019")
d.20<-subset(conts,conts$year=="2020")

outname<- paste("RD2H_SummerSeasonalRanges_95KUD_2020")

#write out the file, "E:/R" is the folder that you want, "test" is the output file name
writeOGR(d.20,"D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs", outname, driver="ESRI Shapefile",overwrite=TRUE)
