
#----------------------#
# set up ####




#ATTN:: Skip to Reproducible below here chunk 


#caution need to run CalculatingFidelityToRouteAndStops lastest version first

rm(list=ls())

require(raster)
require(lubridate)
require(sp)
require(sf)
require(geoR)
require(dplyr)
require(rgeos)
require(tidyverse)
library(parallel)
library(terra)
library(stringr)
library(doSNOW)
#import data
#set your working drive



#bring in data and clip to the timing you need
load("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Migrations/RDH_AllMigrations_2014t2022.RData")
load("C:/Users/lwilde2/Desktop/RDH Database/RD2H_GPSCollarData/RD2H_CleanData.RData")


#determine which have >1 year
Spring.summary$AID <- str_sub(Spring.summary$id_yr,1,3)

keep <- Spring.summary %>% group_by(AID) %>% summarise(n = n_distinct(id_yr)) %>% ungroup() %>% left_join(Spring.summary, by = "AID") %>% filter(n > 1) 

rdh_data.sf <- rdh_data.sf %>% dplyr::filter(id_yr %in% keep$id_yr)





# rm(list=ls()[! ls() %in% c("irg.params", "indvi","dlcparams.wd")])

# load("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Stopover Fidelity/FidelityMetrics_04092023.RData")

#keep only those that are in the fidelity analysis

dataset <- rdh_data.sf %>% arrange(id_yr, POSIXct) #%>% dplyr::filter(id_yr %in% unique(so_fidelity$id_yr))
rm(rdh_data.sf)





#check
length(unique(dataset$id_yr)); length(unique(dataset$AID))


#create max.meters object


#ATTN: make sure this is referencing dist by march 21, not Jan 1 -- ie 80

dataset[which(dataset$JDate == 80),"Date"]

max.meters.vec <- data.frame(dist = NA, id_yr = NA)
for(i in 1:length(unique(dataset$id_yr))){ #length(unique(dataset$id_yr))
  #i = 1
tmp.data <- dataset %>% arrange(id_yr, POSIXct) %>% filter(id_yr == unique(dataset$id_yr)[i]) 
dist <- max(st_distance(x = tmp.data, y = tmp.data[which(tmp.data$JDate %in% 80),], by_element = F) %>% as.numeric())


t <- data.frame(dist = dist, id_yr = unique(dataset$id_yr)[i])

rm(tmp.data)

max.meters.vec <- rbind(max.meters.vec, t)

rm(dist, t)
}

names(max.meters.vec)[1] <- "max.meters"

max.meters.vec <- max.meters.vec[-1,]

head(max.meters.vec)
max.meters.all <- max(max.meters.vec$max.meters)

names(dataset)[names(dataset) == 'AID'] <- 'id'
names(dataset)[names(dataset) == 'POSIXct'] <- 'date'



dataset <- dataset[,c("id_yr","date","Year","geometry", "id")]




# add in long lat column
dataset[,c('x','y')] <- st_coordinates(dataset)



d <- dataset %>% st_drop_geometry() %>% dplyr::select(id_yr, date, x, y, id, Year) %>% as.data.frame()
head(d)
table(is.na(d))

# rename columns
names(d) <- c("id_yr", "date", "x", "y", "id", "year")

# order by id_yr
#d<-d[order(d$id_yr,d$date),]
head(d)

#

#check
length(unique(d$id_yr)); length(unique(d$id))

#mapview::mapview(datamd, zcol = "id_yr")


#----------------------------------------------------------------------#
#----------------------------------------------------------------------#

loc <- "C:/Users/lwilde2/Desktop/RDH Database/RedDesertToHobackMuleDeerResearchProject/RCode/Surfing/Greenscape_SurfingScores/Surfing_Greeenscape_functions.R"
source(loc)

#data is here
d

#IRG NDVI are here
#irg.params
#indvi

r <- rast("Z:/MODIS_NDVI/NDVI2022/DLC_Estimates/DLC_2001_xmidS_Estimate.tif")
r <- raster(r)
#-----------------------------------------#
# Begin the function ####


# START #

 colPalIRG<-colorRampPalette(c("#8c510a", "#bf812d", "#dfc27d",
                               "#f6e8c3", "#f5f5f5", "#c7eae5",
                               "#80cdc1", "#35978f", "#01665e"))(210)

 metaSurf<-data.frame(AnimalID="",
                      Year=2017,
                      AY_ID="",
                      IRG_120=0.9,
                      DFP_120=5.5,
                      IRG_90=0.9,
                      DFP_90=5.5,
                      IRG_60=0.9,
                      DFP_60=5.5,
                      svSlope=1,
                      greenUpDur=34.5,
                      springScale=15.5,
                      SVmax="",
                      midpoint=100)[0, ]



 #par(mfrow=c(2, 2))


#Clarify projected coordinate system
str(d)

#remove duplicate locations
d <- d[ !duplicated(c("x", "y")),]

#remove duplicate times
d <- d[ !duplicated(c("id_yr","date"))]

#d1 <- d %>% filter(id_yr %in% unique(.$id_yr))#202


#check if any animals died in  May or June
#SummerDeath <- d %>% as.data.frame() %>% dplyr::group_by(id_yr) %>% dplyr::summarise(min = min(date), max = max(date)) %>% ungroup() %>% mutate(last.month = month(max)) %>% dplyr::filter(last.month %in% c(5,6))

#d <- d %>% as.data.frame() %>% dplyr::filter(!id_yr %in% SummerDeath$id_yr)

#d <- as.data.frame(d)
d <- st_as_sf(d, coords = c("x", "y"), crs = st_crs(mig.SpringMigration.sf))

projIRG<-st_crs(r)$input

d <- st_transform(d, crs=projIRG)


gps.data<-as(d, "Spatial")
gps.data$Timestamp<-gps.data$date


class(gps.data)


#set cut-off -- Ellen found this to be optimal
cut.off<-20

#define the parameters
gps.data$AY_ID <- gps.data$id_yr
AY_ID<-unique(gps.data$AY_ID)

head(gps.data)
length(unique(gps.data$id_yr)); length(unique(gps.data$id))

unique(AY_ID)
unique(gps.data$AY_ID)

#------------------------------#
# create IRG layers ####
# setwd("Z:/MODIS_NDVI/NDVI2022/DLC_Estimates/")
# xmids14 <- raster(paste(getwd(),"/DLC_2014_xmidS_Estimate.tif", sep = ""))
# xmids15 <- raster(paste(getwd(),"/DLC_2015_xmidS_Estimate.tif", sep = ""))
# xmids16 <- raster(paste(getwd(),"/DLC_2016_xmidS_Estimate.tif", sep = ""))
# xmids17 <- raster(paste(getwd(),"/DLC_2017_xmidS_Estimate.tif", sep = ""))
# xmids18 <- raster(paste(getwd(),"/DLC_2018_xmidS_Estimate.tif", sep = ""))
# xmids19 <- raster(paste(getwd(),"/DLC_2019_xmidS_Estimate.tif", sep = ""))
# xmids20 <- raster(paste(getwd(),"/DLC_2020_xmidS_Estimate.tif", sep = ""))
# xmids21 <- raster(paste(getwd(),"/DLC_2021_xmidS_Estimate.tif", sep = ""))
# xmids22 <- raster(paste(getwd(),"/DLC_2022_xmidS_Estimate.tif", sep = ""))
# 
# xmidsAll <- stack(xmids14,xmids15,xmids16,xmids17,xmids18,xmids19,xmids20,xmids21,xmids22)
# 
# scalS14 <- raster(paste(getwd(),"/DLC_2014_scalS_Estimate.tif", sep = ""))
# scalS15 <- raster(paste(getwd(),"/DLC_2015_scalS_Estimate.tif", sep = ""))
# scalS16 <- raster(paste(getwd(),"/DLC_2016_scalS_Estimate.tif", sep = ""))
# scalS17 <- raster(paste(getwd(),"/DLC_2017_scalS_Estimate.tif", sep = ""))
# scalS18 <- raster(paste(getwd(),"/DLC_2018_scalS_Estimate.tif", sep = ""))
# scalS19 <- raster(paste(getwd(),"/DLC_2019_scalS_Estimate.tif", sep = ""))
# scalS20 <- raster(paste(getwd(),"/DLC_2020_scalS_Estimate.tif", sep = ""))
# scalS21 <- raster(paste(getwd(),"/DLC_2021_scalS_Estimate.tif", sep = ""))
# scalS22 <- raster(paste(getwd(),"/DLC_2022_scalS_Estimate.tif", sep = ""))
# 
# scalSAll <- stack(scalS14,scalS15,scalS16,scalS17,scalS18,scalS19,scalS20,scalS21,scalS22)
# 
# xmids_SD14 <- raster(paste(getwd(),"/DLC_2014_xmids_SD_Estimate.tif", sep = ""))
# xmids_SD15 <- raster(paste(getwd(),"/DLC_2015_xmids_SD_Estimate.tif", sep = ""))
# xmids_SD16 <- raster(paste(getwd(),"/DLC_2016_xmids_SD_Estimate.tif", sep = ""))
# xmids_SD17 <- raster(paste(getwd(),"/DLC_2017_xmids_SD_Estimate.tif", sep = ""))
# xmids_SD18 <- raster(paste(getwd(),"/DLC_2018_xmids_SD_Estimate.tif", sep = ""))
# xmids_SD19 <- raster(paste(getwd(),"/DLC_2019_xmids_SD_Estimate.tif", sep = ""))
# xmids_SD20 <- raster(paste(getwd(),"/DLC_2020_xmids_SD_Estimate.tif", sep = ""))
# xmids_SD21 <- raster(paste(getwd(),"/DLC_2021_xmids_SD_Estimate.tif", sep = ""))
# xmids_SD22 <- raster(paste(getwd(),"/DLC_2022_xmids_SD_Estimate.tif", sep = ""))
# 
# xmids_SDAll <- stack(xmids_SD14,xmids_SD15,xmids_SD16,xmids_SD17,xmids_SD18,xmids_SD19,xmids_SD20,xmids_SD21,xmids_SD22)
# 
# setwd("Z:/MODIS_NDVI/NDVI2022/Bischof_calculations/")
# indvi14 <- raster(paste(getwd(),"/MOD09Q1_2014_IntegratedNDVI.tif", sep = ""))
# indvi15 <- raster(paste(getwd(),"/MOD09Q1_2015_IntegratedNDVI.tif", sep = ""))
# indvi16 <- raster(paste(getwd(),"/MOD09Q1_2016_IntegratedNDVI.tif", sep = ""))
# indvi17 <- raster(paste(getwd(),"/MOD09Q1_2017_IntegratedNDVI.tif", sep = ""))
# indvi18 <- raster(paste(getwd(),"/MOD09Q1_2018_IntegratedNDVI.tif", sep = ""))
# indvi19 <- raster(paste(getwd(),"/MOD09Q1_2019_IntegratedNDVI.tif", sep = ""))
# indvi20 <- raster(paste(getwd(),"/MOD09Q1_2020_IntegratedNDVI.tif", sep = ""))
# indvi21 <- raster(paste(getwd(),"/MOD09Q1_2021_IntegratedNDVI.tif", sep = ""))
# indvi22 <- raster(paste(getwd(),"/MOD09Q1_2022_IntegratedNDVI.tif", sep = ""))
# 
# indviAll <- stack(indvi14,indvi15,indvi16,indvi17,indvi18,indvi19,indvi20,indvi21,indvi22)
# 
# indviAll[[paste0("X", "2014")]]
# 
# 
# load("G:/Shared drives/wyo-coop-wilde/NDVI_Calcs/Aikens_greenscape/dlc_params_2014to2022.RData")
# 
# #create blank raster to downsize rasters
# r <- raster(nrow=dim(xmids_SDAll)[1]/3, ncol=dim(xmids_SDAll)[2]/3, ) #downsize raster by 1.5
# extent(r) <- extent(xmids_SDAll)
# projection(r) <- st_crs(xmids_SDAll)$input
# 
# xmids_SDAll_lite <- resample(xmids_SDAll, r, method = "bilinear")
# xmidsAll_lite <- resample(xmidsAll, r, method = "bilinear")
# scalSAll_lite <- resample(scalSAll, r, method = "bilinear")
# indviAll_lite <- resample(indviAll, r, method = "bilinear")
# 

# save(indviAll_lite, xmidsAll_lite, xmids_SDAll_lite, scalSAll_lite, file = "G:/Shared drives/wyo-coop-wilde/NDVI_Calcs/Aikens_greenscape/dlc_params_2014to2022_lite.RData")

# xmids_SDAll_lite <- projectRaster(xmids_SDAll_lite, crs=CRS("+proj=longlat +ellps=WGS84"))
# xmidsAll_lite <- projectRaster(xmidsAll_lite, crs=CRS("+proj=longlat +ellps=WGS84"))
# scalSAll_lite <- projectRaster(scalSAll_lite, crs=CRS("+proj=longlat +ellps=WGS84"))
# indviAll_lite <- projectRaster(indviAll_lite, crs=CRS("+proj=longlat +ellps=WGS84"))
# 
# save(indviAll_lite, xmidsAll_lite, xmids_SDAll_lite, scalSAll_lite, file = "G:/Shared drives/wyo-coop-wilde/NDVI_Calcs/Aikens_greenscape/dlc_params_2014to2022_lite_latlon.RData")



#save.image("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Greenscape/InputDataExcept_20230414.RData")


#this is where the bischof calc layers are stored!


#--------------------------------------------#
# Reproducible below here ####

#this is where the bischof calc layers are stored!
BischCalcs <- "G:/Shared drives/wyo-coop-wilde/NDVI_Calcs/Aikens_greenscape/dlc_params_2014to2022.RData" #_lite for subsampled by 3

BischCalcs.ll <- "G:/Shared drives/wyo-coop-wilde/NDVI_Calcs/Aikens_greenscape/dlc_params_2014to2022_lite_latlon.RData"


#these are the data

load("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Greenscape/InputData_20230414.RData")


load(BischCalcs) 



# check on animals that it didnt return == did it have an error? Run it
# mapview in all projections and see that it is where you think it is relative to the rasters too

#maybe its skipping animals


library(mapview)
#mapview(st_as_sf(gps.data.ay))



#------------------------------------#
# #parallel version - turns out there is a problem running parallel with Rasters bc it can struggle to recollect tmp files ####
# library(parallel)
# library(snow)
# library(doParallel)
# library(doSNOW)
# library(doRNG)
# 
# #Create unique list of animals
# #missed <- gps.data %>% as.data.frame() %>% filter(!AY_ID %in% metasurf_all$AY_ID)
# id.list <- AY_ID#unique(missed$AY_ID)#unique(gps.data$AY_ID)
# names(id.list) <- NULL
# row.names(id.list) <- NULL
# 
# #Loop over animals and extract greenscape metrics
# #First, set up parallel processing
# no_cores <- detectCores() - 1
# cat(no_cores, " cores detected.")
# 
# cl <- makeSOCKcluster(no_cores)
# registerDoSNOW(cl) #Use doSNOW instead of doParallel to create progress bar
# 
# #Set up progress bar parameters and function
# #Progress bar will clarify the percentage of complete BBMMs
# pb <- txtProgressBar(min=1, max=length(AY_ID), style=3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress=progress)
# 
# showConnections()
# 
# gps.data <- gps.data %>% st_as_sf() %>% dplyr::filter(AY_ID %in% AY_ID) %>% as(., "Spatial")
# 
# metasurf_miss2 <- foreach(i = 1:length(AY_ID), #1:length(AY_ID)
#                          .packages = c('raster','rgeos','sf','sp',
#                                        'rgdal','stringr','geoR','lubridate','dplyr','tidyverse','terra'),
#                          .combine = 'rbind',
#                          .errorhandling='remove',
#                          .options.multicore=list(set.seed=TRUE),
#                          .options.snow=opts
# 
# ) %dorng% {
# 
# #i=135 #check on AY_ID[15]
# source(loc)
# 
# 
#   year<-strsplit(as.character(AY_ID[i]), "_")[[1]][2]
#   ny <- as.numeric(strsplit(as.character(AY_ID[i]), "_")[[1]][2])-2013
#   #yearn <- as.numeric(year)
# 
#   #if(year>=2001 & year<=2022){
# 
#   load(BischCalcs) #here is where you can load in the raster stacks
# 
#   #get IRG data
#   irg.max.day<-xmidsAll[[ny]]
#   spring.scale<-scalSAll[[ny]]
#   sd.irgDate<-xmids_SDAll[[ny]]
#   ndvi<-indviAll[[ny]]
#   #names(ndvi)
#   irg.max.day #inspect
# 
#   #rm(xmidsAll,scalSAll,xmids_SDAll,indviAll) #remove to clear memory space
# 
#   #subset gps data
#   gps.data.ay <- gps.data[gps.data$AY_ID==AY_ID[i], ]
# 
# 
# 
#   #extract the phenology data
# 
#   sd.vals<-raster::extract(x=sd.irgDate, gps.data.ay)
#   irg.vals<-raster::extract(x=irg.max.day, gps.data.ay)
#   scale.vals<-raster::extract(x=spring.scale, gps.data.ay)
#   indvi.vals<-raster::extract(x=ndvi, gps.data.ay)
# 
#   rm(ndvi)
# 
# 
# 
#   #filter out really poor fitting IRG curves
#   idx<-which(sd.vals>cut.off)
#   if(length(idx)>0){
#     gps.data.ay<-gps.data.ay[-idx, ]
#     sd.vals<-sd.vals[-idx]
#     irg.vals<-irg.vals[-idx]
#     scale.vals<-scale.vals[-idx]
#     indvi.vals<-indvi.vals[-idx]
#   }
# 
#   #reproject to lat long, important for semivariogram
#   gps.data.ay.ll <- spTransform(gps.data.ay,CRS("+proj=longlat +ellps=WGS84 +no_defs"))
# 
#   #calculate greenscape
#   greenscapeVals<-getSvGreenscapeVals(spDF=gps.data.ay.ll,
#                                       irg.vals=irg.vals,
#                                       scale.vals=scale.vals,
#                                       PlotVariogram=F,
#                                       saveJPEG=F,
#                                       dir="C:/Users/lwilde2/Documents/AdaptiveFidelity/Outputs")
# 
#   #calculate surfing
#   mid.i<-greenscapeVals[6]
#   #used points greenscape
# 
#   surfingVals120<-getSurfingMetricsVals(spDF=gps.data.ay.ll,
#                                         irg.vals=irg.vals,
#                                         scale.vals=scale.vals,
#                                         startJDay=mid.i-60, endJDay=mid.i+60)
# 
# 
#   #make "availability" buffer
#   # sp.buffer<-max.meters[which(names(max.meters)==species)]
#   # sp.buffer<-max.meters.vec[which(max.meters.vec$id_yr == AY_ID[i]),"max.meters"] #could make this iterate for each animal
#   # centroid<-apply(coordinates(gps.data.ay.ll), MARGIN=2, FUN=mean, na.rm=T)
#   # xy<-data.frame(x=centroid[1], y=centroid[2])
#   # row.names(xy)<-AY_ID[i]
#   # centroid.sp<-SpatialPointsDataFrame(SpatialPoints(xy), data.frame(ID=as.character(AY_ID[i])))
#   # proj4string(centroid.sp)<-CRS("+proj=longlat +ellps=WGS84")
#   # centroid.sp<-spTransform(centroid.sp, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
#   # circ.sp<-gBuffer(centroid.sp, width=sp.buffer, byid=F)
#   # circ.sp<-spTransform(circ.sp, st_crs(irg.max.day)$input)  #here bring it back to sinu projection
#   #was CRS("+proj=longlat +ellps=WGS84")
# 
#   # #create blank raster to downsize rasters
#   # r <- raster(nrow=dim(irg.max.day)[1]/3, ncol=dim(irg.max.day)[2]/3) #downsize raster by 1.5
#   # extent(r) <- extent(irg.max.day)
#   # projection(r) <- st_crs(irg.max.day)$input
#   #
#   # #each takes 34 seconds!!
#   # #whole code below takes 20 min
#   #
#   # #landscape-scale greenscape
#   # system.time( sp.greensscape<-getSvGreenscapeCirle(circles=circ.sp,
#   #                                                   irgRaster=resample(irg.max.day, r, method = "bilinear"),
#   #                                                   scaleRaster=resample(spring.scale, r, method = "bilinear"),
#   #                                                   sdRaster=resample(sd.irgDate, r, method = "bilinear"), sdCutoff=cut.off,
#   #                                                   PlotVariogram=F, plotColors=colPalIRG, saveJPEG=F,
#   #                                                   dir="C:/Users/lwilde2/Documents/AdaptiveFidelity/Outputs")) #make sure the plot and save are FALSE
# 
# 
#   #put it all together and save it
#   # add in available greenscape metrics in the dataframe, once species-specific bufers are calculated
#   metaSurf.i<-data.frame(AnimalID=gps.data.ay$id[1],
#                          Year=gps.data.ay$year[1],
#                          AY_ID=AY_ID[i],
#                          IRG_120=surfingVals120[1],
#                          DFP_120=surfingVals120[2],
#                          svSlope=greenscapeVals[1],
#                          greenUpDur=greenscapeVals[4],
#                          springScale=greenscapeVals[5],
#                          SVmax=greenscapeVals[14],
#                          midpoint=mid.i)
# 
#   #maybe remove this
#   #metaSurf<-metaSurf.i
# 
#   #return(rbind(metaSurf, metaSurf.i))
# }
# 
# #stop the cluster
# stopCluster(cl)
# showConnections()


#save(metasurf_miss, file = "C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Greenscape/Output/Greenscape_Calcs_1of5.RData")



#This is how it means

#IRG_120 and DFP_120 -- the average, DFP is absolute
# svSlope -- is order (higher means more ordered; 0-1 bounded; Spearmans rank correlation)
# greenupDur -- duration, 37 days 
#midpoint -- median date of IRG, what all is based off
# SVmax -- higher number is heterogeneity? (need to look into it)
# springScale -- greenup rate; inverse for graphing; higher number is more rapid; rapid = wavelike

#save.image("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/GreenscapeOutput_02152023.RData")


#the way to graph is similar to - longer, rapid, consecutive leads to higher foraging penalty, better surfers and predicts animals that migrate



#----------------------------------#
# Start For Loop ####

load("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Greenscape/InputDataExcept_20230414.RData")






#memory allocation issues #
#-----------------------------------#
#clear memory space
gc()

rm(list=ls()[! ls() %in% c("d", "dataset", "gps.data", "AY_ID", "metaSurf.i", "metaSurf", "max.meters.vec", "cut.off","loc","BischCalcs","BischCalcs.ll","colPalIRG","max.meters.all")])
gc()
source(loc)

load("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Greenscape/Output/Greenscape_Calcs_1of5.RData")
load("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Greenscape/Output/Greenscape_Calcs_2of5.RData")
load("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Greenscape/Output/Greenscape_Calcs_3of5.RData")

#meta1 <- metaSurf
rm(metaSurf)
AY_ID2 <- AY_ID[which(!AY_ID %in% c(metasurf_miss$AY_ID, meta1$AY_ID, meta_miss2$AY_ID))]

#meta_miss2 <- meta_miss2 %>% mutate(AY_ID = paste(AnimalID,"_",Year,sep=""))

metaSurf<-data.frame(AnimalID="",
                     Year=2017,
                     AY_ID="",
                     IRG_120=0.9,
                     DFP_120=5.5,
                     IRG_90=0.9,
                     DFP_90=5.5,
                     IRG_60=0.9,
                     DFP_60=5.5,
                     svSlope=1,
                     greenUpDur=34.5,
                     springScale=15.5,
                     SVmax="",
                     midpoint=100)[0, ]

for(i in 1:length(AY_ID2)){#
#i = 10

tryCatch({

# dlcparams.wd<-"Z:/MODIS_NDVI/Bischof_calculations"
# setwd(dlcparams.wd)

year<-strsplit(as.character(AY_ID2[i]), "_")[[1]][2]
ny <- as.numeric(strsplit(as.character(AY_ID2[i]), "_")[[1]][2])-2013
#yearn <- as.numeric(year)

#if(year>=2001 & year<=2022){

load(BischCalcs) #here is where you can load in the raster stacks

#get IRG data
irg.max.day<-xmidsAll[[ny]]
spring.scale<-scalSAll[[ny]]
sd.irgDate<-xmids_SDAll[[ny]]
ndvi<-indviAll[[ny]]
#names(ndvi)
irg.max.day #inspect

#rm(xmidsAll,scalSAll,xmids_SDAll,indviAll) #remove to clear memory space

#subset gps data 
gps.data.ay <- gps.data[gps.data$AY_ID==AY_ID2[i], ]



#extract the phenology data 

sd.vals<-raster::extract(x=sd.irgDate, gps.data.ay)
irg.vals<-raster::extract(x=irg.max.day, gps.data.ay)
scale.vals<-raster::extract(x=spring.scale, gps.data.ay)
indvi.vals<-raster::extract(x=ndvi, gps.data.ay)

rm(ndvi)



#filter out really poor fitting IRG curves
idx<-which(sd.vals>cut.off) #this removes when the fits are crappy, sensitivity analysis says this was best for most
if(length(idx)>0){
  gps.data.ay<-gps.data.ay[-idx, ]
  sd.vals<-sd.vals[-idx]
  irg.vals<-irg.vals[-idx]
  scale.vals<-scale.vals[-idx]
  indvi.vals<-indvi.vals[-idx]
}

#reproject to lat long, important for semivariogram
gps.data.ay.ll <- spTransform(gps.data.ay,CRS("+proj=longlat +ellps=WGS84 +no_defs"))

#calculate greenscape 
greenscapeVals<-getSvGreenscapeVals(spDF=gps.data.ay.ll, 
                                    irg.vals=irg.vals, 
                                    scale.vals=scale.vals,
                                    PlotVariogram=F, 
                                    saveJPEG=F, 
                                    dir="C:/Users/lwilde2/Documents/AdaptiveFidelity/Outputs")

#calculate surfing 
mid.i<-greenscapeVals[6] #midpoint of date of peak green-up, allows you to set the period
#midpoint is related to the distribution of peak green-up, not the double logistic curves

#used points greenscape 

surfingVals120<-getSurfingMetricsVals(spDF=gps.data.ay.ll, 
                                      irg.vals=irg.vals, 
                                      scale.vals=scale.vals, 
                                      startJDay=mid.i-60, endJDay=mid.i+60)


#make "availability" buffer
# sp.buffer<-max.meters[which(names(max.meters)==species)]
# sp.buffer<-max.meters.vec[which(max.meters.vec$id_yr == AY_ID[i]),"max.meters"] #could make this iterate for each animal
# centroid<-apply(coordinates(gps.data.ay.ll), MARGIN=2, FUN=mean, na.rm=T)
# xy<-data.frame(x=centroid[1], y=centroid[2])
# row.names(xy)<-AY_ID[i]
# centroid.sp<-SpatialPointsDataFrame(SpatialPoints(xy), data.frame(ID=as.character(AY_ID[i])))
# proj4string(centroid.sp)<-CRS("+proj=longlat +ellps=WGS84")
# centroid.sp<-spTransform(centroid.sp, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
# circ.sp<-gBuffer(centroid.sp, width=sp.buffer, byid=F)
# circ.sp<-spTransform(circ.sp, st_crs(irg.max.day)$input)  #here bring it back to sinu projection
# #was CRS("+proj=longlat +ellps=WGS84")


#landscape-scale greenscape
#you only need this if -- to assess whether a landscape should select for residency or migration beh, when you just use the above its just what the animal was experiencing
#asks what is the greenscape of this landscape in this year (could do this for just each year once) -- if you wanted to look at landscape scale diffs across pop or years then control but just need it once per year
#because the RDH use all the same route, its better to just do one

#create blank raster to downsize rasters
# r <- raster(nrow=dim(irg.max.day)[1]/3, ncol=dim(irg.max.day)[2]/3, ) #downsize raster by 1.5
# extent(r) <- extent(irg.max.day)
# projection(r) <- st_crs(irg.max.day)$input

# system.time( sp.greensscape<-getSvGreenscapeCirle(circles=circ.sp, 
#                                      irgRaster=resample(irg.max.day, r, method = "bilinear"), 
#                                      scaleRaster=resample(spring.scale, r, method = "bilinear"),
#                                      sdRaster=resample(sd.irgDate, r, method = "bilinear"), sdCutoff=cut.off,
#                                      PlotVariogram=F, plotColors=colPalIRG, saveJPEG=F, 
#                                      dir="C:/Users/lwilde2/Documents/AdaptiveFidelity/Outputs")) #make sure the plot and save are FALSE


#put it all together and save it 
# add in available greenscape metrics in the dataframe, once species-specific bufers are calculated 
metaSurf.i<-data.frame(AnimalID=gps.data.ay$id[1], 
                       Year=gps.data.ay$year[1], 
                       AY_ID=AY_ID2[i], 
                       IRG_120=surfingVals120[1], 
                       DFP_120=surfingVals120[2], 
                       svSlope=greenscapeVals[1], 
                       greenUpDur=greenscapeVals[4], 
                       springScale=greenscapeVals[5], 
                       SVmax=greenscapeVals[14], 
                       midpoint=mid.i)



metaSurf<-rbind(metaSurf, metaSurf.i)


}, error=function(e){})


}#i


head(metaSurf)

meta3 <- metaSurf

#save(meta3, file = "C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/Processed Data/Greenscape/Output/Greenscape_Calcs_4of5.RData")


metaSurf

sessionInfo()

# 
# no_cores <- detectCores()-1
# AY_ID <- AY_ID
# 
# # Setup cluster
# clust <- makeCluster(no_cores)
# 
# # export the objects you need in the loop
# clusterExport(clust, varlist=c("d", "dataset", "gps.data", "AY_ID", "max.meters.vec", "cut.off","BischCalcs","colPalIRG"))
# showConnections()#check
# 
# ## run the parallel
# green.out.120 <- do.call(rbind, clusterApplyLB(clust, c(1,25,33,54), function(i){
# 
#   #i = 3
#   # dlcparams.wd<-"Z:/MODIS_NDVI/Bischof_calculations"
#   # setwd(dlcparams.wd)
#   
#   year<-strsplit(as.character(AY_ID[i]), "_")[[1]][2]
#   ny <- as.numeric(strsplit(as.character(AY_ID[i]), "_")[[1]][2])-2013
#   #yearn <- as.numeric(year)
#   
#   if(year>=2001 & year<=2022){
#     
#     load(BischCalcs) #here is where you can load in the raster stacks
#     
#     #get IRG data
#     irg.max.day<-xmidsAll[[ny]]
#     spring.scale<-scalSAll[[ny]]
#     sd.irgDate<-xmids_SDAll[[ny]]
#     ndvi<-indviAll[[ny]]
#     #names(ndvi)
#     irg.max.day #inspect
#     
#     rm(xmidsAll,scalSAll,xmids_SDAll,indviAll) #remove to clear memory space
#     
#     #subset gps data 
#     gps.data.ay <- gps.data[gps.data$AY_ID==AY_ID[i], ]
#     
#     
#     
#     #extract the phenology data 
#     
#     sd.vals<-raster::extract(x=sd.irgDate, gps.data.ay)
#     irg.vals<-raster::extract(x=irg.max.day, gps.data.ay)
#     scale.vals<-raster::extract(x=spring.scale, gps.data.ay)
#     indvi.vals<-raster::extract(x=ndvi, gps.data.ay)
#     
#     rm(ndvi)
#     
#     #reproject to lat long, important for semivariogram
#     gps.data.ay.ll <- spTransform(gps.data.ay,CRS("+proj=longlat +ellps=WGS84 +no_defs"))
#     
#     #filter out really poor fitting IRG curves
#     idx<-which(sd.vals>cut.off)
#     if(length(idx)>0){
#       gps.data.ay<-gps.data.ay[-idx, ]
#       sd.vals<-sd.vals[-idx]
#       irg.vals<-irg.vals[-idx]
#       scale.vals<-scale.vals[-idx]
#       indvi.vals<-indvi.vals[-idx]
#     }
#     
#     
#     #calculate greenscape 
#     greenscapeVals<-getSvGreenscapeVals(spDF=gps.data.ay.ll, 
#                                         irg.vals=irg.vals, 
#                                         scale.vals=scale.vals,
#                                         PlotVariogram=F, 
#                                         saveJPEG=F, 
#                                         dir="C:/Users/lwilde2/Documents/AdaptiveFidelity/Outputs")
#     
#     #calculate surfing 
#     mid.i<-greenscapeVals[6]
#     #used points greenscape 
#     
#     surfingVals120<-getSurfingMetricsVals(spDF=gps.data.ay.ll, 
#                                           irg.vals=irg.vals, 
#                                           scale.vals=scale.vals, 
#                                           startJDay=mid.i-60, endJDay=mid.i+60)
#     
#     
#     #make "availability" buffer
#     # sp.buffer<-max.meters[which(names(max.meters)==species)]
#     sp.buffer<-max.meters.vec[which(max.meters.vec$id_yr == AY_ID[i]),"max.meters"] #could make this iterate for each animal
#     centroid<-apply(coordinates(gps.data.ay.ll), MARGIN=2, FUN=mean, na.rm=T)
#     xy<-data.frame(x=centroid[1], y=centroid[2])
#     row.names(xy)<-AY_ID[i]
#     centroid.sp<-SpatialPointsDataFrame(SpatialPoints(xy), data.frame(ID=as.character(AY_ID[i])))
#     proj4string(centroid.sp)<-CRS("+proj=longlat +ellps=WGS84")
#     centroid.sp<-spTransform(centroid.sp, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
#     circ.sp<-gBuffer(centroid.sp, width=sp.buffer, byid=F)
#     circ.sp<-spTransform(circ.sp, st_crs(irg.max.day)$input)  #here bring it back to sinu projection
#     #was CRS("+proj=longlat +ellps=WGS84")
#     
#     #create blank raster to downsize rasters
#     r <- raster(nrow=dim(irg.max.day)[1]/3, ncol=dim(irg.max.day)[2]/3, ) #downsize raster by 1.5
#     extent(r) <- extent(irg.max.day)
#     projection(r) <- st_crs(irg.max.day)$input
#     
#     #each takes 34 seconds!!
#     #whole code below takes 20 min
#     
#     #landscape-scale greenscape
#     system.time( sp.greensscape<-getSvGreenscapeCirle(circles=circ.sp, 
#                                                       irgRaster=resample(irg.max.day, r, method = "bilinear"), 
#                                                       scaleRaster=resample(spring.scale, r, method = "bilinear"),
#                                                       sdRaster=resample(sd.irgDate, r, method = "bilinear"), sdCutoff=cut.off,
#                                                       PlotVariogram=F, plotColors=colPalIRG, saveJPEG=F, 
#                                                       dir="C:/Users/lwilde2/Documents/AdaptiveFidelity/Outputs")) #make sure the plot and save are FALSE
#     
#     
#     #put it all together and save it 
#     # add in available greenscape metrics in the dataframe, once species-specific bufers are calculated 
#     metaSurf.i<-data.frame(AnimalID=gps.data.ay$id[1], 
#                            Year=gps.data.ay$year[1], 
#                            AY_ID=AY_ID[i], 
#                            IRG_120=surfingVals120[1], 
#                            DFP_120=surfingVals120[2], 
#                            svSlope=greenscapeVals[1], 
#                            greenUpDur=greenscapeVals[4], 
#                            springScale=greenscapeVals[5], 
#                            SVmax=greenscapeVals[14], 
#                            midpoint=mid.i)
#     
#     
#     
#     # metaSurf<-rbind(metaSurf, metaSurf.i)
#     
#     gc()
#     rm(list=ls()[! ls() %in% c("d", "dataset", "gps.data", "AY_ID", "metaSurf.i", "metaSurf", "max.meters.vec", "cut.off","loc","BischCalcs","colPalIRG")])
#     gc()
# 
# }#year filter
# }#loop
# ))
# 
# stopCluster(clust)
# showConnections()
# 
# green.out.120




#This is how it means

#IRG_120 and DFP_120 -- the average, DFP is absolute
# svSlope -- is order (higher means more ordered; 0-1 bounded; Spearmans rank correlation)
# greenupDur -- duration, 37 days 
#midpoint -- median date of IRG, what all is based off
# SVmax -- higher number is heterogeneity? (need to look into it)
# springScale -- greenup rate; inverse for graphing; higher number is more rapid; rapid = wavelike

#save.image("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/GreenscapeOutput_02152023.RData")


#the way to graph is similar to - longer, rapid, consecutive leads to higher foraging penalty, better surfers and predicts animals that migrate



