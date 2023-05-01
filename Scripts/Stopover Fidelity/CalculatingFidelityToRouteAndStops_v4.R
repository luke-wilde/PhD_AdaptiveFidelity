# as of 04/09/23: still having data dropout, where NA for the stop numbers is causing animals to drop (21 id_yr's).


#---------------------------------#
#---------------------------------#
#per Anna and Van der Kerk 2021, stop is more than 1 d and less than 23
#Anna defines high use if >3 d

# onstop_yr_fin <- onstop_yr_fin %>% filter(TimeStopped <= (24*23)) #TimeStopped >= (24*3) & 
#---------------------------------#
#---------------------------------#

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

#load functions
source("C:/Users/lwilde2/Documents/AdaptiveFidelity/Scripts/Quantifying Fidelity/CalculatingFidelityToRoute_function.R")



#load in the datasets
# Stopovers <- st_read("Data/RD2H_HighUseStopovers_Polygons.shp")
load("C:/Users/lwilde2/Desktop/RDH Database/Processed Data/RDH_AllMigrations_Bischof_2014t2022_20230425.RData")

load("C:/Users/lwilde2/Desktop/RDH Database/Processed Data/RDH_AllStopovers_14t22_20230425.RData")

sa <- st_read("C:/Users/lwilde2/Documents/Chapter2_StopoverFidelity/Chapter1_StopoverFidelity/Data/GIS/StudyArea_Polygon_TableMountains.kml") #_TableMountains if you want to include the first SO

setwd("C:/Users/lwilde2/Desktop/RDH Database/RD2H_GPSCollarData/")
DeerID <- read.table("RD2H_AllDeerIDs_17APRIL2023.csv", header=TRUE, sep=",")
#Format date and time columns
DeerID$StartDate <- as.POSIXct(DeerID$StartDate, format="%m/%d/%Y", origin = "1970-01-01")
DeerID$StopDate <- as.POSIXct(DeerID$StopDate,format="%m/%d/%Y", origin = "1970-01-01")

#prep the data

#keep only long migrations for now, will change later
DeerID_long <- DeerID %>% filter(Mgtry == "Long") %>% distinct(AID, Mgtry, .keep_all = T) %>% arrange(AID)
data <- mig.SpringMigration.sf %>% filter(AID %in% DeerID_long$AID)


#data <- st_as_sf(x = data, coords = c("UTM_E", "UTM_N"), crs = "+proj=utm +zone=12 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

#get sa to the same
sa <- st_transform(sa, crs = st_crs(data)$input)

#data <- data %>% st_transform(crs = "+proj=utm +zone=12 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

#create test data
#deer108 <- data %>% filter(AID == "108") %>% mutate(id_yr = paste(AID,"_",Year,sep=""))

data <- st_crop(data, sa)

#which animals with >1 year, keep
keep <- data %>% st_drop_geometry() %>% group_by(AID) %>% summarize(n = n_distinct(Year)) %>% ungroup() %>% filter(n>1)

length(unique(DeerID_long$AID))
length(unique(data$AID))
length(unique(keep$AID))


data <- data %>% filter(AID %in% keep$AID)



#Calculate distance from winter range for each point using NSDs ####
#Create nsd object
d <- data.frame()
for(i in 1:length(unique(data$id_yr))){
  #i = 15
  df <- data %>% filter(id_yr == unique(data$id_yr)[i])
  df <- df %>% arrange(POSIXct) %>% mutate(km_mark = as.numeric(st_distance(x = df[,], y = df[1,], by_element = F)/1000))
  d <- rbind(d, df)
}

data <- d

nrow(data); length(unique(data$AID)); length(unique(data$id_yr))


#------------------------------------------------------#
#### Stopover Status ####

# writes a series of for loops to go through create a blank dataframe and iterate through each id to create a matrix, then can extract what proportion of each of the years points on its own migration stopovers are in the other ones too

#make the point layer and the SO layer within same crs
SO_all <- st_transform(SO_all, crs = st_crs(data)$input)

SO_all_sep <- SO_all %>% mutate(AID = str_sub(id_yr, 0,3)) %>% filter(layer == 1) %>% st_cast(., "POLYGON")

SO_all_sep <- SO_all_sep %>% group_by(id_yr) %>% mutate(stop.n = row_number(), year = as.numeric(str_sub(id_yr,5,8))) %>% ungroup()

length(unique(SO_all_sep$AID)); length(unique(SO_all_sep$id_yr))

table(unique(data$id_yr) %in% unique(SO_all_sep$id_yr))

check <- SO_all_sep[which(is.na(SO_all_sep$stop.n)),]
length(unique(check$AID)); length(unique(check$id_yr))

#get rid of animals missing SO info
missingSO <- data %>% filter(!id_yr %in% SO_all_sep$id_yr) %>% arrange(id_yr) %>% st_drop_geometry() %>% dplyr::select(id_yr)
unique(missingSO$id_yr)


data1 <- data %>% filter(!id_yr %in% missingSO$id_yr) %>% arrange(id_yr)


data1 <- data1[!duplicated(data1[, c("POSIXct", "AID")]), ]
nrow(data1); length(unique(data1$AID)); length(unique(data1$id_yr))




# Step 1: Extract points on stops within same year, for loop

onstop <- data.frame()

for(i in 1:length(unique(data1$id_yr))){
  #i = 2
  #
  pt <- data1 %>% filter(id_yr %in% unique(data1$id_yr)[i])
  st <- SO_all_sep %>% filter(id_yr %in% unique(data1$id_yr)[i])
  
  #extract points that are inside a stopover at the time
  stopped_pt <- st_join(pt, st %>% dplyr::select(id_yr, stop.n, year), join = st_within)
  onstop <- rbind(onstop, stopped_pt)
  
  
} #i

#rename the columns and shed useless
onstop_yr <- onstop %>% rename("id_yr" = id_yr.x) %>% dplyr::select(-id_yr.y)

nrow(onstop_yr); length(unique(onstop_yr$AID)); length(unique(onstop_yr$id_yr))

#remove duplicates, NOT SURE WHERE THESE COME FROM
onstop_yr <- onstop_yr[!duplicated(onstop_yr[, c("POSIXct", "AID")]), ]
nrow(onstop_yr); length(unique(onstop_yr$AID)); length(unique(onstop_yr$id_yr))

#first prep the stopovers, need to isolate to the AID, then for previous years, calculate the IYD between a given stopover and the stopovers in another year


onstop_yr_fin <- onstop_yr %>% filter(!is.na(stop.n))  %>% mutate(stop.n.c = paste(id_yr,"_",stop.n,sep="")) #
nrow(onstop_yr_fin); length(unique(onstop_yr_fin$AID)); length(unique(onstop_yr_fin$id_yr))

onstop_yr_fin <- onstop_yr_fin %>% group_by(stop.n.c) %>% mutate(minT = min(POSIXct), maxT = max(POSIXct), TimeStopped = (maxT - minT)/3600) %>% ungroup()
nrow(onstop_yr_fin); length(unique(onstop_yr_fin$AID)); length(unique(onstop_yr_fin$id_yr))

summary(as.numeric(onstop_yr_fin$TimeStopped))
onstop_yr_fin$TimeStopped <- as.numeric(onstop_yr_fin$TimeStopped)

setwd("C:/Users/lwilde2/Documents/Chapter2_StopoverFidelity/Chapter1_StopoverFidelity")
save(onstop_yr_fin, file = "Data/Stopover/RDH_PointsOnStop_14t22_20230425.RData")

#---------------------------------#
#---------------------------------#
## ATTN:: data drop solved!
#---------------------------------#
#---------------------------------#







# 
# t <- SO_all %>% filter(id_yr == "108_2017")
# tt <- onstop_yr %>% filter(id_yr == "108_2014")
# ttt <- data %>% filter(id_yr == "108_2014")
# z <- mapview(ttt, col.regions = "black", size = .5) + mapview(t, alpha.regions = .1) + mapview(tt, col.regions = "red", size = .5); z

#mapshot(z, url = "C:/Users/lwilde2/Downloads/FidelityAnalysis_OverlapMissesThePoint.html")





#--------------------------------------------#
# Morrison Fidelity -- Fidelity to each stop ####

#for loop version

ids <- unique(onstop_yr_fin$id_yr)

df <- data.frame()


for(i in 1:length(ids)){#length(ids)
  #i = 157
  
  # need library() here for the packages your calculations require for your calculations
  library(lubridate)
  library(sf)
  library(zoo)
  
  #define previous id's
  previous.id.1 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) - 1, "")
  previous.id.2 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -2, "")
  previous.id.3 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -3, "")
  
  #subset dataset based on current and previous 3 id_yrs
  current.id <- onstop_yr_fin %>% filter(id_yr == ids[i]) 
  
  if(previous.id.1 %in% ids){
    
    temp.1 <- onstop_yr_fin %>% filter(id_yr==previous.id.1)
    
    
    x.1 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.1 = as.numeric(min(sf::st_distance(geometry, temp.1, by_element = F)))) %>% ungroup() %>% dplyr::left_join(current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame(), by = "stop.n.c")
    
    o.1 <- current.id  %>% st_join(temp.1 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.1 <- (nrow(o.1)/nrow(current.id))*100
    
    df1 = data.frame("iyd_1" = x.1$min.IYD.1,"SpatOver_1" = rep(o.1, nrow(x.1)),"stop.n_1" = x.1$stop.n.c, "km_mark_1" = x.1$km_mark, "TimeStopped_1" = x.1$TimeStopped,"curr_Y" = rep(unique(current.id$Year), nrow(x.1)), "prev_Y_1" = rep(unique(temp.1$Year), nrow(x.1)), "AID" = rep(unique(current.id$AID) , nrow(x.1)))
    
  } else{
    df1 = data.frame("iyd_1" = NA,"SpatOver_1" = NA,"stop.n_1" = NA, "km_mark_1" = NA, "TimeStopped_1" = NA,"curr_Y" = unique(current.id$Year), "prev_Y_1" = NA, "AID" = unique(current.id$AID))
    
  }#2
  
  #second year
  if(previous.id.2 %in% ids){
    
    temp.2 <- onstop_yr_fin %>% filter(id_yr==previous.id.2)
    
    
    x.2 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.2 = as.numeric(min(sf::st_distance(geometry, temp.2, by_element = F)))) %>% ungroup() %>% dplyr::left_join(current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame(), by = "stop.n.c")
    
    o.2 <- current.id  %>% st_join(temp.2 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.2 <- (nrow(o.2) / nrow(current.id))*100
    
    df2 = data.frame("iyd_2" = x.2$min.IYD.2,"SpatOver_2" = rep(o.2, nrow(x.2)),"stop.n_2" = x.2$stop.n.c, "km_mark_2" = x.2$km_mark, "TimeStopped_2" = x.2$TimeStopped,"curr_Y" = rep(unique(current.id$Year), nrow(x.2)), "prev_Y_2" = rep(unique(temp.2$Year), nrow(x.2)), "AID" = rep(unique(current.id$AID) , nrow(x.2)))
    
  } else{
    df2 = data.frame("iyd_2" = NA,"SpatOver_2" = NA,"stop.n_2" = NA, "km_mark_2" = NA, "TimeStopped_2" = NA,"curr_Y" = unique(current.id$Year), "prev_Y_2" = NA, "AID" = unique(current.id$AID))
    
  }#2
  
  #third year
  if(previous.id.3 %in% ids){
    
    temp.3 <- onstop_yr_fin %>% filter(id_yr==previous.id.3)
    
    
    x.3 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.3 = as.numeric(min(sf::st_distance(geometry, temp.3, by_element = F)))) %>% ungroup() %>% dplyr::left_join(current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame(), by = "stop.n.c")
    
    o.3 <- current.id  %>% st_join(temp.3 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.3 <- (nrow(o.3) / nrow(current.id))*100
    
    df3 = data.frame("iyd_3" = x.3$min.IYD.3,"SpatOver_3" = rep(o.3, nrow(x.3)),"stop.n_3" = x.3$stop.n.c, "km_mark_3" = x.3$km_mark, "TimeStopped_3" = x.3$TimeStopped,"curr_Y" = rep(unique(current.id$Year), nrow(x.3)), "prev_Y_3" = rep(unique(temp.3$Year), nrow(x.3)), "AID" = rep(unique(current.id$AID) , nrow(x.3)))
    
  } else{
    df3 = data.frame("iyd_3" = NA,"SpatOver_3" = NA,"stop.n_3" = NA, "km_mark_3" = NA, "TimeStopped_3" = NA,"curr_Y" = unique(current.id$Year), "prev_Y_3" = NA, "AID" = unique(current.id$AID))
    
  }#2
  
  df.full <- cbind(df1[,c("iyd_1","SpatOver_1","stop.n_1","km_mark_1","TimeStopped_1","curr_Y","prev_Y_1","AID")], df2[,c("iyd_2","SpatOver_2","stop.n_2","km_mark_2","TimeStopped_2","curr_Y","prev_Y_2","AID")], df3[,c("iyd_3","SpatOver_3","stop.n_3","km_mark_3","TimeStopped_3","curr_Y","prev_Y_3","AID")])
  
  rm(df1, df2, df3)
  
  df <- rbind(df, df.full)
}

IYD_stop_fidelity <- df

head(IYD_stop_fidelity)

length(unique(IYD_stop_fidelity$AID))

stfid1 <-IYD_stop_fidelity %>% drop_na(prev_Y_1) %>% dplyr::select("iyd_1","SpatOver_1","stop.n_1","km_mark_1","TimeStopped_1","curr_Y","prev_Y_1","AID") %>% mutate(id_yr = paste(AID,"_",curr_Y,sep=""))
stfid2 <- IYD_stop_fidelity %>% drop_na(prev_Y_2) %>% dplyr::select("iyd_2","SpatOver_2","stop.n_2","km_mark_2","TimeStopped_2","curr_Y","prev_Y_2","AID") %>% mutate(id_yr = paste(AID,"_",curr_Y,sep=""))
stfid3 <- IYD_stop_fidelity %>% drop_na(prev_Y_3) %>% dplyr::select("iyd_3","SpatOver_3","stop.n_3","km_mark_3","TimeStopped_3","curr_Y","prev_Y_3","AID") %>% mutate(id_yr = paste(AID,"_",curr_Y,sep=""))

IYD_stop_fidelity_list <- list(stfid1, stfid2, stfid3)


#--------------------------------------------#
# Morrison Fidelity -- Mean Fidelity ####

ids <- unique(data1$id_yr)
df <- data.frame()

for(i in 1:length(ids)){#length(ids)
  #i = 31
  
  # need library() here for the packages your calculations require for your calculations
  library(lubridate)
  library(sf)
  library(zoo)
  
  #define previous id's
  previous.id.1 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) - 1, "")
  previous.id.2 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -2, "")
  previous.id.3 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -3, "")
  
  #subset dataset based on current and previous 3 id_yrs
  current.id <- onstop_yr_fin %>% filter(id_yr == ids[i]) 
  
  if(previous.id.1 %in% ids){
    
    temp.1 <- onstop_yr_fin %>% filter(id_yr==previous.id.1)
    
    
    x.1 <- current.id %>% group_by(stop.n.c) %>% mutate(min.IYD.1 = as.numeric(min(sf::st_distance(geometry, temp.1, by_element = F)))) %>% ungroup()
    mean.1 <- mean(x.1$min.IYD.1, na.rm = T)
    sd.1 <- sd(x.1$min.IYD.1, na.rm = T)
    df1 = data.frame("mean_fid_1" = mean.1, "var_fid_1" = sd.1, "curr_Y" = unique(current.id$Year), "prev_Y_1" = unique(temp.1$Year), "AID" = unique(current.id$AID))
  } else{
    df1 = data.frame("mean_fid_1" = NA, "var_fid_1" = NA, "curr_Y" = unique(current.id$Year), "prev_Y_1" = NA, "AID" = unique(current.id$AID))
    
  }#2
  
  #second year
  if(previous.id.2 %in% ids){
    
    temp.2 <- onstop_yr_fin %>% filter(id_yr==previous.id.2)
    
    x.2 <- current.id %>% group_by(stop.n.c) %>% mutate(min.IYD.2 = as.numeric(min(sf::st_distance(geometry, temp.2, by_element = F)))) %>% ungroup()
    mean.2 <- mean(x.2$min.IYD.2, na.rm = T)
    sd.2 <- sd(x.2$min.IYD.2, na.rm = T)
    df2 = data.frame("mean_fid_2" = mean.2, "var_fid_2" = sd.2, "curr_Y" = unique(current.id$Year), "prev_Y_2" = unique(temp.2$Year), "AID" = unique(current.id$AID))
  } else{
    df2 = data.frame("mean_fid_2" = NA, "var_fid_2" = NA, "curr_Y" = unique(current.id$Year), "prev_Y_2" = NA, "AID" = unique(current.id$AID))
    
  }#2
  
  #third year
  if(previous.id.3 %in% ids){
    
    temp.3 <- onstop_yr_fin %>% filter(id_yr==previous.id.3)
    
    x.3 <- current.id %>% group_by(stop.n.c) %>% mutate(min.IYD.3 = as.numeric(min(sf::st_distance(geometry, temp.3, by_element = F)))) %>% ungroup()
    mean.3 <- mean(x.3$min.IYD.3, na.rm = T)
    sd.3 <- sd(x.3$min.IYD.3, na.rm = T)
    df3 = data.frame("mean_fid_3" = mean.3, "var_fid_3" = sd.3, "curr_Y" = unique(current.id$Year), "prev_Y_3" = unique(temp.3$Year), "AID" = unique(current.id$AID))
  } else{
    df3 = data.frame("mean_fid_3" = NA, "var_fid_3" = NA, "curr_Y" = unique(current.id$Year), "prev_Y_3" = NA, "AID" = unique(current.id$AID))
    
  }#2
  
  df.full <- cbind(df1[,c("AID","curr_Y","mean_fid_1","var_fid_1","prev_Y_1")],df2[,c("mean_fid_2","var_fid_2","prev_Y_2")], df3[,c("mean_fid_3","var_fid_3","prev_Y_3")])
  
  rm(df1, df2, df3)
  
  df <- rbind(df, df.full)
}
IYD_mean_fidelity <- df





IYD_mean_fidelity <- IYD_mean_fidelity %>% mutate(id_yr = paste(AID,"_",curr_Y,sep=""))
head(IYD_mean_fidelity)

length(unique(IYD_mean_fidelity$AID)); length(unique(IYD_mean_fidelity$id_yr))

fid1 <-IYD_mean_fidelity %>% drop_na(prev_Y_1) %>% dplyr::select(AID, curr_Y, id_yr, mean_fid_1, var_fid_1, prev_Y_1)
fid2 <- IYD_mean_fidelity %>% drop_na(prev_Y_2) %>% dplyr::select(AID, curr_Y, id_yr, mean_fid_2, var_fid_2, prev_Y_2)
fid3 <- IYD_mean_fidelity %>% drop_na(prev_Y_3) %>% dplyr::select(AID, curr_Y, id_yr, mean_fid_3, var_fid_3, prev_Y_3)

IYD_mean_fidelity_list <- list(fid1, fid2, fid3)


#------------------------------------------------------#
#### fidelity to route ####

# COMMENTED FOR NOW, MOVING OVER TO DO THE STOPOVER FIDELITY

#could easily write a local function -- perhaps at which scale?

lst <- list(length(unique(data1$AID)), NA) #length(unique(data$AID))

for(j in 1:length(unique(data1$AID))){
  CalculateRouteFidelity(gps_name = data1 %>% filter(AID == unique(data1$AID)[j]) ,date_name = "POSIXct", id_name = "id_yr", buffer = 400)
  lst[[j]] <- mat
}

lst


#--------------------------#
# Save progress ####

#save(IYD_stop_fidelity_list, IYD_mean_fidelity_list, file = "C:/Users/lwilde2/Documents/Chapter2_StopoverFidelity/Chapter1_StopoverFidelity/Data/Stopover/FidelityMetrics_20230425.RData")









# NOT RUN BELOW ###


#save(lst, file = "C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/FidelityMetrics/IndividualRouteOverlap.RData")












load("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/FidelityMetrics/IndividualRouteOverlap.RData")

#keep only animals with 3 yrs data
lst_3yr <- lst[sapply(lst, nrow) >= 2 ]

df <- data.frame(matrix(nrow = length(lst_3yr), ncol = 3, data = NA))



for(i in 1:length(lst_3yr)){

  df[i,1] <- lst_3yr[[i]]$AID[1]
  df[i,2] <- median(as.numeric(lst_3yr[[i]]$Overlap))
  df[i,3] <- sd(as.numeric(lst_3yr[[i]]$Overlap))
}

names(df) <- c("AID", "median", "sd")
head(df)

ggplot(df) + geom_histogram(aes(x = median), binwidth = .02, fill = NA, color = "black")
summary(df)
bins <- sort(c(seq(0.2,0.76,by=0.02)))
RouteFidelity <- hist(df$median,breaks=bins,plot = F)
RouteFidelity <- RouteFidelity$density / 0.02

plot(RouteFidelity, type="l")

RF_med <- ggplot(df) + geom_density(aes(x = median, y = ..count..), fill = "blue", alpha = .3, color = "purple", size = 1.2, bounds = c(0,1)) + lims(x = c(0,1)) + labs(y = "Density", x = "% Fidelity") + theme_classic()
RF_sd <- ggplot(df) + geom_density(aes(x = sd, y = ..count..), fill = "green", alpha = .3, color = "orange", size = 1.2, bounds = c(0,1)) + lims(x = c(0,1)) + labs(y = "Density",x = "% Fidelity") + theme_classic()

plot1 <- RF_med + RF_sd
#ggsave(plot1, file = "C:/Users/lwilde2/Documents/AdaptiveFidelity/Outputs/RouteFidelity_20230227.png")



#parallel version, in case for loop is too slow (write to file not working, no matter!!!)
# 
# no_cores <- detectCores()-1
# ids <- unique(onstop_yr_fin$id_yr)
# 
# # Setup cluster
# clust <- makeCluster(no_cores) 
# 
# # export the objects you need in the loop
# clusterExport(clust, varlist=c("onstop_yr_fin","ids"))
# 
# #start loop
# MorrisonFidelity <- do.call(rbind, clusterApplyLB(clust, 1:length(ids), function(i){
#   #i = 5
#   # need library() here for the packages your calculations require for your calculations
#   library(lubridate)
#   library(sf)
#   library(zoo)
#   
#   #define previous id's
#   previous.id.1 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) - 1, "")
#   previous.id.2 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -2, "")
#   
#   #subset dataset based on current and previous 2 ids
#   current.id <- onstop_yr_fin %>% filter(id_yr == ids[i]) 
#   previous.1 <- onstop_yr_fin %>% filter(id_yr == previous.id.1)
#   previous.2 <- onstop_yr_fin %>% filter(id_yr == previous.id.2)
#   
#  
#   #try the subset, and block error messages
#   try(temp.1 <- onstop_yr_fin[onstop_yr_fin$id_yr==previous.id.1,], silent=TRUE)
#   try(temp.2 <- onstop_yr_fin[onstop_yr_fin$id_yr==previous.id.2,], silent=TRUE)
#   
#   #ifelse where mean, minimum distance to t-1 year will be estimated
#   if(any(class(temp.1)=="try-error")){ #if id does not exist
#     
#   } else { 
#     
#     if(previous.id.1 %in% ids){
#       
#     temp.1 <- onstop_yr_fin[onstop_yr_fin$id_yr_seas==previous.id.1,]
#       
#     
#     x.1 <- current.id %>% group_by(stop.n.c) %>% mutate(min.IYD.1 = as.numeric(min(sf::st_distance(geometry, temp.1, by_element = F)))) %>% ungroup()
#     mean.1 <- mean(x.1$min.IYD.1, na.rm = T)
#     sd.1 <- sd(x.1$min.IYD.1, na.rm = T)
#     df1 = data.frame("mean_fid_1" = mean.1, "var_fid_1" = sd.1, "curr_Y" = unique(current.id$Year), "prev_Y_1" = unique(temp.1$Year), "AID" = unique(current.id$AID))
#     } else{
#       df1 = data.frame("mean_fid_1" = NA, "var_fid_1" = NA, "curr_Y" = unique(current.id$Year), "prev_Y_1" = NA, "AID" = unique(current.id$AID))
#  
#    }#2
#   }#1
#   
#   df1
#   
#   #ifelse where mean, minimum distance to t-2 year will be estimated
#   if(any(class(temp.2)=="try-error")){ #if id does not exist
#     
#   } else { 
#     
#     if(previous.id.2 %in% ids){
#       
#       temp.2 <- onstop_yr_fin[onstop_yr_fin$id_yr_seas==previous.id.2,]
#       
#       x.2 <- current.id %>% group_by(stop.n.c) %>% mutate(min.IYD.2 = as.numeric(min(sf::st_distance(geometry, temp.2, by_element = F)))) %>% ungroup()
#       mean.2 <- mean(x.2$min.IYD.2, na.rm = T)
#       sd.2 <- sd(x.2$min.IYD.2, na.rm = T)
#       df2 = data.frame("mean_fid_2" = mean.2, "var_fid_2" = sd.2, "curr_Y" = unique(current.id$Year), "prev_Y_2" = unique(temp.2$Year), "AID" = unique(current.id$AID))
#     } else{
#       df2 = data.frame("mean_fid_2" = NA, "var_fid_2" = NA, "curr_Y" = unique(current.id$Year), "prev_Y_2" = NA, "AID" = unique(current.id$AID))
#       
#     }#2
#   }#1
#     
#     df.full <- cbind(df1[,c("AID","curr_Y","mean_fid_1","var_fid_1","prev_Y_1")],df2[,c("mean_fid_2","var_fid_2","prev_Y_2")])
#     saveRDS(df.full, paste0("C:/Users/lwilde2/Documents/AdaptiveFidelity/Data/FidelityFiles/", ids[i], ".rds"))
#   
# }))
# 
# stopCluster(clust)   # you must stop the parallelization framework
# 
# 
# 
# 





