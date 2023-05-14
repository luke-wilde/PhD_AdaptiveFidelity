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
library(MoveTools)

#load functions
source("C:/Users/lwilde2/Documents/AdaptiveFidelity/Scripts/Quantifying Fidelity/CalculatingFidelityToRoute_function.R")



#load in the datasets
# Stopovers <- st_read("Data/RD2H_HighUseStopovers_Polygons.shp")
load("C:/Users/lwilde2/Desktop/RDH Database/Processed Data/RDH_AllMigrations_Bischof_2014t2022_20230425.RData")

load("C:/Users/lwilde2/Desktop/RDH Database/Processed Data/RDH_AllStopovers_14t22_20230425.RData")

#sa <- st_read("C:/Users/lwilde2/Documents/Chapter2_StopoverFidelity/Data_out/Data/GIS/StudyArea_Polygon_TableMountains.kml") #_TableMountains if you want to include the first SO

setwd("C:/Users/lwilde2/Desktop/RDH Database/RD2H_GPSCollarData/")
DeerID <- read.table("RD2H_AllDeerIDs_17APRIL2023.csv", header=TRUE, sep=",")
#Format date and time columns
DeerID$StartDate <- as.POSIXct(DeerID$StartDate, format="%m/%d/%Y", origin = "1970-01-01")
DeerID$StopDate <- as.POSIXct(DeerID$StopDate,format="%m/%d/%Y", origin = "1970-01-01")

#prep the data

#keep only long migrations for now, will change later
DeerID_long <- DeerID %>% filter(Mgtry == "Long") %>% distinct(AID, Mgtry, .keep_all = T) %>% arrange(AID)
data <- mig.SpringMigration.sf %>% filter(AID %in% DeerID_long$AID & Mgtry == "Long")


#data <- st_as_sf(x = data, coords = c("UTM_E", "UTM_N"), crs = "+proj=utm +zone=12 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

#get sa to the same
#sa <- st_transform(sa, crs = st_crs(data)$input)

#data <- data %>% st_transform(crs = "+proj=utm +zone=12 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

#create test data
#deer108 <- data %>% filter(AID == "108") %>% mutate(id_yr = paste(AID,"_",Year,sep=""))

#data <- st_crop(data, sa)

#which animals with >1 year, keep
keep <- data %>% st_drop_geometry() %>% group_by(AID) %>% dplyr::summarise(n = n_distinct(Year)) %>% ungroup() %>% filter(n>1)

range(keep$n)

length(unique(DeerID_long$AID))
length(unique(data$AID))
length(unique(keep$AID))


data <- data %>% filter(AID %in% keep$AID)



#Calculate distance from winter range for each point using NSDs ####
#Create nsd object

# 

# dl <- Points2Lines(data = data,
#                    date_name = "POSIXct",
#                    id_name = "id_yr",
#                    byid = TRUE,
#                    no_cores = 4)
# 
# 
# c <- dl[1,] %>% group_by(id_yr) %>% st_simplify(dTolerance = 2000)
# cp <- st_line_sample(c, density = 0.001) %>% st_cast("POINT")
# mapview(dl[1,]) + mapview(c, col.regions = "red") + mapview(cp)
# 


#note here what we can do is make an average line for each AID, then use that to assign the km



datalines <- Points2Lines(data = data,
                   date_name = "POSIXct",
                   id_name = "id_yr",
                   byid = TRUE,
                   no_cores = 4)


#now distance is references per AID, so that it is comparable across years
d <- data[0,0]
for(i in 1:length(unique(data$AID))){
 # i = 15
  dl <- datalines %>% mutate(AID = str_sub(id_yr,0,3)) %>% filter(AID == unique(data$AID)[i])
  c <- dl %>% group_by(AID) %>% st_simplify(dTolerance = 1500)
  
  lst <- list(length(unique(c$id_yr)), NA)
 
  dat <- data %>% ungroup() %>% dplyr::filter(AID %in% unique(data$AID)[i])
  
  for(j in 1:length(st_geometry(c))){
  #j = 3
  
  cp <- st_as_sf(c %>% filter(id_yr == c$id_yr[j]) %>% dplyr::summarise() %>% st_line_sample(., density = 0.001) %>% st_cast("POINT"))
  cp$dist <- seq(1,length(st_geometry(cp)),1)
  st_geometry(cp) <- "geometry"
  
  lst[[j]] <- cp
  rm(cp)
 
  }#j

  all <- do.call(rbind, lst) %>% group_by(dist) %>% summarise() %>% st_centroid()
  
  dat <- dat %>% ungroup() %>% group_by(POSIXct) %>% mutate(km_mark = all$dist[sf::st_nearest_feature(geometry, all, pairwise = F)])
  
  d <- rbind(d, dat)
}#i

data <- d

nrow(data); length(unique(data$AID)); length(unique(data$id_yr))


#------------------------------------------------------#
#### Stopover Status ####

load("C:/Users/lwilde2/Desktop/RDH Database/Processed Data/RDH_AllStopovers_14t22_MergedAt5k_20230505.RData")
# 
# # writes a series of for loops to go through create a blank dataframe and iterate through each id to create a matrix, then can extract what proportion of each of the years points on its own migration stopovers are in the other ones too
# 
# #make the point layer and the SO layer within same crs
# SO_all <- st_transform(SO_all, crs = st_crs(data)$input)
# 
# SO_all_sep <- SO_all %>% mutate(AID = str_sub(id_yr, 0,3)) %>% filter(layer == 1) %>% st_cast(., "POLYGON")
# 
# #combine stops that are within 5 km of each other, takes two steps to get rid of overlap and dual counts
# # ATTN:: This code combines stops that are within 5000 n if 
# SO_all_sep_comb <- SO_all_sep[0,0]
# 
#
#
# for(i in unique(SO_all_sep$id_yr)){
#   #i = unique(SO_all_sep$id_yr)[3]
#   x <- SO_all_sep %>% filter(id_yr == i) 
#   
#   for(j in 1:nrow(x)){
#    # j = 11
#   
#     x$grp[j] <- st_is_within_distance(x$geometry[j], x, dist = 5000)
#     
#     
# }#j
#   x$grp.c <- as.numeric(factor(as.character(x$grp)))
#   
#   xx <- x %>% group_by(grp.c) %>% summarise(geometry = st_union(geometry), layer = 1, id_yr = unique(id_yr), AID = unique(AID)) %>% ungroup()
#   xx$grp1 <- NA
#   
#   rm(x)
#   
#   for(h in 1:nrow(xx)){
#      #h = 1
#     
#     xx$grp1[h] <- st_is_within_distance(xx$geometry[h], xx, dist = 5000)
#     
#   }#h
# 
#   xx$grp1.c <- as.numeric(factor(as.character(xx$grp1)))
#   xx <- xx %>% dplyr::select(-c(grp1,grp.c))
#   xx <- xx %>% group_by(grp1.c) %>% summarise(geometry = st_union(geometry), layer = 1, id_yr = unique(id_yr), AID = unique(AID)) %>% ungroup()
#   xx <- xx %>% dplyr::select(layer, id_yr, AID, grp1.c, geometry)
#   SO_all_sep_comb <- rbind(SO_all_sep_comb, xx)
#   rm(xx)
#   }#i
# 
# SO_all_sep_comb



#inserted it here V
#SO_all_sep <- SO_all_sep_comb %>% group_by(id_yr) %>% mutate(stop.n = row_number(), year = as.numeric(str_sub(id_yr,5,8))) %>% ungroup()

length(unique(SO_all_sep$AID)); length(unique(SO_all_sep$id_yr))
length(unique(data$AID)); length(unique(data$id_yr))


table(unique(data$id_yr) %in% unique(SO_all_sep$id_yr))

check <- SO_all_sep[which(is.na(SO_all_sep$stop.n)),]
length(unique(check$AID)); length(unique(check$id_yr))

#get rid of animals missing SO info
missingSO <- data %>% filter(!id_yr %in% SO_all_sep$id_yr) %>% arrange(id_yr) %>% st_drop_geometry() %>% dplyr::select(id_yr)
unique(missingSO$id_yr)





data1 <- data %>% filter(!id_yr %in% missingSO$id_yr) %>% arrange(id_yr)


data1 <- data1[!duplicated(data1[, c("POSIXct", "AID")]), ]
nrow(data1); length(unique(data1$AID)); length(unique(data1$id_yr))


SO_all_sep

# Step 1: Extract points on stops within same year, for loop

onstop <- data.frame()

ids <- unique(data1$id_yr)

for(i in 1:length(ids)){ #
 # i = 201
  #
  pt <- data1 %>% dplyr::filter(id_yr == ids[i]) #[which(data1$id_yr == ids[i]),]#no idea why this only works in base now...
  st <- SO_all_sep %>% dplyr::filter(id_yr == ids[i])
  
  #can keep this or not...
  st <- st %>% group_by(stop.n) %>% dplyr::summarise(id_yr = unique(id_yr), year = unique(year)) %>% st_convex_hull()
  st_crs(pt) == st_crs(st)
  
  #extract points that are inside a stopover at the time
  pt <- st_join(pt, st)
  stopped_pt <- pt %>% drop_na(stop.n)
  onstop <- rbind(onstop, stopped_pt)
  
  
} #i

View(onstop)

#rename the columns and shed useless
onstop_yr <- onstop %>% rename("id_yr" = id_yr.x) %>% dplyr::select(-id_yr.y)


nrow(onstop_yr); length(unique(onstop_yr$AID)); length(unique(onstop_yr$id_yr))

#remove duplicates, NOT SURE WHERE THESE COME FROM
onstop_yr <- onstop_yr[!duplicated(onstop_yr[, c("POSIXct", "AID")]), ]
nrow(onstop_yr); length(unique(onstop_yr$AID)); length(unique(onstop_yr$id_yr))

#first prep the stopovers, need to isolate to the AID, then for previous years, calculate the IYD between a given stopover and the stopovers in another year


onstop_yr_fin <- onstop_yr %>% filter(!is.na(stop.n))  %>% mutate(stop.n.c = paste(id_yr,"_",stop.n,sep="")) #
nrow(onstop_yr_fin); length(unique(onstop_yr_fin$AID)); length(unique(onstop_yr_fin$id_yr))




onstop_yr_fin <- onstop_yr_fin %>% group_by(stop.n.c) %>% mutate(minT = min(POSIXct), maxT = max(POSIXct), TimeStopped = (maxT - minT)/3600, km_mark = mean(km_mark)) %>% ungroup()
nrow(onstop_yr_fin); length(unique(onstop_yr_fin$AID)); length(unique(onstop_yr_fin$id_yr))

summary(as.numeric(onstop_yr_fin$TimeStopped))
summary(as.numeric(onstop_yr_fin$km_mark))

onstop_yr_fin$TimeStopped <- as.numeric(onstop_yr_fin$TimeStopped)


#repeat the below for each one of these
onhigh_yr_fin <- onstop_yr_fin %>% filter(TimeStopped >= 72)
onstop_yr_fin <- onstop_yr_fin %>% filter(TimeStopped >= 24)

setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
save(onstop_yr_fin, data1, onhigh_yr_fin, file = "Data_out/Data/Stopover/RDH_StopoverBundle_14t22_20230507.RData") #crap accidentally overwrote on 20230503

#---------------------------------#
#---------------------------------#
## ATTN:: data drop solved!
#---------------------------------#
#---------------------------------#


setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
load("Data_out/Data/Stopover/RDH_StopoverBundle_14t22_20230507.RData")



# 
# t <- SO_all %>% filter(id_yr == "108_2017")
# tt <- onstop_yr %>% filter(id_yr == "108_2014")
# ttt <- data %>% filter(id_yr == "108_2014")
# z <- mapview(ttt, col.regions = "black", size = .5) + mapview(t, alpha.regions = .1) + mapview(tt, col.regions = "red", size = .5); z

#mapshot(z, url = "C:/Users/lwilde2/Downloads/FidelityAnalysis_OverlapMissesThePoint.html")

# t <- onstop_yr_fin[which(onstop_yr_fin$id_yr == "108_2020"),] %>% group_by(stop.n.c) %>% dplyr::summarise() %>%
#   st_convex_hull()
# mapview(t)


#correct for differences in starting position? Is that the issue with the km_mark?

##ATTN -- ONLY WORKS FOR LONG AND SOME MEDIUM!!

SuperiorExit <- st_as_sf(x = data.frame(x = 672301.29, y = 4615297.52), coords = c("x", "y"), crs = st_crs(onstop_yr_fin)$input)
#mapview(SuperiorExit)
# 
# set <- onstop_yr_fin %>% mutate(set = as.numeric(st_distance(., SuperiorExit))) %>% group_by(id_yr) %>% dplyr::summarise(sup = min(km_mark - (set/1000)))
# 
# x <- onstop_yr_fin %>% left_join(set %>% st_drop_geometry(), by = "id_yr") %>% mutate(., km_mark = km_mark + sup)
# 
# x[,c("km_mark", "sup")]
# 
# onstop_yr_fin <-  x %>% dplyr::select(-c(sup))


#STILL NEED TO correct for the km differences -- tmw


#--------------------------------------------#
# Morrison Fidelity -- Fidelity to each stop ####

#for loop version

ids <- unique(onstop_yr_fin$id_yr)

df <- data.frame()


for(i in 1:length(ids)){#length(ids)
  #i = 6
  
  # need library() here for the packages your calculations require for your calculations
  library(lubridate)
  library(sf)
  library(zoo)
  
  #define previous id's
  previous.id.1 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) - 1, "")
  previous.id.2 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -2, "")
  previous.id.3 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -3, "")
  previous.id.4 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) - 4, "")
  previous.id.5 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -5, "")
  previous.id.6 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -6, "")
  previous.id.7 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -7, "")
  
  #subset dataset based on current and previous 3 id_yrs
  current.id <- onstop_yr_fin %>% filter(id_yr == ids[i]) 
  
  
  
  if(previous.id.1 %in% ids){
    
    temp.1 <- onstop_yr_fin %>% filter(id_yr==previous.id.1)
  
    #temp.1 <- temp.1 %>% group_by(stop.n.c) %>% mutate(km_mark = mean(km_mark))
    
    # temp.1 <- temp.1 %>% mutate(km_mark_corr = current.id$km_mark[sf::st_nearest_feature(geometry, current.id, pairwise = F)])
    
    x.1 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.1 = as.numeric(min(sf::st_distance(temp.1, geometry, by_element = F)))) %>% ungroup()
    xx.1 <- current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame()
    
    x.1 <- x.1 %>% left_join(xx.1, by = "stop.n.c")
    
   s.1 <- current.id %>% group_by(stop.n.c) %>% dplyr::reframe(stop_prev = temp.1$stop.n.c[sf::st_nearest_feature(geometry, temp.1, pairwise = F)], km_prev = temp.1$km_mark[sf::st_nearest_feature(geometry, temp.1, pairwise = F)]) %>% distinct(stop.n.c, stop_prev, km_prev) %>%  as.data.frame()
   
   
   
   x.1 <- x.1 %>% dplyr::left_join(s.1, by = c("stop.n.c"))
    
    o.1 <- current.id  %>% st_join(temp.1 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.1 <- (nrow(o.1)/nrow(current.id))*100
    
    df1 = data.frame("iyd_1" = x.1$min.IYD.1, "stop_prev_1" = x.1$stop_prev, "km_prev_1" = x.1$km_prev, "SpatOver_1" = rep(o.1, nrow(x.1)),"stop.n_1" = x.1$stop.n.c, "km_mark_1" = x.1$km_mark, "TimeStopped_1" = x.1$TimeStopped,"curr_Y" = rep(unique(current.id$Year), nrow(x.1)), "prev_Y_1" = rep(unique(temp.1$Year), nrow(x.1)), "AID" = rep(unique(current.id$AID) , nrow(x.1)))
    
    # t1 <- temp.1 %>% group_by(stop.n.c) %>% dplyr::summarise(km_mark = km_mark[1]) %>%
    #   st_convex_hull() #%>% dplyr::select(stop.n.c, km_mark)
    # t2 <- current.id
    # 
    # mapview(t1, col.regions = "red") + mapview(t2)
    
  } else{
    df1 = data.frame("iyd_1" = NA,"SpatOver_1" = NA,"stop.n_1" = NA, "km_mark_1" = NA, "TimeStopped_1" = NA,"curr_Y" = unique(current.id$Year), "prev_Y_1" = NA,"stop_prev_1" = NA, "km_prev_1" = NA, "AID" = unique(current.id$AID))
    
  }#2
  
  #second year
  if(previous.id.2 %in% ids){
    
    
    temp.2 <- onstop_yr_fin %>% filter(id_yr==previous.id.2)
    
    #temp.2 <- temp.2 %>% group_by(stop.n.c) %>% mutate(km_mark = mean(km_mark))
    
    # temp.2 <- temp.2 %>% mutate(km_mark_corr = current.id$km_mark[sf::st_nearest_feature(geometry, current.id, pairwise = F)])
    
    x.2 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.2 = as.numeric(min(sf::st_distance(temp.2, geometry, by_element = F)))) %>% ungroup()
    xx.2 <- current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame()
    
    x.2 <- x.2 %>% left_join(xx.2, by = "stop.n.c")
    
    s.2 <- current.id %>% group_by(stop.n.c) %>% dplyr::reframe(stop_prev = temp.2$stop.n.c[sf::st_nearest_feature(geometry, temp.2, pairwise = F)], km_prev = temp.2$km_mark[sf::st_nearest_feature(geometry, temp.2, pairwise = F)]) %>% distinct(stop.n.c, stop_prev, km_prev) %>%  as.data.frame()
    
    
    
    x.2 <- x.2 %>% dplyr::left_join(s.2, by = c("stop.n.c"))
    
    o.2 <- current.id  %>% st_join(temp.2 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.2 <- (nrow(o.2)/nrow(current.id))*100
    
    df2 = data.frame("iyd_2" = x.2$min.IYD.2, "stop_prev_2" = x.2$stop_prev, "km_prev_2" = x.2$km_prev, "SpatOver_2" = rep(o.2, nrow(x.2)),"stop.n_2" = x.2$stop.n.c, "km_mark_2" = x.2$km_mark, "TimeStopped_2" = x.2$TimeStopped,"curr_Y.2" = rep(unique(current.id$Year), nrow(x.2)), "prev_Y_2" = rep(unique(temp.2$Year), nrow(x.2)), "AID.2" = rep(unique(current.id$AID) , nrow(x.2)))
    
    
  } else{
    df2 = data.frame("iyd_2" = NA,"SpatOver_2" = NA,"stop.n_2" = NA, "km_mark_2" = NA, "TimeStopped_2" = NA,"curr_Y.2" = unique(current.id$Year), "prev_Y_2" = NA,"stop_prev_2" = NA, "km_prev_2" = NA, "AID.2" = unique(current.id$AID))
    
  }#2
  
  #third year
  if(previous.id.3 %in% ids){
    
    
    temp.3 <- onstop_yr_fin %>% filter(id_yr==previous.id.3)
    
    #temp.3 <- temp.3 %>% group_by(stop.n.c) %>% mutate(km_mark = mean(km_mark))
    
    # temp.3 <- temp.3 %>% mutate(km_mark_corr = current.id$km_mark[sf::st_nearest_feature(geometry, current.id, pairwise = F)])
    
    x.3 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.3 = as.numeric(min(sf::st_distance(temp.3, geometry, by_element = F)))) %>% ungroup()
    xx.3 <- current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame()
    
    x.3 <- x.3 %>% left_join(xx.3, by = "stop.n.c")
    
    s.3 <- current.id %>% group_by(stop.n.c) %>% dplyr::reframe(stop_prev = temp.3$stop.n.c[sf::st_nearest_feature(geometry, temp.3, pairwise = F)], km_prev = temp.3$km_mark[sf::st_nearest_feature(geometry, temp.3, pairwise = F)]) %>% distinct(stop.n.c, stop_prev, km_prev) %>%  as.data.frame()
    
    
    
    x.3 <- x.3 %>% dplyr::left_join(s.3, by = c("stop.n.c"))
    
    o.3 <- current.id  %>% st_join(temp.3 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.3 <- (nrow(o.3)/nrow(current.id))*100
    
    df3 = data.frame("iyd_3" = x.3$min.IYD.3, "stop_prev_3" = x.3$stop_prev, "km_prev_3" = x.3$km_prev, "SpatOver_3" = rep(o.3, nrow(x.3)),"stop.n_3" = x.3$stop.n.c, "km_mark_3" = x.3$km_mark, "TimeStopped_3" = x.3$TimeStopped,"curr_Y.3" = rep(unique(current.id$Year), nrow(x.3)), "prev_Y_3" = rep(unique(temp.3$Year), nrow(x.3)), "AID.3" = rep(unique(current.id$AID) , nrow(x.3)))
    
    
  } else{
    df3 = data.frame("iyd_3" = NA,"SpatOver_3" = NA,"stop.n_3" = NA, "km_mark_3" = NA, "TimeStopped_3" = NA,"curr_Y.3" = unique(current.id$Year), "prev_Y_3" = NA,"stop_prev_3" = NA, "km_prev_3" = NA, "AID.3" = unique(current.id$AID))
    
  }#2
  
  
  if(previous.id.4 %in% ids){
    
    temp.4 <- onstop_yr_fin %>% filter(id_yr==previous.id.4)
    
    #temp.4 <- temp.4 %>% group_by(stop.n.c) %>% mutate(km_mark = mean(km_mark))
    
    # temp.4 <- temp.4 %>% mutate(km_mark_corr = current.id$km_mark[sf::st_nearest_feature(geometry, current.id, pairwise = F)])
    
    x.4 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.4 = as.numeric(min(sf::st_distance(temp.4, geometry, by_element = F)))) %>% ungroup()
    xx.4 <- current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame()
    
    x.4 <- x.4 %>% left_join(xx.4, by = "stop.n.c")
    
    s.4 <- current.id %>% group_by(stop.n.c) %>% dplyr::reframe(stop_prev = temp.4$stop.n.c[sf::st_nearest_feature(geometry, temp.4, pairwise = F)], km_prev = temp.4$km_mark[sf::st_nearest_feature(geometry, temp.4, pairwise = F)]) %>% distinct(stop.n.c, stop_prev, km_prev) %>%  as.data.frame()
    
    
    
    x.4 <- x.4 %>% dplyr::left_join(s.4, by = c("stop.n.c"))
    
    o.4 <- current.id  %>% st_join(temp.4 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.4 <- (nrow(o.4)/nrow(current.id))*100
    
    df4 = data.frame("iyd_4" = x.4$min.IYD.4, "stop_prev_4" = x.4$stop_prev, "km_prev_4" = x.4$km_prev, "SpatOver_4" = rep(o.4, nrow(x.4)),"stop.n_4" = x.4$stop.n.c, "km_mark_4" = x.4$km_mark, "TimeStopped_4" = x.4$TimeStopped,"curr_Y.4" = rep(unique(current.id$Year), nrow(x.4)), "prev_Y_4" = rep(unique(temp.4$Year), nrow(x.4)), "AID.4" = rep(unique(current.id$AID) , nrow(x.4)))
    
    # t1 <- temp.4 %>% group_by(stop.n.c) %>% dplyr::summarise(km_mark = km_mark[1]) %>%
    #   st_convex_hull() #%>% dplyr::select(stop.n.c, km_mark)
    # t2 <- current.id
    # 
    # mapview(t1, col.regions = "red") + mapview(t2)
    
  } else{
    df4 = data.frame("iyd_4" = NA,"SpatOver_4" = NA,"stop.n_4" = NA, "km_mark_4" = NA, "TimeStopped_4" = NA,"curr_Y.4" = unique(current.id$Year), "prev_Y_4" = NA,"stop_prev_4" = NA, "km_prev_4" = NA, "AID.4" = unique(current.id$AID))
    
  }#2
  
  
  if(previous.id.5 %in% ids){
    
    temp.5 <- onstop_yr_fin %>% filter(id_yr==previous.id.5)
    
    #temp.5 <- temp.5 %>% group_by(stop.n.c) %>% mutate(km_mark = mean(km_mark))
    
    # temp.5 <- temp.5 %>% mutate(km_mark_corr = current.id$km_mark[sf::st_nearest_feature(geometry, current.id, pairwise = F)])
    
    x.5 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.5 = as.numeric(min(sf::st_distance(temp.5, geometry, by_element = F)))) %>% ungroup()
    xx.5 <- current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame()
    
    x.5 <- x.5 %>% left_join(xx.5, by = "stop.n.c")
    
    s.5 <- current.id %>% group_by(stop.n.c) %>% dplyr::reframe(stop_prev = temp.5$stop.n.c[sf::st_nearest_feature(geometry, temp.5, pairwise = F)], km_prev = temp.5$km_mark[sf::st_nearest_feature(geometry, temp.5, pairwise = F)]) %>% distinct(stop.n.c, stop_prev, km_prev) %>%  as.data.frame()
    
    
    
    x.5 <- x.5 %>% dplyr::left_join(s.5, by = c("stop.n.c"))
    
    o.5 <- current.id  %>% st_join(temp.5 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.5 <- (nrow(o.5)/nrow(current.id))*100
    
    df5 = data.frame("iyd_5" = x.5$min.IYD.5, "stop_prev_5" = x.5$stop_prev, "km_prev_5" = x.5$km_prev, "SpatOver_5" = rep(o.5, nrow(x.5)),"stop.n_5" = x.5$stop.n.c, "km_mark_5" = x.5$km_mark, "TimeStopped_5" = x.5$TimeStopped,"curr_Y.5" = rep(unique(current.id$Year), nrow(x.5)), "prev_Y_5" = rep(unique(temp.5$Year), nrow(x.5)), "AID.5" = rep(unique(current.id$AID) , nrow(x.5)))
    
    # t1 <- temp.5 %>% group_by(stop.n.c) %>% dplyr::summarise(km_mark = km_mark[1]) %>%
    #   st_convex_hull() #%>% dplyr::select(stop.n.c, km_mark)
    # t2 <- current.id
    # 
    # mapview(t1, col.regions = "red") + mapview(t2)
    
  } else{
    df5 = data.frame("iyd_5" = NA,"SpatOver_5" = NA,"stop.n_5" = NA, "km_mark_5" = NA, "TimeStopped_5" = NA,"curr_Y.5" = unique(current.id$Year), "prev_Y_5" = NA,"stop_prev_5" = NA, "km_prev_5" = NA, "AID.5" = unique(current.id$AID))
    
  }#2
  
  
  if(previous.id.6 %in% ids){
    
    temp.6 <- onstop_yr_fin %>% filter(id_yr==previous.id.6)
    
    #temp.6 <- temp.6 %>% group_by(stop.n.c) %>% mutate(km_mark = mean(km_mark))
    
    # temp.6 <- temp.6 %>% mutate(km_mark_corr = current.id$km_mark[sf::st_nearest_feature(geometry, current.id, pairwise = F)])
    
    x.6 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.6 = as.numeric(min(sf::st_distance(temp.6, geometry, by_element = F)))) %>% ungroup()
    xx.6 <- current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame()
    
    x.6 <- x.6 %>% left_join(xx.6, by = "stop.n.c")
    
    s.6 <- current.id %>% group_by(stop.n.c) %>% dplyr::reframe(stop_prev = temp.6$stop.n.c[sf::st_nearest_feature(geometry, temp.6, pairwise = F)], km_prev = temp.6$km_mark[sf::st_nearest_feature(geometry, temp.6, pairwise = F)]) %>% distinct(stop.n.c, stop_prev, km_prev) %>%  as.data.frame()
    
    
    
    x.6 <- x.6 %>% dplyr::left_join(s.6, by = c("stop.n.c"))
    
    o.6 <- current.id  %>% st_join(temp.6 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.6 <- (nrow(o.6)/nrow(current.id))*100
    
    df6 = data.frame("iyd_6" = x.6$min.IYD.6, "stop_prev_6" = x.6$stop_prev, "km_prev_6" = x.6$km_prev, "SpatOver_6" = rep(o.6, nrow(x.6)),"stop.n_6" = x.6$stop.n.c, "km_mark_6" = x.6$km_mark, "TimeStopped_6" = x.6$TimeStopped,"curr_Y.6" = rep(unique(current.id$Year), nrow(x.6)), "prev_Y_6" = rep(unique(temp.6$Year), nrow(x.6)), "AID.6" = rep(unique(current.id$AID) , nrow(x.6)))
    
    # t1 <- temp.6 %>% group_by(stop.n.c) %>% dplyr::summarise(km_mark = km_mark[1]) %>%
    #   st_convex_hull() #%>% dplyr::select(stop.n.c, km_mark)
    # t2 <- current.id
    # 
    # mapview(t1, col.regions = "red") + mapview(t2)
    
  } else{
    df6 = data.frame("iyd_6" = NA,"SpatOver_6" = NA,"stop.n_6" = NA, "km_mark_6" = NA, "TimeStopped_6" = NA,"curr_Y.6" = unique(current.id$Year), "prev_Y_6" = NA,"stop_prev_6" = NA, "km_prev_6" = NA, "AID.6" = unique(current.id$AID))
    
  }#2
  
  
  if(previous.id.7 %in% ids){
    
    temp.7 <- onstop_yr_fin %>% filter(id_yr==previous.id.7)
    
    #temp.7 <- temp.7 %>% group_by(stop.n.c) %>% mutate(km_mark = mean(km_mark))
    
    # temp.7 <- temp.7 %>% mutate(km_mark_corr = current.id$km_mark[sf::st_nearest_feature(geometry, current.id, pairwise = F)])
    
    x.7 <- current.id %>% group_by(stop.n.c) %>% summarise(min.IYD.7 = as.numeric(min(sf::st_distance(temp.7, geometry, by_element = F)))) %>% ungroup()
    xx.7 <- current.id %>% st_drop_geometry() %>% group_by(stop.n.c) %>% summarise(km_mark = mean(km_mark), TimeStopped = mean(TimeStopped)) %>% dplyr::select(TimeStopped, km_mark, stop.n.c) %>% as.data.frame()
    
    x.7 <- x.7 %>% left_join(xx.7, by = "stop.n.c")
    
    s.7 <- current.id %>% group_by(stop.n.c) %>% dplyr::reframe(stop_prev = temp.7$stop.n.c[sf::st_nearest_feature(geometry, temp.7, pairwise = F)], km_prev = temp.7$km_mark[sf::st_nearest_feature(geometry, temp.7, pairwise = F)]) %>% distinct(stop.n.c, stop_prev, km_prev) %>%  as.data.frame()
    
    
    
    x.7 <- x.7 %>% dplyr::left_join(s.7, by = c("stop.n.c"))
    
    o.7 <- current.id  %>% st_join(temp.7 %>% dplyr::group_by(stop.n.c) %>% dplyr::summarise(n = n()) %>% filter(n > 3) %>% st_cast("POLYGON"), join = st_within) %>% drop_na(stop.n.c.y)
    o.7 <- (nrow(o.7)/nrow(current.id))*100
    
    df7 = data.frame("iyd_7" = x.7$min.IYD.7, "stop_prev_7" = x.7$stop_prev, "km_prev_7" = x.7$km_prev, "SpatOver_7" = rep(o.7, nrow(x.7)),"stop.n_7" = x.7$stop.n.c, "km_mark_7" = x.7$km_mark, "TimeStopped_7" = x.7$TimeStopped,"curr_Y.7" = rep(unique(current.id$Year), nrow(x.7)), "prev_Y_7" = rep(unique(temp.7$Year), nrow(x.7)), "AID.7" = rep(unique(current.id$AID) , nrow(x.7)))
    
    # t1 <- temp.7 %>% group_by(stop.n.c) %>% dplyr::summarise(km_mark = km_mark[1]) %>%
    #   st_convex_hull() #%>% dplyr::select(stop.n.c, km_mark)
    # t2 <- current.id
    # 
    # mapview(t1, col.regions = "red") + mapview(t2)
    
  } else{
    df7 = data.frame("iyd_7" = NA,"SpatOver_7" = NA,"stop.n_7" = NA, "km_mark_7" = NA, "TimeStopped_7" = NA,"curr_Y.7" = unique(current.id$Year), "prev_Y_7" = NA,"stop_prev_7" = NA, "km_prev_7" = NA, "AID.7" = unique(current.id$AID))
    
  }#2
  
  df.full <- cbind(df1[,c("iyd_1","SpatOver_1","stop.n_1","km_mark_1","TimeStopped_1","curr_Y","prev_Y_1","stop_prev_1", "km_prev_1", "AID")], df2[,c("iyd_2","SpatOver_2","stop.n_2","km_mark_2","TimeStopped_2","curr_Y.2","prev_Y_2","stop_prev_2", "km_prev_2","AID.2")], df3[,c("iyd_3","SpatOver_3","stop.n_3","km_mark_3","TimeStopped_3","curr_Y.3","prev_Y_3","stop_prev_3", "km_prev_3","AID.3")], df4[,c("iyd_4","SpatOver_4","stop.n_4","km_mark_4","TimeStopped_4","curr_Y.4","prev_Y_4","stop_prev_4", "km_prev_4", "AID.4")],df5[,c("iyd_5","SpatOver_5","stop.n_5","km_mark_5","TimeStopped_5","curr_Y.5","prev_Y_5","stop_prev_5", "km_prev_5", "AID.5")],df6[,c("iyd_6","SpatOver_6","stop.n_6","km_mark_6","TimeStopped_6","curr_Y.6","prev_Y_6","stop_prev_6", "km_prev_6", "AID.6")],df7[,c("iyd_7","SpatOver_7","stop.n_7","km_mark_7","TimeStopped_7","curr_Y.7","prev_Y_7","stop_prev_7", "km_prev_7", "AID.7")])
  
  rm(df1, df2, df3,df4, df5, df6,df7)
  
  df <- rbind(df, df.full)
}

IYD_stop_fidelity <- df

#names(IYD_stop_fidelity)[c(6,26,16,10,20,30)] <- c("curr_Y", "curr_Y.2", "curr_Y.3","AID", "AID.2", "AID.3")


head(IYD_stop_fidelity)






length(unique(IYD_stop_fidelity$AID))

stfid1 <-IYD_stop_fidelity %>% drop_na(prev_Y_1) %>% dplyr::select("iyd_1","SpatOver_1","stop.n_1","km_mark_1","km_prev_1","TimeStopped_1","stop_prev_1","curr_Y","prev_Y_1","AID") %>% mutate(id_yr = paste(AID,"_",curr_Y,sep=""))
stfid2 <- IYD_stop_fidelity %>% drop_na(prev_Y_2) %>% dplyr::select("iyd_2","SpatOver_2","stop.n_2","km_mark_2","km_prev_2","TimeStopped_2","stop_prev_2","curr_Y.2","prev_Y_2","AID.2") %>% mutate(id_yr = paste(AID.2,"_",curr_Y.2,sep=""))
stfid3 <- IYD_stop_fidelity %>% drop_na(prev_Y_3) %>% dplyr::select("iyd_3","SpatOver_3","stop.n_3","km_mark_3","km_prev_3","TimeStopped_3","stop_prev_3","curr_Y.3","prev_Y_3","AID.3") %>% mutate(id_yr = paste(AID.3,"_",curr_Y.3,sep=""))
stfid4 <- IYD_stop_fidelity %>% drop_na(prev_Y_4) %>% dplyr::select("iyd_4","SpatOver_4","stop.n_4","km_mark_4","km_prev_4","TimeStopped_4","stop_prev_4","curr_Y.4","prev_Y_4","AID.4") %>% mutate(id_yr = paste(AID.4,"_",curr_Y.4,sep=""))
stfid5 <- IYD_stop_fidelity %>% drop_na(prev_Y_5) %>% dplyr::select("iyd_5","SpatOver_5","stop.n_5","km_mark_5","km_prev_5","TimeStopped_5","stop_prev_5","curr_Y.5","prev_Y_5","AID.5") %>% mutate(id_yr = paste(AID.5,"_",curr_Y.5,sep=""))
stfid6 <- IYD_stop_fidelity %>% drop_na(prev_Y_6) %>% dplyr::select("iyd_6","SpatOver_6","stop.n_6","km_mark_6","km_prev_6","TimeStopped_6","stop_prev_6","curr_Y.6","prev_Y_6","AID.6") %>% mutate(id_yr = paste(AID.6,"_",curr_Y.6,sep=""))
stfid7 <- IYD_stop_fidelity %>% drop_na(prev_Y_7) %>% dplyr::select("iyd_7","SpatOver_7","stop.n_7","km_mark_7","km_prev_7","TimeStopped_7","stop_prev_7","curr_Y.7","prev_Y_7","AID.7") %>% mutate(id_yr = paste(AID.7,"_",curr_Y.7,sep=""))

IYD_stop_fidelity_list <- list(stfid1, stfid2, stfid3,stfid4,stfid5,stfid6,stfid7)

num <- sapply(list(sapply(list(stfid1$id_yr,stfid2$id_yr,stfid3$id_yr,stfid4$id_yr,stfid5$id_yr,stfid6$id_yr,stfid7$id_yr),unique)),as.vector)




save(IYD_stop_fidelity_list, num, file = "Data_out/Data/Stopover/FidelityMetrics_20230507.RData")

#--------------------------------------------#
# Morrison Fidelity -- Mean Fidelity ####

## NOT RUN, NEED TO UPDATE TO BE CORRECT ##
# 
# ids <- unique(data1$id_yr)
# df <- data.frame()
# 
# for(i in 1:length(ids)){#length(ids)
#   #i = 31
#   
#   tryCatch({
#   # need library() here for the packages your calculations require for your calculations
#   library(lubridate)
#   library(sf)
#   library(zoo)
#   
#   #define previous id's
#   previous.id.1 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) - 1, "")
#   previous.id.2 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -2, "")
#   previous.id.3 <- paste0(sapply(strsplit(ids[i], "_"), "[", 1), "_", as.numeric(sapply(strsplit(ids[i], "_"), "[", 2)) -3, "")
#   
#   #subset dataset based on current and previous 3 id_yrs
#   current.id <- onhigh_yr_fin %>% filter(id_yr == ids[i]) 
#   
#   if(previous.id.1 %in% ids){
#     
#     temp.1 <- onhigh_yr_fin %>% filter(id_yr==previous.id.1)
#     
#     
#     x.1 <- current.id %>% group_by(stop.n.c) %>% mutate(min.IYD.1 = as.numeric(min(sf::st_distance(geometry, temp.1, by_element = F)))) %>% ungroup()
#     mean.1 <- mean(x.1$min.IYD.1, na.rm = T)
#     sd.1 <- sd(x.1$min.IYD.1, na.rm = T)
#     df1 = data.frame("mean_fid_1" = mean.1, "var_fid_1" = sd.1, "curr_Y" = unique(current.id$Year), "prev_Y_1" = unique(temp.1$Year), "AID" = unique(current.id$AID))
#   } else{
#     df1 = data.frame("mean_fid_1" = NA, "var_fid_1" = NA, "curr_Y" = unique(current.id$Year), "prev_Y_1" = NA, "AID" = unique(current.id$AID))
#     
#   }#2
#   
#   #second year
#   if(previous.id.2 %in% ids){
#     
#     temp.2 <- onhigh_yr_fin %>% filter(id_yr==previous.id.2)
#     
#     x.2 <- current.id %>% group_by(stop.n.c) %>% mutate(min.IYD.2 = as.numeric(min(sf::st_distance(geometry, temp.2, by_element = F)))) %>% ungroup()
#     mean.2 <- mean(x.2$min.IYD.2, na.rm = T)
#     sd.2 <- sd(x.2$min.IYD.2, na.rm = T)
#     df2 = data.frame("mean_fid_2" = mean.2, "var_fid_2" = sd.2, "curr_Y" = unique(current.id$Year), "prev_Y_2" = unique(temp.2$Year), "AID" = unique(current.id$AID))
#   } else{
#     df2 = data.frame("mean_fid_2" = NA, "var_fid_2" = NA, "curr_Y" = unique(current.id$Year), "prev_Y_2" = NA, "AID" = unique(current.id$AID))
#     
#   }#2
#   
#   #third year
#   if(previous.id.3 %in% ids){
#     
#     temp.3 <- onhigh_yr_fin %>% filter(id_yr==previous.id.3)
#     
#     x.3 <- current.id %>% group_by(stop.n.c) %>% mutate(min.IYD.3 = as.numeric(min(sf::st_distance(geometry, temp.3, by_element = F)))) %>% ungroup()
#     mean.3 <- mean(x.3$min.IYD.3, na.rm = T)
#     sd.3 <- sd(x.3$min.IYD.3, na.rm = T)
#     df3 = data.frame("mean_fid_3" = mean.3, "var_fid_3" = sd.3, "curr_Y" = unique(current.id$Year), "prev_Y_3" = unique(temp.3$Year), "AID" = unique(current.id$AID))
#   } else{
#     df3 = data.frame("mean_fid_3" = NA, "var_fid_3" = NA, "curr_Y" = unique(current.id$Year), "prev_Y_3" = NA, "AID" = unique(current.id$AID))
#     
#   }#2
#   
#   df.full <- cbind(df1[,c("AID","curr_Y","mean_fid_1","var_fid_1","prev_Y_1")],df2[,c("mean_fid_2","var_fid_2","prev_Y_2")], df3[,c("mean_fid_3","var_fid_3","prev_Y_3")])
#   
#   rm(df1, df2, df3)
#   
#   df <- rbind(df, df.full)
#   }, error=function(e){})
# }
# IYD_HighUse_mean_fidelity <- df
# 
# 
# 
# 
# 
# IYD_HighUse_mean_fidelity <- IYD_HighUse_mean_fidelity %>% mutate(id_yr = paste(AID,"_",curr_Y,sep=""))
# head(IYD_mean_fidelity)
# 
# length(unique(IYD_mean_fidelity$AID)); length(unique(IYD_mean_fidelity$id_yr))
# 
# fid1 <-IYD_mean_fidelity %>% drop_na(prev_Y_1) %>% dplyr::select(AID, curr_Y, id_yr, mean_fid_1, var_fid_1, prev_Y_1)
# fid2 <- IYD_mean_fidelity %>% drop_na(prev_Y_2) %>% dplyr::select(AID, curr_Y, id_yr, mean_fid_2, var_fid_2, prev_Y_2)
# fid3 <- IYD_mean_fidelity %>% drop_na(prev_Y_3) %>% dplyr::select(AID, curr_Y, id_yr, mean_fid_3, var_fid_3, prev_Y_3)
# 
# IYD_mean_fidelity_list <- list(fid1, fid2, fid3)


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

save(IYD_stop_fidelity_list, file = "C:/Users/lwilde2/Documents/Chapter2_StopoverFidelity/Data_out/Data/Stopover/FidelityMetrics_20230506.RData")

#, IYD_HighUse_mean_fidelity_list, IYD_mean_fidelity_list







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


#-------------------------------#
# how much is IYD just keep moving ####

library(lme4)
library(MuMIn)
library(mgcv)

IYD_stop_fidelity$km_diff_1 <- abs(IYD_stop_fidelity$km_mark_1-IYD_stop_fidelity$km_prev_1)*1000
IYD_stop_fidelity$km_diff_2 <- abs(IYD_stop_fidelity$km_mark_2-IYD_stop_fidelity$km_prev_2)*1000
IYD_stop_fidelity$km_diff_3 <- abs(IYD_stop_fidelity$km_mark_3-IYD_stop_fidelity$km_prev_3)*1000


#glm
m1 <- lmer(scale(iyd_1) ~ scale(km_diff_1) + (1|curr_Y) + (1+scale(km_diff_1)|AID), data = IYD_stop_fidelity %>% filter(iyd_1 < 10000))
m2 <- lmer(scale(iyd_2) ~ scale(km_diff_2) + (1|curr_Y) + (1+scale(km_diff_2)|AID), data = IYD_stop_fidelity%>% filter(iyd_2 < 10000))
m3 <- lmer(scale(iyd_3) ~ scale(km_diff_3) + (1|curr_Y) + (1+scale(km_diff_3)|AID), data = IYD_stop_fidelity%>% filter(iyd_3 < 10000))

m11 <-lmer(scale(iyd_1) ~ scale(km_diff_1) + (1|curr_Y) + (1+scale(km_diff_1)|AID),  data = IYD_stop_fidelity %>% filter(iyd_1 >= 10000))
m22 <- lmer(scale(iyd_2) ~ scale(km_diff_2) + (1|curr_Y) + (1+scale(km_diff_2)|AID),  data = IYD_stop_fidelity%>% filter(iyd_2 >= 10000))
m33 <-lmer(scale(iyd_3) ~ scale(km_diff_3) + (1|curr_Y) + (1+scale(km_diff_3)|AID),  data = IYD_stop_fidelity%>% filter(iyd_3 >= 10000))

lapply(list(m1,m2,m3, m11, m22, m33), r.squaredGLMM)

#gams

head(IYD_stop_fidelity)



y1 <- ggplot( ) + geom_point(aes(x = (km_diff_1/1000),color = "#8B4513", y = (iyd_1/1000)), size = 3.6, data = IYD_stop_fidelity[,-c(6,16,26)] %>% filter(iyd_1 < 10000)) + geom_point(aes(x = (km_diff_1/1000), color = "#4682B4", y = (iyd_1/1000)), size = 3.6, data = IYD_stop_fidelity[,-c(6,16,26)] %>% filter(iyd_1 > 10000)) + geom_abline(slope = 1, intercept = 0, linewidth = 1.3) + scale_color_manual(labels = c("> 10km", "< 10km"), values = c("#4682B4","#8B4533")) + theme_classic() + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + labs(y = "Inter-year Distance", x = "Distance Along Route", color = "Inter-year Distance bin") + scale_x_continuous(expand = c(0, 0),limits = c(0,150), breaks = seq(0,150,50)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150), breaks = seq(0,150,50))

y2 <- ggplot( ) + geom_point(aes(x = (km_diff_2/1000), color = "#8B4523",y = (iyd_2/1000)), size = 3.6, data = IYD_stop_fidelity[,-c(6,26,26)] %>% filter(iyd_2 < 10000)) + geom_point(aes(x = (km_diff_2/1000), color = "#4682B4",y = (iyd_2/1000)), size = 3.6, data = IYD_stop_fidelity[,-c(6,26,26)] %>% filter(iyd_2 > 10000)) + geom_abline(slope = 1, intercept = 0, linewidth = 1.3) + theme_classic() + scale_color_manual(labels = c("> 10km", "< 10km"), values = c("#4682B4","#8B4533")) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"), legend.position = "none") + labs(x = "Route distance", y = "") + coord_cartesian(ylim = c(0,150), xlim = c(0,150)) + scale_x_continuous(expand = c(0, 0),limits = c(0,150), breaks = seq(0,150,50)) + scale_y_continuous(expand = c(0, 0),limits = c(0,150), breaks = seq(0,150,50))

y3 <- ggplot( ) + geom_point(aes(x = (km_diff_3/1000), color = "#8B4533", y = (iyd_3/1000)), size = 3.6, data = IYD_stop_fidelity[,-c(6,36,26)] %>% filter(iyd_3 < 10000)) + geom_point(aes(x = (km_diff_3/1000), color = "#4682B4",y = (iyd_3/1000)), size = 3.6, data = IYD_stop_fidelity[,-c(6,36,26)] %>% filter(iyd_3 > 10000)) + geom_abline(slope = 1, intercept = 0, linewidth = 1.3) + theme_classic() + scale_color_manual(labels = c("> 10km", "< 10km"), values = c("#4682B4","#8B4533")) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + labs(y = "", x = "", color = "Inter-year Distance bin") + scale_x_continuous(expand = c(0, 0),limits = c(0,150), breaks = seq(0,150,50)) + scale_y_continuous(expand = c(0, 0), limits = c(0,150), breaks = seq(0,150,50))

y1 #+ y2 + y3

ggsave(filename = "Figures/RouteDistvsIYD_20230507.jpg", width = 30, height = 24, units = "cm", dpi = 600)
