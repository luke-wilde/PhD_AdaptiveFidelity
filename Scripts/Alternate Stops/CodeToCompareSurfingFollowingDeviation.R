# code to investigate the origins of fidelity
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
library(rptR)
library(glmmTMB)
#install.packages("sm")
library(sm)


load("REnvs/EnvForDeviationsAnalysis_20230509.RDS")

#----------------------#
# Load data ####
setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
load("Data_out/Data/Stopover/RDH_StopoverBundle_14t22_20230507.RData")

load("Data_out/Data/Stopover/FidelityMetrics_20230507.RData")

#
onstop_yr_fin <- onstop_yr_fin %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))


#pull out the used stops and those that were skipped, the nearest neighbor
used <- onstop_yr_fin %>% filter(stop.n.c %in% IYD_stop_fidelity_list[[1]]$stop.n_1)
alt <- onstop_yr_fin %>% filter(stop.n.c %in% IYD_stop_fidelity_list[[1]]$stop_prev_1)

#find centroid of used
alt.cent <- alt %>% group_by(stop.n.c) %>% summarise() %>% st_centroid()
alt.cent[,c("xalt","yalt")] <- st_coordinates(alt.cent)
attr(alt.cent$stop.n.c, "ATT") <- NULL
names(alt.cent)[1] <- "stop_prev_1"

current <- used %>% left_join(IYD_stop_fidelity_list[[1]] %>% rename(iyd = iyd_1, stop.n.c = stop.n_1), by = "stop.n.c")

current <- current %>% left_join(alt.cent %>% st_drop_geometry() %>% as.data.frame(), by = "stop_prev_1")

#calculate current centroid
curr.cent <- current %>% ungroup() %>% group_by(stop.n.c) %>% summarise() %>% st_centroid()
curr.cent[,c("xcurr","ycurr")] <- st_coordinates(curr.cent)

current <- current %>% left_join(curr.cent %>% st_drop_geometry() %>% as.data.frame(), by = "stop.n.c")

names(current)

current <- current %>% mutate(xshift = xcurr - xalt, yshift = ycurr - yalt)
b <- st_as_sf(current  %>% st_drop_geometry() %>% dplyr::select(xshift, yshift), coords = c("xshift", "yshift"), crs = st_crs(current)$input)
alternate <- st_as_sf(current$geometry + b$geometry, crs = st_crs(current)$input)
alternate <- cbind(alternate, current %>% st_drop_geometry() %>% dplyr::select(-c("MaxNDVIDay", "SpringStartDay", "SpringEndDay", "IntegratedNDVI", "NDVI_scaled", "MaxBrownDownDate", "SE_springDLCfit","SpringScale", "sumIRG")))


st_geometry(alternate) <- "geometry"
alternate[,c("x","y")] <- st_coordinates(alternate)

#extract ndvi information
source("Z:/MODIS_NDVI/info/ndviEXTRACT2.R")


folder2 <- "Z:/MODIS_NDVI/NDVI2022"#"//DESKTOP-ug87ror.uwyo.edu//GIS_WyomingPlus//MODIS_NDVI//NDVI2022"
proj <- "+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs+towgs84=0,0,0"
metrics22 <- c("MaxNDVIDay", "MaxIRGday", "SpringStartDay", "SpringEndDay", 
               "IntegratedNDVI", "NDVI_scaled", "IRG_scaled", "MaxBrownDownDate", "SE_springDLCfit","SpringScale", "SpringLength", "sumIRG")





f<-alternate %>% st_drop_geometry() %>% as.data.frame(); head(f)


for(i in c(2,7,11)){
  #i = 2000
  #f$date <- as.POSIXct(paste(i,"-01-01 00:00:00",sep=""), format = "%Y-%m-%d %H:%M:%S", tz = "UTC") + f$J
  
  f$temp <- ndviEXTRACT2(XYdata=f, NDVImetric=metrics22[i], NDVIfolder=folder2,
                         maxcpus= detectCores()-1, xname="x", yname="y", datesname="POSIXct",
                         xyCRS=CRS(proj))
  names(f)[names(f)=="temp"] <- paste(metrics22[i],"_alt",sep="")
  
}#i

head(f)
#beep(sound = 1)
showConnections(); sfStop(); showConnections()


data <- f

data <- data %>% ungroup() %>% mutate(DFP_alt = JDate - MaxIRGday_alt, absDFP_alt = abs(DFP_alt),diffDFP = DFP - DFP_alt, diffabsDFP = absDFP - absDFP_alt, diffIRGday = MaxIRGday - MaxIRGday_alt, diffSpring = SpringLength - SpringLength_alt) %>% ungroup()  

datasum <- data %>% group_by(stop.n.c) %>% dplyr::summarise(diffDFP = mean(diffDFP), diffabsDFP = mean(diffabsDFP), diffIRGday = mean(diffIRGday), diffSpring = mean(diffSpring), iyd = mean(iyd),km_mark_1 = unique(km_mark_1), km_prev_1 = unique(km_prev_1), diffIRGsum = sum(IRG_scaled) - sum(IRG_scaled_alt))



head(datasum)

datasum$diffkm = abs(datasum$km_mark_1 - datasum$km_prev_1)

summary(datasum$diffkm)

#rules
# within 5 km then the same stop
# > 5 km diff


datasum <- datasum %>% mutate(iydc= ifelse(iyd < 5000, 1, ifelse(iyd > 5000 & iyd < 10000, 2, 3)), kmc = ifelse(diffkm < 5, 1, ifelse(diffkm > 5 & diffkm < 10, 2, 3)), kmcut = ifelse(km_mark_1 <= 80, "start", ifelse(km_mark_1 > 80 & km_mark_1 < 160, "middle", "end")), status = iydc + kmc)

#codes
#2 = same , 3 = shift , 4 = shift, 5 = shift, 6 = skip

unique(datasum$iydc)
unique(datasum$kmc)
unique(datasum$status)

datasum$status <- ifelse(datasum$status == 2, "use", ifelse(datasum$status < 6 & datasum$status > 2, "alt", "skip"))

# datasum <- datasum %>% mutate(status = ifelse(iyd < 10000 & diffkm < 10000, ifelse(iyd < 1000,"use","alt"),"skip"), kmcut = ifelse(km_mark_1 <= 80, "start", ifelse(km_mark_1 > 80 & km_mark_1 < 160, "middle", "end")))

datasum <- datasum %>% mutate( AID = str_sub(stop.n.c, 0,3), id_yr = str_sub(stop.n.c, 0,8), Year = str_sub(stop.n.c, 5, 8), stop = as.numeric(str_split(stop.n.c, "_", simplify = TRUE)[,3]))

library(ggridges)

ggplot(datasum) + theme_classic() +  stat_density_ridges(aes(x = diffDFP, y = factor(kmcut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)
ggplot(datasum) + geom_density_ridges2(aes(x = diffabsDFP, y = factor(kmcut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic()
ggplot(datasum) + geom_density_ridges2(aes(x = diffIRGday, y = factor(kmcut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic()
ggplot(datasum) + geom_density_ridges2(aes(x = diffSpring, y = factor(kmcut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic()
ggplot(datasum) + geom_density_ridges2(aes(x = diffIRGsum, y = factor(kmcut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic()


idyrsum <- datasum %>% group_by(id_yr,status, kmcut) %>% summarise_at(c("diffDFP", "diffabsDFP" ,"diffIRGday" ,"diffSpring","diffIRGsum"),mean)




#------------------------#
# chi sq of profiles ####
chitest <- datasum %>% dplyr::select(AID, id_yr, Year, stop, km_mark_1)


#calculate migration distance
moving <- data1 %>% filter(as.numeric(dist) > 2000)

moving <- st_as_sf(moving %>% ungroup() %>% as.data.frame(), coords = c("xend", "yend"), crs = st_crs(onstop_yr_fin)$input)

miglin <- Points2Lines(data = moving, date_name = "POSIXct", id_name = "id_yr", byid = T, no_cores = 4)
miglin <- miglin %>% mutate(Year = year(firstdate))
miglin$length = as.numeric(st_length(miglin))


chitest <- chitest %>% left_join(miglin %>% mutate(length = length / 1000) %>% st_drop_geometry() %>% select(id_yr, length), by = "id_yr")

chitest <- chitest %>% mutate(kmcut = 10*round(km_mark_1 / length, 1))

chitest2 <- chitest %>% group_by(id_yr,kmcut) %>% summarise(count = as.numeric(n())) %>% mutate(freq = count / sum(count))
chitest2 <- lapply(chitest2, function(x) { attributes(x) <- NULL; x })
chitest2 <- as.data.frame(do.call(cbind, chitest2))

str(chitest2)

chitest2$kmcut <- as.numeric(chitest2$kmcut)
chitest2$freq <- as.numeric(chitest2$freq)


chitest2 <- tidyr::complete(chitest2,id_yr, kmcut = 1:10, fill = list(freq = 0))

chitest2 <- chitest2 %>% mutate( AID = str_sub(id_yr, 0,3), Year = str_sub(id_yr, 5, 8))

head(chitest2)

chitest2$count <- as.numeric(chitest2$count)
chitest2[which(is.na(chitest2$count)),"count"] <- 0

#sequences

reptab <- data.frame(kmcut = rep(chitest2$kmcut,times= chitest2$count), AID = rep(chitest2$AID,times= chitest2$count),id_yr = rep(chitest2$id_yr,times= chitest2$count), Year = rep(chitest2$Year,times= chitest2$count))

intersect_all <- function(a,b,...){
  all_data <- c(a,b,...)
  require(plyr)
  count_data<- length(list(a,b,...))
  freq_dist <- count(all_data)
  intersect_data <- freq_dist[which(freq_dist$freq==count_data),"x"]
  intersect_data
}

intersect_all(reptab[which(reptab$Year == 2022),"kmcut"])
intersect_all(reptab[which(reptab$AID == 2022),"kmcut"])



library(stringdist)
seq_sim(a = reptab[1:7,1], b = reptab[8:18,1], method="cosine", q = 2)

wilcox.test(reptab[1:7,1], reptab[8:18,1], alternative = "two.sided", exact = FALSE)
wilcox.test(reptab[1:7,1], reptab[214:227,1], alternative = "two.sided", exact = FALSE)
wilcox.test(reptab[228:236,1], reptab[214:227,1], alternative = "two.sided", exact = FALSE)



d108 <- ggplot(chitest2 %>% filter(AID == "108"), aes(x = kmcut*10, fill = id_yr, y = freq))  + scale_x_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,10))  + stat_smooth( geom = 'area', method = 'loess', alpha = 1/2, aes(fill = id_yr)) + scale_y_continuous(expand = c(0,0), limits = c(0,.3), breaks = seq(0,.3,.1)) + labs(x = "Progress along Route (%)", y = "Relative Frequency", fill = "Animal-Year") + theme_classic()  # + geom_area(aes(color = id_yr), alpha = .3) + geom_bar(stat = "identity", position = position_dodge2(width = .2, preserve = "single"))

d255 <- ggplot(chitest2 %>% filter(AID == "255"), aes(x = kmcut*10, fill = id_yr, y = freq))  + scale_x_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,10))  + stat_smooth( geom = 'area', method = 'loess', alpha = 1/2, aes(fill = id_yr)) + scale_y_continuous(expand = c(0,0), limits = c(0,.3), breaks = seq(0,.3,.1)) + labs(x = "Progress along Route (%)", y = "Relative Frequency", fill = "Animal-Year") + theme_classic()

d108 + d255


# sm.density.compare(x = chitest2$kmcut,
#                    group = chitest2$id_yr,
#                    model = "equal")
