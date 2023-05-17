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

load("Data_out/Data/Stopover/FidelityDataset_HighUse_20230513.RData")
load("Data_out/Data/Stopover/FidelityMetrics_HighUse_20230513.RData")

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
beep(sound = 1)
showConnections(); sfStop(); showConnections()


data <- f



setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
load("C:/Users/lwilde2/Desktop/RDH Database/Processed Data/RDH_AllMigrations_Bischof_2014t2022_20230425.RData")



Spring.summary <- Spring.summary %>% mutate(Year = year(start), J = yday(start))  %>% drop_na(J)

Spring.sum.c <- Spring.summary %>% filter(Year > 2015) %>% group_by(Year) %>% mutate(cut = cut(J, breaks = quantile(J, c(0,.15,.45,1)),include.lowest = TRUE, labels = FALSE))


# summary(Spring.summary[which(Spring.summary$Year == 2014),"J"])
# summary(Spring.summary[which(Spring.summary$Year == 2015),"J"])
# summary(Spring.summary[which(Spring.summary$Year == 2016),"J"])
# summary(Spring.summary[which(Spring.summary$Year == 2017),"J"])
# summary(Spring.summary[which(Spring.summary$Year == 2018),"J"])
# summary(Spring.summary[which(Spring.summary$Year == 2019),"J"])
# summary(Spring.summary[which(Spring.summary$Year == 2020),"J"])
# summary(Spring.summary[which(Spring.summary$Year == 2021),"J"])
# summary(Spring.summary[which(Spring.summary$Year == 2022),"J"])


Spring.sum.c1 <- Spring.summary %>% filter(Year <= 2015) %>% group_by(Year) %>% mutate(cut = cut(J, breaks = quantile(J, c(0,.5,1)),include.lowest = TRUE, labels = FALSE))

sumall <- rbind(Spring.sum.c, Spring.sum.c1)

Spring.summary[which(is.na(Spring.summary$J)),]      


sumall %>% group_by(cut) %>% summarise(n = n())




data <- data %>% ungroup() %>% mutate(DFP_alt = JDate - MaxIRGday_alt, absDFP_alt = abs(DFP_alt),diffDFP = DFP - DFP_alt, diffabsDFP = absDFP - absDFP_alt, diffIRGday = MaxIRGday - MaxIRGday_alt, diffSpring = SpringLength - SpringLength_alt) %>% ungroup()  

datasum <- data %>% group_by(stop.n.c) %>% dplyr::summarise(diffDFP = mean(diffDFP), diffabsDFP = mean(diffabsDFP), diffIRGday = mean(diffIRGday), diffSpring = mean(diffSpring), iyd = mean(iyd),km_mark_1 = unique(km_mark_1), km_prev_1 = unique(km_prev_1), diffIRGsum = sum(IRG_scaled) - sum(IRG_scaled_alt))


datasum

datasum <- datasum %>% mutate( AID = str_sub(stop.n.c, 0,3), id_yr = str_sub(stop.n.c, 0,8), Year = str_sub(stop.n.c, 5, 8), stop = as.numeric(str_split(stop.n.c, "_", simplify = TRUE)[,3]))

datasum <- datasum %>% left_join(sumall %>% mutate(Year = as.character(Year)) %>% dplyr::select(id_yr, cut), by="id_yr")

datasum$cut <- ifelse(datasum$cut == 1, "early", ifelse(datasum$cut == 2, "mid", "late"))

datasum %>% group_by(cut) %>% summarise(n = n())


# ggplot(datasum, aes(x = km_mark_1, y = iyd, group = cut, color = cut)) + geom_point() + geom_smooth(method = "loess") + coord_cartesian(xlim = c(35,240),ylim = c(0,120000))
# 
# 
# datasum[which(datasum$iyd == 0), "iyd"] <- 0.00000001
# 
# mod.diffsites <- mgcv::gamm(round(iyd,0) ~ s(km_mark_1, by = factor(cut), bs = "cs", k = 5), random = list(Year.x = ~1, AID = ~1, id_yr = ~1), family = poisson(link = "log"), data = datasum )
# 
# summary(mod.diffsites$gam)
# 
# library(gammit)
# 
# newdata = data.frame(expand.grid(km_mark_1 = seq(0, 240, length = 240), cut = c("early", "mid", "late"), AID = 108, Year.x = 2019, id_yr = "108_2019"))
# pred.diffsite <- gammit::predict_gamm(mod.diffsites$gam , newdata = newdata, se = T)
# 
# ggplot(cbind(newdata, pred.diffsite)) + geom_line(aes(x = km_mark_1, y = exp(prediction), color = cut)) + coord_cartesian(ylim = c(0, 80000))


head(datasum)

datasum$diffkm = abs(datasum$km_mark_1 - datasum$km_prev_1)
datasum$direction = ifelse(datasum$km_mark_1 - datasum$km_prev_1 > 0, 1, -1)
summary(datasum$diffkm)

#rules
# within 5 km then the same stop
# > 5 km diff



datasum <- datasum %>% mutate(iyd.signed = iyd * direction, iyd.signed.c = ifelse(iyd > 5000, 1, NA), kmcut = ifelse(km_mark_1 <= 80, "start", ifelse(km_mark_1 > 80 & km_mark_1 < 160, "middle", "end")), status = iyd.signed.c ) #, kmc = ifelse(diffkm < 5, 1, ifelse(diffkm > 5 & diffkm < 10, 2, 3)) + kmc

#codes
#2 = same , 3 = shift , 4 = shift, 5 = shift, 6 = skip

unique(datasum$status)

datasum <- datasum %>% drop_na(iyd.signed.c)





# datasum <- datasum %>% mutate(status = ifelse(iyd < 10000 & diffkm < 10000, ifelse(iyd < 1000,"use","alt"),"skip"), kmcut = ifelse(km_mark_1 <= 80, "start", ifelse(km_mark_1 > 80 & km_mark_1 < 160, "middle", "end")))



library(ggridges)

ggplot(datasum %>% filter(kmcut == "start")) + theme_classic() +  stat_density_ridges(aes(x = diffDFP, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9) +
ggplot(datasum %>% filter(kmcut == "middle")) + theme_classic() +  stat_density_ridges(aes(x = diffDFP, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+
ggplot(datasum %>% filter(kmcut == "end")) + theme_classic() +  stat_density_ridges(aes(x = diffDFP, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)


ggplot(datasum %>% filter(kmcut == "start")) + theme_classic() +  stat_density_ridges(aes(x = diffabsDFP, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9) +
  ggplot(datasum %>% filter(kmcut == "middle")) + theme_classic() +  stat_density_ridges(aes(x = diffabsDFP, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+
  ggplot(datasum %>% filter(kmcut == "end")) + theme_classic() +  stat_density_ridges(aes(x = diffabsDFP, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)

head(datasum)
unique(datasum$status)

stat_func <- function(data) {
  group_values <- data
  # Change the statistic of interest as needed
  mean(group_values)
}

stat_func(datasum$diffabsDFP)

# Use group_by and summarize to calculate bootstrap confidence intervals per group
library(dplyr)
library(Hmisc)

boot <- datasum %>% 
  mutate(km_rnd = round(km_mark_1/20,0)*20) %>% 
  dplyr::select(km_rnd, cut, status, diffabsDFP) %>% 
  group_by(km_rnd, cut, status) %>% 
  group_map(~ smean.cl.boot(., conf.int = .95, B = 100, na.rm = TRUE)) %>%
  bind_rows()

dat <- datasum %>% 
  mutate(km_rnd = round(km_mark_1/20,0)*20) %>%
  dplyr::select(km_rnd, cut, status, diffabsDFP) %>% 
  group_by(km_rnd, cut, status) %>% summarise(Mean = mean(diffabsDFP))

dat <- dat %>% left_join(boot, by = "Mean")

data_new <- dat                             # Replicate data
data_new$group <- factor(data_new$cut,     
                         levels = c("early", "mid", "late")) # Reordering group factor levels
#data_new$x <- factor(data_new$kmcut, levels = c("start", "middle", "end")) 

ggplot(data_new %>% filter(km_rnd < 240)) + geom_pointrange(aes(x = km_rnd, y = Mean, ymin = Lower, ymax = Upper, group = factor(status), color = factor(status)), position = position_dodge(width = .5), size = 1.75, linewidth = 1.3) + geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_line(aes(x = km_rnd, y = Mean, group = factor(status), color = factor(status)), position = position_dodge(width = .5), size = 2.15) + geom_hline(aes(yintercept = 0), linetype = "dashed") + labs(y = "Surfing Improvement From Remembered Site", x = "Distance from Winter Range (km)", color = "Stop Status") + scale_color_manual(values = c("#3C5488FF"), labels = c("Short")) + theme_classic() + facet_wrap(~ group) + theme(axis.title.x = element_text(size = 48,color = "grey18"), axis.title.y = element_text(size = 48,color = "grey18"),axis.text.x = element_text(size = 40,color = "grey18"),axis.text.y = element_text(size = 40,color = "grey18"),legend.text = element_text(size = 30,color = "grey18"), legend.title = element_text(size = 40,color = "grey18"),legend.key.size = unit(2, 'cm'), legend.key.height = unit(1, 'cm'), legend.position = "none", legend.key.width = unit(2, 'cm'), strip.text = element_text(size = 30,color = "grey18")) + coord_cartesian(ylim = c(-40,40)) + scale_y_continuous(limits = c(-60,60), breaks = seq(-40,40,10)) + scale_x_continuous(limits = c(0,240), breaks = seq(0,240,60))

ggsave(filename = "C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/Figures/ShortAndSkippedSites_20230516.jpg", width = 72, height = 42, units = "cm", dpi = 600)








datasum %>% mutate(iyd, )

ggplot(datasum) + geom_point(aes(x = km_mark_1, y = diffabsDFP, size = iyd, color = cut)) + geom_smooth(aes(x = km_mark_1, y = diffabsDFP, group = cut, color = cut), method = "loess") + theme_classic()  + theme(axis.title.x = element_text(size = 40,color = "grey18"), axis.title.y = element_text(size = 40,color = "grey18"),axis.text.x = element_text(size = 30,color = "grey18"),axis.text.y = element_text(size = 30,color = "grey18"),legend.text = element_text(size = 30,color = "grey18"), legend.title = element_text(size = 40,color = "grey18"),legend.key.size = unit(2, 'cm'), legend.key.height = unit(1, 'cm'), legend.key.width = unit(2, 'cm'), strip.text = element_text(size = 30,color = "grey18")) + scale_y_continuous(limits = c(-30,30), breaks = seq(-30,30,10)) + geom_hline(aes(yintercept = 0), linetype = "dashed") + labs(y = "Surfing Improvement From Remembered Site", x = "Section of Corridor", color = "Stop Status") + scale_color_manual(values = c("#3C5488FF","#DC0000FF"), labels = c("Short", "Skipped")) + theme_classic() + facet_wrap(~ group) + theme(axis.title.x = element_text(size = 60,color = "grey18"), axis.title.y = element_text(size = 60,color = "grey18"),axis.text.x = element_text(size = 50,color = "grey18"),axis.text.y = element_text(size = 50,color = "grey18"),legend.text = element_text(size = 50,color = "grey18"), legend.title = element_text(size = 40,color = "grey18"),legend.key.size = unit(2, 'cm'), legend.key.height = unit(1, 'cm'), legend.key.width = unit(2, 'cm'), strip.text = element_text(size = 30,color = "grey18")) + scale_y_continuous(limits = c(-30,30), breaks = seq(-30,30,10))






ggplot(datasum %>% filter(kmcut == "start")) + geom_density_ridges2(aes(x = diffIRGday, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic() +
  ggplot(datasum %>% filter(kmcut == "middle")) + geom_density_ridges2(aes(x = diffIRGday, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic() +
  ggplot(datasum %>% filter(kmcut == "end")) + geom_density_ridges2(aes(x = diffIRGday, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic()


ggplot(datasum %>% filter(kmcut == "start")) + geom_density_ridges2(aes(x = diffSpring, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic() +
  ggplot(datasum %>% filter(kmcut == "middle")) + geom_density_ridges2(aes(x = diffSpring, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic() +
  ggplot(datasum %>% filter(kmcut == "end")) + geom_density_ridges2(aes(x = diffSpring, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic()


ggplot(datasum %>% filter(kmcut == "start")) + geom_density_ridges2(aes(x = diffIRGsum, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic()+
ggplot(datasum %>% filter(kmcut == "middle")) + geom_density_ridges2(aes(x = diffIRGsum, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic()+
ggplot(datasum %>% filter(kmcut == "end")) + geom_density_ridges2(aes(x = diffIRGsum, y = factor(cut), fill = factor(status), color = factor(status)), quantile_lines = TRUE, quantiles = 2, alpha = .3, scale = 0.9)+ theme_classic()


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
