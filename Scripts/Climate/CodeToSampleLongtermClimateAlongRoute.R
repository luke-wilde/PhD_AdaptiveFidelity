#Code to simulate the alternative sites that could have been used, and extract iNDVI to them

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

#----------------------#
# Load data ####
setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
load("Data_out/Data/Stopover/RDH_StopoverBundle_14t22_20230507.RData")

load("Data_out/Data/Stopover/FidelityMetrics_20230507.RData")

`#----------------------#
# join with stops ####

head(data1)

data2 <- data1 %>% left_join(onstop_yr_fin %>% dplyr::select(id_yr, POSIXct, TimeStopped, stop.n.c) %>% mutate(OnStop = 1) %>% st_drop_geometry(), by = c("id_yr", "POSIXct"))
data2[which(is.na(data2$OnStop)),"OnStop"] <- 0
data2 <- data2 %>% mutate(km_rnd = round(km_mark,0))
head(data2)

clim <- data2 %>% group_by(km_rnd) %>% summarise(sd = sd(IntegratedNDVI), n = n())


ggplot(clim, aes(x=km_rnd, y=sd)) + geom_line()

#----------------------#
# turn to lines and sample climatology ####

moving <- data2 %>% filter(as.numeric(dist) > 2000)

moving <- st_as_sf(moving %>% ungroup() %>% as.data.frame(), coords = c("xend", "yend"), crs = st_crs(onstop_yr_fin)$input)

miglin <- Points2Lines(data = moving, date_name = "POSIXct", id_name = "id_yr", byid = T, no_cores = 4)
miglin <- miglin %>% mutate(Year = year(firstdate))
miglin

df <- data1[0,0]
for(i in 1:nrow(miglin)){
#i = 2
ex <- st_line_sample(x = miglin[i,], density = 1/1000, type = "regular")

ex <- st_as_sf(st_cast(ex, "POINT"))
st_geometry(ex) <- "geometry"
ex$id_yr <- rep(miglin[i,"id_yr"] %>% st_drop_geometry(),nrow(ex))
ex$Year <- rep(miglin[i,"Year"] %>% st_drop_geometry(),nrow(ex))
ex$km <- row_number(ex)
df <- rbind(df,ex)
}

df$J <- 132 # midpoint

mapview(miglin[1,]) + mapview(ex)

source("Z:/MODIS_NDVI/info/ndviEXTRACT2.R")


folder2 <- "Z:/MODIS_NDVI/NDVI2022"#"//DESKTOP-ug87ror.uwyo.edu//GIS_WyomingPlus//MODIS_NDVI//NDVI2022"
proj <- "+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs+towgs84=0,0,0"
metrics22 <- c("MaxNDVIDay", "MaxIRGday", "SpringStartDay", "SpringEndDay", 
               "IntegratedNDVI", "NDVI_scaled", "IRG_scaled", "MaxBrownDownDate", "SE_springDLCfit","SpringScale", "SpringLength", "sumIRG")

df[,c("x","y")] <- st_coordinates(df)



f<-as.data.frame(df); head(f)


for(i in 2000:2022){
  #i = 2000
  f$date <- as.POSIXct(paste(i,"-01-01 00:00:00",sep=""), format = "%Y-%m-%d %H:%M:%S", tz = "UTC") + f$J

  f$temp <- ndviEXTRACT2(XYdata=f, NDVImetric=metrics22[5], NDVIfolder=folder2,
                         maxcpus= detectCores()-1, xname="x", yname="y", datesname="date",
                         xyCRS=CRS(proj))
  names(f)[names(f)=="temp"] <- paste(i,"_indvi",sep="")
  
}#i

head(f)
#beep(sound = 1)
showConnections(); sfStop(); showConnections()

f$sd <- apply(f[,c(9:31)], 1, sd, na.rm=TRUE)
f$mean <- apply(f[,c(9:31)], 1, mean, na.rm=TRUE)

ggplot(data = f, aes(sd)) + geom_density()
summary(f$sd)
f$id_yr <- do.call(rbind, f$id_yr)
f$Year <- do.call(rbind, f$Year)
str(f)

attr(f$id_yr, "ATT") <- NULL
attr(f$Year, "ATT") <- NULL
attr(f$sd, "ATT") <- NULL

#clim.var.data <- data2 %>% rename(km = km_rnd) %>% left_join(f, by = c("id_yr", "km"))

table(clim.var.data$OnStop) #use cloglog bc I dont have even #s of 0 and 1

#run the model
mod.stop.status <- glmmTMB(OnStop ~ sd + (1|AID) + (1|AID) + (1|id_yr), data = clim.var.data, family = binomial(link = "cloglog"))
res <- DHARMa::simulateResiduals(mod.stop.status)
plot(res)
summary(mod.stop.status)

#predict
new1 <- data.frame(sd = seq(min(na.omit(clim.var.data$sd)),max(na.omit(clim.var.data$sd)),length = 50), AID = "108" , id_yr = "108_2014" , Year = 2014)
P1 <- predict(mod.stop.status, newdata = new1, type = "link", se=T)

clog.inv<-function(values) {
  e1<-(-1)*(exp(values))
  1-exp(e1)
}

lw <- clog.inv(P1$fit - 1.96*P1$se)
up <- clog.inv(P1$fit + 1.96*P1$se)
mn <- clog.inv(P1$fit)

#linear response, but not strong influence 
ggplot() + geom_point(aes(x = sd, y = OnStop), data = clim.var.data %>% sample_n(5000)) + geom_line(aes(y = as.vector(mn), x = seq(min(na.omit(clim.var.data$sd)),max(na.omit(clim.var.data$sd)),length = 50)), color = "blue") + theme_bw()

# slight response 


#----------------------#
# test for fidelity response to variability ####

head(f)

dat <- IYD_stop_fidelity_list[[1]] %>% mutate(km = round(km_mark_1, 0)) %>% left_join(f %>% dplyr::select(id_yr, km, sd, mean), by = c("id_yr", "km"))

dat <- dat %>% mutate(fidc = cut(iyd_1, breaks = quantile(iyd_1, c(0,.5,1.00)),include.lowest = TRUE, labels = FALSE), sdc = scale(sd), kmcut = ifelse(km_mark_1 <= 80, "start", ifelse(km_mark_1 > 80 & km_mark_1 < 160, "middle", "end")))

library(ggridges)

ggplot(dat) + stat_density_ridges(aes(x = sd, fill = factor(fidc), y = kmcut, height = after_stat(density)), color = "black", linewidth = 1.2, alpha = .25, quantile_lines = TRUE, quantiles = 2, scale = .95) + theme_classic()+ labs(y = "Position En-route", x = "St. Dev. 22-year integrated-NDVI (scaled)", fill = "Fidelity", color = "Fidelity", stat = "density")  + scale_x_continuous(limits = c(-1.05,8), breaks = seq(-1,8,1)) + theme(axis.title.x = element_text(size = 36,color = "grey18"), axis.title.y = element_text(size = 36,color = "grey18"),axis.text.x = element_text(size = 30,color = "grey18"),axis.text.y = element_text(size = 30,color = "grey18"), legend.text = element_text(size = 30,color = "grey18"), legend.title = element_text(size = 32,color = "grey18"))  + geom_text(x=6, y=3.5, label="t(203)=0.71, p=0.48", family = "serif", size = 12) + geom_text(x=6, y=2.5, label="t(369)=0.25, p=0.79", family = "serif", size = 12) + geom_text(x=6, y=1.5, label="t(157)=0.34, p=0.73", family = "serif", size = 12) + scale_fill_manual(values = c("#FF4500","#008080"), labels = c("Above Average", "Below Average")) + scale_color_manual(values = c("#FF4500","#008080"), labels = c("Above Average", "Below Average"))

# route1 <- ggplot(dat %>% filter(km < 80)) + geom_vline(aes(xintercept = sd, group = factor(fidc), color = factor(fidc)),  linewidth = 1.2, data = dat %>% filter(km < 80) %>% group_by(fidc) %>% summarise(sd = mean(na.omit(sd)))) + geom_density(aes(x = sd, fill = factor(fidc), y = after_stat(scaled)), color = "black", linewidth = 1.2, alpha = .25)  + theme_classic()+ theme(legend.position = "none") + labs(y = "Density", x = "")  + scale_x_continuous(limits = c(-1.05,8), breaks = seq(-1,8,1)) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 14,color = "grey18"), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18"))  + geom_text(x=-.25, y=.8, label="t(203)=0.71, p=0.48", family = "serif", size = 10) + scale_fill_manual(values = c("#FF4500","#008080"), labels = c("Above Average", "Below Average")) + scale_color_manual(values = c("#FF4500","#008080"), labels = c("Above Average", "Below Average"))
# 
# route2 <-ggplot(dat %>% filter(km > 80 & km < 160)) + geom_vline(aes(xintercept = sd, group = factor(fidc), color = factor(fidc)),  linewidth = 1.2, data = dat %>% filter(km > 80 & km < 160) %>% group_by(fidc) %>% summarise(sd = mean(na.omit(sd)))) + geom_density(aes(x = sd, fill = factor(fidc), y = after_stat(scaled)), color = "black", linewidth = 1.2, alpha = .3) + theme_classic() + labs(y = "Density", x = "", fill = "Fidelity", color = "Fidelity") + scale_x_continuous(limits = c(-1.05,8), breaks = seq(-1,8,1)) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 14,color = "grey18"), legend.text = element_text(size = 26,color = "grey18"), title = element_text(size = 20,color = "grey18")) + geom_text(x=-.25, y=.8, label="t(369)=0.25, p=0.79", family = "serif", size = 10) + scale_fill_manual(values = c("#FF4500","#008080"), labels = c("Above Average", "Below Average")) + scale_color_manual(values = c("#FF4500","#008080"), labels = c("Above Average", "Below Average"))
# 
# route3 <-ggplot(dat %>% filter(km > 160 & km < 240)) + geom_vline(aes(xintercept = sd, group = factor(fidc), color = factor(fidc)),  linewidth = 1.2, data = dat %>% filter(km > 160 & km < 240) %>% group_by(fidc) %>% summarise(sd = mean(na.omit(sd)))) + geom_density(aes(x = sd, fill = factor(fidc), y = after_stat(scaled)), color = "black", linewidth = 1.2, alpha = .3) + theme_classic() + labs(y = "Density", x = "St. Dev. 22-year integrated-NDVI (scaled)") + theme(legend.position = "none") + scale_x_continuous(limits = c(-1.05,8), breaks = seq(-1,8,1)) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 14,color = "grey18"), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + geom_text(x=-.25, y=.8, label="t(157)=0.34, p=0.73", family = "serif", size = 10) + scale_fill_manual(values = c("#FF4500","#008080"), labels = c("Above Average", "Below Average")) + scale_color_manual(values = c("#FF4500","#008080"), labels = c("Above Average", "Below Average"))
# 
# route1 + route2 + route3 + plot_layout(ncol = 1, nrow = 3)

ggsave(filename = "Figures/LongTermBiomassProd_ridge_20230507.jpg", width = 48, height = 36, units = "cm", dpi = 600)

library(lme4)

summary(glm(sd ~ factor(fidc), family = Gamma(link = "log"), data = dat %>% filter(km < 80)))
summary(glm(sd ~ factor(fidc), family = Gamma(link = "log"), data = dat %>% filter(km > 80 & km < 160)))
summary(glm(sd ~ factor(fidc), family = Gamma(link = "log"), data = dat %>% filter(km > 160 & km < 240)))

summary(lm(mean ~ factor(fidc), data = dat %>% filter(km < 80)))
summary(lm(mean ~ factor(fidc), data = dat %>% filter(km > 80 & km < 160)))
summary(lm(mean ~ factor(fidc), data = dat %>% filter(km > 160 & km < 240)))




high_1 <- IYD_HighUse_stop_fidelity_list[[1]] %>% mutate(km = round(km_mark_1, 0))
all_1 <- IYD_stop_fidelity_list[[1]] %>% mutate(km = round(km_mark_1, 0))

fidH1 <- high_1 %>% left_join(f %>% st_drop_geometry() %>% dplyr::select(sd, id_yr, km), by = c("id_yr", "km"))
fidA1 <- all_1 %>% left_join(f %>% st_drop_geometry() %>% dplyr::select(sd, id_yr, km), by = c("id_yr", "km"))

high_2 <- IYD_HighUse_stop_fidelity_list[[2]] %>% mutate(km = round(km_mark_2, 0))
all_2 <- IYD_stop_fidelity_list[[2]] %>% mutate(km = round(km_mark_2, 0))

fidH2 <- high_2 %>% left_join(f %>% st_drop_geometry() %>% dplyr::select(sd, id_yr, km), by = c("id_yr", "km"))
fidA2 <- all_2 %>% left_join(f %>% st_drop_geometry() %>% dplyr::select(sd, id_yr, km), by = c("id_yr", "km"))

high_3 <- IYD_HighUse_stop_fidelity_list[[3]] %>% mutate(km = round(km_mark_3, 0))
all_3 <- IYD_stop_fidelity_list[[3]] %>% mutate(km = round(km_mark_3, 0))

fidH3 <- high_3 %>% left_join(f %>% st_drop_geometry() %>% dplyr::select(sd, id_yr, km), by = c("id_yr", "km"))
fidA3 <- all_3 %>% left_join(f %>% st_drop_geometry() %>% dplyr::select(sd, id_yr, km), by = c("id_yr", "km"))


#models
head(fidH1)
fidH1[which(fidH1$iyd_1 == 0),] <- 0.000000001
mod.FidHighStopClim.Y1 <- glmmTMB(iyd_1 ~ sd + (1|AID) + (1|curr_Y) + (1|id_yr), data = fidH1, family = Gamma(link="log"))

plot(DHARMa::simulateResiduals(mod.FidHighStopClim.Y1))
summary(mod.FidHighStopClim.Y1)


new <- data.frame()
predict(type = "response")

head(fidA1)
fidA1[which(fidA1$iyd_1 < 1),] <- 1
mod.FidAllStopClim.Y1 <- glmmTMB(iyd_1 ~ sd + (1|AID) + (1|curr_Y) + (1|id_yr), data = fidA1, family = Gamma(link="log"))

plot(DHARMa::simulateResiduals(mod.FidAllStopClim.Y1))
summary(mod.FidAllStopClim.Y1)

new <- data.frame(sd = seq(min(na.omit(fidA1$sd)),max(na.omit(fidA1$sd)),length=50), AID = unique(fidA1$AID)[13], curr_Y = unique(fidA1$curr_Y)[13], id_yr = unique(fidA1$id_yr)[13])
new$fit <- predict(mod.FidAllStopClim.Y1, new, type = "response", se.fit = T,allow.new.levels=TRUE)$fit
new$sefit <- predict(mod.FidAllStopClim.Y1, new, type = "response", se.fit = T,allow.new.levels=TRUE)$se.fit


ggplot(data = fidA1, aes(x = sd, y = 1/iyd_1)) + geom_point() + geom_line(aes(x = sd, y = 1/fit), data = new)
ggplot() + geom_point() + geom_line(aes(x = sd, y = 1/fit), data = new)

fidA1$FidIndex <- ifelse(fidA1$iyd_1 > 5000, "0", "1")


plot1 <- ggplot(fidA1) + geom_density(aes(x = sd, fill = FidIndex, y = after_stat(density)), alpha = .5, color = "grey45") + theme_classic() + scale_fill_manual(values = c("#87CEEB","#B22222"), labels = c("High Fidelity", "Low Fidelity")) + labs(x = expression(paste("22-year integrated-NDVI St. Dev. ", km^{-1})), y = "Density", fill = "") + coord_cartesian(xlim = c(0,22), ylim = c(0,0.5)) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"))

#ggsave(filename = "./Figures/StopoverResponseToLongterm_iNDVI.svg", dpi = 500, width = 15, height = 10, units = "cm")

#save.image("./REnvs/ImageOfClimateSampling_20230428.RDS")

#----------------------#
# biomass binned ####

readRDS("REnvs/ImageOfClimateSampling_20230428.rds")

head(f)

fsimpl <- f %>% group_by(km) %>%  sample_n(1)

ggplot(fsimpl) + geom_line(aes(x = km, y = mean), color = "") + geom_errorbar(aes(x = km, ymin = mean - sd, ymax = mean + sd))
