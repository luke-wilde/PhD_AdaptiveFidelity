# Script for modelling

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
library(survival)
library(MASS)
library(car)
library(dummies)
library(ggeffects)
library(mgcv)


#--------------------------#
# prep datasets ####

setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
load("Data_out/Data/Stopover/RDH_StopoverBundle_14t22_20230425.RData")
load("Data_out/Data/Stopover/FidelityMetrics_20230428.RData")
load("Data_out/Data/Stopover/RDH_AlternateStops_Bischof_20230428.RData")
load("Data_out/Data/Greenscape/RDH_14t22_greenscapes_20230425.RData")

# 
#combine these data
onstop_yr_fin$flag <- 1
# 
# AlternateStops_Bischof <- AlternateStops_Bischof %>% st_as_sf(coords = c("x","y"), crs = st_crs(onstop_yr_fin)$input)
# 
# AlternateStops_Bischof <- AlternateStops_Bischof %>% group_by(stop.n.c) %>% mutate(km = round(mean(km_mark, 0))) %>% ungroup() %>% dplyr::select(names(onstop_yr_fin %>% dplyr::select(-IncompMig)))
# 
# stop_data <- rbind(onstop_yr_fin %>% dplyr::select(-IncompMig), AlternateStops_Bischof)
# stop_data <- stop_data %>% arrange(id_yr, POSIXct, flag)
# rm(AlternateStops_Bischof)
# proj <- st_crs(onstop_yr_fin)$input
# 
# stop_data <- stop_data %>% st_as_sf(coords = c("xend","yend"), crs = proj)
# 
# save(stop_data, file = "Data_out/Data/Stopover/RDH_UsedAvailStops_14t22_20230428.RData")
# 
# #load data
# load("Data/Stopover/RDH_UsedAvailStops_14t22_20230426.RData")
# head(stop_data)
# stop_data$flag <- as.character(stop_data$flag)
# 
# mapview(stop_data[which(stop_data$id_yr == unique(stop_data$id_yr)[1]),], zcol = "flag")
# 
# #load greenscape calcs
# load("Data/Greenscape/RDH_14t22_greenscapes_20230425.RData")
# head(RDH_14t22_greenscapes)
# 
# #load fidelity to stops
# load("Data/Stopover/FidelityMetrics_20230425.RData")
# 
# #per stop
# stop1 <- IYD_stop_fidelity_list[[1]]
# head(stop1)
# 
# #route level
# mean1 <- IYD_mean_fidelity_list[[1]]
# head(mean1)

setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
load("Data_out/Data/Stopover/RDH_UsedAvailStops_14t22_20230428.RData")

# readRDS("./REnvs/Analysis_04282023.RDS")

# extract fidelity
meanfid1 <- IYD_mean_fidelity_list[[1]]
meanfid2 <- IYD_mean_fidelity_list[[2]]
meanfid3 <- IYD_mean_fidelity_list[[3]]

fid1 <- IYD_stop_fidelity_list[[1]]
fid2 <- IYD_stop_fidelity_list[[2]]
fid3 <- IYD_stop_fidelity_list[[3]]

himeanfid1 <- IYD_HighUse_mean_fidelity_list[[1]]
himeanfid2 <- IYD_HighUse_mean_fidelity_list[[2]]
himeanfid3 <- IYD_HighUse_mean_fidelity_list[[3]]

hifid1 <- IYD_HighUse_stop_fidelity_list[[1]]
hifid2 <- IYD_HighUse_stop_fidelity_list[[2]]
hifid3 <- IYD_HighUse_stop_fidelity_list[[3]]

# Hypotheses

#Greenscape quality and fidelity affect surfing score
## absDFP ~ Duration*IYD + Order*IYD at 1, 2, 3 years

#Stops with greater fidelity are because of longterm stability
## IYD_km ~ sdiNDVI (random = id_yr, km)

#information gathering along route
## IYD ~ km, gam for nonlinear

# how have conditions changed along route
## hist of mean iNDVI


#----------------------------#
# Information hypothesis ####
# route level

head(RDH_14t22_greenscapes)
head(IYD_mean_fidelity_list[[1]])




RDH_14t22_greenscapes = RDH_14t22_greenscapes %>% rename(id_yr = AY_ID, AID = AnimalID)

test <- meanfid1 %>% left_join(RDH_14t22_greenscapes %>% dplyr::select(-c(Year, AID)), by = "id_yr")


#does the greenscape affect fidelity? all stops
mod.GWFidelity.all <- glmmTMB(DFP_120 ~ svSlope*mean_fid_1 + greenUpDur*mean_fid_1 + (1|AID) + (1|curr_Y) + (1|id_yr), data = test, family = poisson(link = "log"), na.action = "na.fail")
plot(DHARMa::simulateResiduals(mod.GWFidelity.all))
summary(mod.GWFidelity.all)

newa <- data.frame(svSlope = seq(min(na.omit(test$svSlope)),max(na.omit(test$svSlope)),length = 100), greenUpDur = seq(min(na.omit(test2$greenUpDur)),max(na.omit(test2$greenUpDur)),length = 100), AID = unique(test2$AID)[10], curr_Y = unique(test2$curr_Y)[10], id_yr = unique(test2$id_yr)[10])
newa$fit <- predict(mod.GWFidelity.all, new, type = "response", se.fit = T)$fit
newa$se <- predict(mod.GWFidelity.all, new, type = "response", se.fit = T)$se.fit
# newa$real <- exp(new$fit)
# newa$se <- exp(1.96*new$se)

test2 <- IYD_HighUse_mean_fidelity_list[[1]] %>% left_join(RDH_14t22_greenscapes %>% dplyr::select(-c(Year, AID)), by = "id_yr")

#does the greenscape affect fidelity, high use only
mod.GWFidelity.hi <- glmmTMB(mean_fid_1 ~ svSlope + greenUpDur + (1|AID) + (1|curr_Y) + (1|id_yr), data = test2, family = Gamma(link = "log"))
plot(DHARMa::simulateResiduals(mod.GWFidelity.hi))
summary(mod.GWFidelity.hi)

new <- data.frame(svSlope = seq(min(na.omit(test2$svSlope)),max(na.omit(test2$svSlope)),length = 100), greenUpDur = seq(min(na.omit(test2$greenUpDur)),max(na.omit(test2$greenUpDur)),length = 100), AID = unique(test2$AID)[10], curr_Y = unique(test2$curr_Y)[10], id_yr = unique(test2$id_yr)[10])
new$fit <- predict(mod.GWFidelity.hi, new, type = "response", se.fit = T)$fit
new$se <- predict(mod.GWFidelity.hi, new, type = "response", se.fit = T)$se.fit
# new$real <- exp(new$fit)
# new$upp <- exp(new$fit + 1.96*new$se)
# new$low <- exp(new$fit - 1.96*new$se)

#plot of high use effect
ggplot(new) + geom_line(aes(x = greenUpDur, y = fit), data = newa, color = "#CD853F", size = 1.2) + geom_line(aes(x = greenUpDur, y = fit), color = "#008B8B", size = 1.2) + geom_line(aes(x = greenUpDur, y = fit + se), linetype = "dashed", color = "#008B8B", size = 1.2) + geom_line(aes(x = greenUpDur, y = fit - se), linetype = "dashed", color = "#008B8B", size = 1.2)   + labs(x = "Green-up Duration (d)", y = "Inter-year Distance (m)", fill = "") + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + theme_classic() #+ coord_cartesian(xlim = c(0,22), ylim = c(0,0.5))

#ggsave(filename = "./Figures/ResponseToGreenscape_GreenUp Order.svg", dpi = 500, width = 15, height = 10, units = "cm")

#--------------------------------#
# information along route ####


head(IYD_HighUse_stop_fidelity_list[[1]])


onstop_yr_fin <- onstop_yr_fin %>% group_by(stop.n.c) %>% mutate(km = round(mean(km_mark), 0))

test <- IYD_stop_fidelity_list[[1]] %>% rename(stop.n.c = stop.n_1) %>% left_join(onstop_yr_fin %>% dplyr::select(-AID), by = c("id_yr", "stop.n.c")) %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

test[which(test$iyd_1 == 0),"iyd_1"] <- 1

#does fidelity get better as they go? all
mod.SampAlongRoute.all <- glmmTMB(iyd_1 ~ km + absDFP + (1|AID) + (1|id_yr), data =test, family = Gamma(link = "log"))
plot(DHARMa::simulateResiduals(mod.SampAlongRoute.all))
summary(mod.SampAlongRoute.all)

newS <- data.frame(km = seq(min(na.omit(test$km)),max(na.omit(test$km)),length = 200), absDFP = seq(min(na.omit(test$absDFP)),max(na.omit(test$absDFP)),length = 200), AID = unique(test$AID)[10], id_yr = unique(test$id_yr)[10])
newS$fit <- predict(mod.SampAlongRoute.all, newS, type = "response", se.fit = T)$fit
newS$se <- predict(mod.SampAlongRoute.all, newS, type = "response", se.fit = T)$se.fit

ggplot(newS) + geom_line(aes(x = km, y = fit), color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = fit + se), linetype = "dashed", color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = fit - se), linetype = "dashed", color = "#008B8B", size = 1.2)

testa <- IYD_HighUse_stop_fidelity_list[[1]] %>% rename(stop.n.c = stop.n_1) %>% left_join(onstop_yr_fin %>% dplyr::select(-AID), by = c("id_yr", "stop.n.c")) %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

testa[which(testa$iyd_1 == 0),"iyd_1"] <- 1
head(testa)

#does fidelity get better as they go? hi use only
testa1 <- testa %>% filter(!is.na(km) | !is.na(absDFP) | !is.na(IRG_scaled)) %>% mutate(fY = as.character(curr_Y), fAID = as.character(AID))


head(testa1)

# table(is.na(testa1$iyd_1))
# table(is.na(testa1$km))
# table(is.na(testa1$absDFP))
# table(is.na(testa1$IRG_scaled))
# table(is.na(testa1$fY))
# table(is.na(testa1$fAID))
# table(is.na(testa1$id_yr))

#

#mod <- gamm(absDFP ~ s(iyd_1, bs="cs", k = 5) + s(km, bs="cs", k = 5), random=list(id_yr=~1, fY=~1, fAID=~1), method = "ML", na.action=na.exclude, data = testa1, family = poisson(link = "log"))
#plot(DHARMa::simulateResiduals(mod$lme))
summary(mod$gam)

mod.SampAlongRoute.hi.gam <- gamm(round(iyd_1,0) ~ s(km, bs = "cs", k = 5), random=list(id_yr=~km|id_yr, fY=~1, fAID=~1), method = "ML", na.action=na.exclude, data = testa1, family = poisson(link = "log")) #convergence issues, but doesnt work with Gamma either
summary(mod.SampAlongRoute.hi.gam$gam)
summary(mod.SampAlongRoute.hi.gam$lme)

newSH <- data.frame(km = seq(min(na.omit(testa$km)),max(na.omit(testa$km)),length = 200), absDFP = seq(min(na.omit(testa$absDFP)),43,length = 200), AID = unique(testa$AID)[10], id_yr = unique(testa$id_yr)[10])
newSH$fit <- predict(mod.SampAlongRoute.hi.gam$gam, newSH,re.form=NA,se.fit = TRUE, type="link")$fit
newSH$se <- predict(mod.SampAlongRoute.hi.gam$gam, newSH, re.form=NA,se.fit = TRUE, type="link")$se.fit

newSH$real <- exp(newSH$fit)
newSH$upp <- exp(newSH$fit + 1.96*newSH$se)
newSH$low <- exp(newSH$fit - 1.96*newSH$se)

ggplot(newSH) + geom_line(aes(x = km, y = real), color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = upp), linetype = "dashed", color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = low), linetype = "dashed", color = "#008B8B", size = 1.2) + theme_classic()


#save.image("./REnvs/Analysis_04282023.RDS")



# mod.SampAlongRoute.hi <- glmmTMB(iyd_1 ~ km + absDFP + (1|AID) + (1|id_yr), data =testa1, family = Gamma(link = "log")) # "log"
# plot(DHARMa::simulateResiduals(mod.SampAlongRoute.hi))
# summary(mod.SampAlongRoute.hi)
# 
# newSH <- data.frame(km = seq(min(na.omit(testa$km)),max(na.omit(testa$km)),length = 200), absDFP = seq(min(na.omit(testa$absDFP)),43,length = 200), AID = unique(testa$AID)[10], id_yr = unique(testa$id_yr)[10])
# newSH$fit <- predict(mod.SampAlongRoute.hi, newSH, type = "response", se.fit = T)$fit
# newSH$se <- predict(mod.SampAlongRoute.hi, newSH, type = "response", se.fit = T)$se.fit
# 
# ggplot(newSH) + geom_line(aes(x = km, y = fit), color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = fit + se), linetype = "dashed", color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = fit - se), linetype = "dashed", color = "#008B8B", size = 1.2)
# 
# ggplot(newSH) + geom_line(aes(x = absDFP, y = fit), color = "#008B8B", size = 1.2) + geom_line(aes(x = absDFP, y = fit + se), linetype = "dashed", color = "#008B8B", size = 1.2) + geom_line(aes(x = absDFP, y = fit - se), linetype = "dashed", color = "#008B8B", size = 1.2)


#--------------------------#
# Cue affects fidelity ####
data1 <- mean1 %>% left_join(RDH_14t22_greenscapes %>% rename(id_yr = AY_ID), by = "id_yr")
head(data1)

mod1 <- glmmTMB(mean_fid_1 ~ scale(greenUpDur) + scale(springScale) + (1|AnimalID) + (1|Year), data = data1, family = poisson())
summary(mod1)

#--------------------------#
# Mismatch affects fidelity ####
data2 <- stop_data %>% filter(flag == "1") %>% left_join(stop1, by = "id_yr", relationship = "many-to-many")
head(data2)
data2 <- data2 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

mod2 <- glmmTMB(iyd_1 ~ absDFP + DFP + (1|AID.x) + (1|Year), data = data2, family = poisson())


#-----------------------#
# Signal of Stop ####

stop_data <- stop_data %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP)) %>% group_by(stop.n.c, flag) %>% mutate(varPeak = sd(MaxIRGday), varSpring = sd(SpringLength)) %>% ungroup()

str(stop_data)

TMBStruc = glmmTMB(flag ~ absDFP + varPeak + varSpring +   
                (1|stop.n.c)+(0 + absDFP| id_yr),family=poisson,
               doFit=F,data=stop_data)

# # Fix the standard deviation of the first random term, which is the (1|strata) component
# in the above model equation
TMBStruc$parameters$theta[1] = log(1e3)  # see Muff et al. for details
TMBStruc$mapArg = list(theta=factor(c(NA,1)))   # I have 3 random effects, thus the theta is NA,1,2,3. If there is 4 ranef, then c(NA,1:4)


# Fit the model
mTMB = glmmTMB:::fitTMB(TMBStruc)
summary(mTMB)
AIC(mTMB)
x <- as.data.frame(confint(mTMB))
names(x) <- c("low", "upp", "mean")

ggplot(data = x[c(2:4),]) + geom_bar(aes(x = mean, y = row.names(x[c(2:4),])), stat= "identity", fill = "#008B8B", color = "grey40") + geom_linerange(aes(x = mean, xmin = low, xmax = upp,y = row.names(x[c(2:4),])), size = 1.1)

#stop_data$predval <- predict(mTMB, type="risk")

new <- data.frame(absDFP = seq(min(stop_data$absDFP), max(stop_data$absDFP), length = 100), varPeak = seq(min(stop_data$varPeak), max(stop_data$varPeak), length = 100),varSpring = seq(min(stop_data$varSpring), max(stop_data$varSpring), length = 100), stop.n.c = unique(stop_data$stop.n.c)[130], id_yr = unique(stop_data$id_yr)[30])
new$fit <- predict(mTMB, new, se.fit = T, type = "response")$fit
new$se.fit <- predict(mTMB, new, se.fit = T, type = "response")$se.fit



#--------------------------#
# Age/Cohort affects fidelity ####

#--------------------------#
# Repeatable for animals ####


#--------------------------#
# Landscape affects fidelity ####






#--------------------------#
# Repeatability of Stopover Behavior ####
load("Data/Stopover/FidelityMetrics_20230425.RData")

stop1 <- IYD_stop_fidelity_list[[1]]
head(stop1)

stop1$iyd_1 <- round(stop1$iyd_1,0)

#need to calculate the quartiles to round the kms


rep1 <- rpt(iyd_1 ~  (1|AID), grname = c("AID"), data = stop1, datatype = "Poisson", nboot = 100, 
            npermut = 0)
print(rep1)

#--------------------------#
# Investigate need for PCA ####
# 
# tmp <- stop_data %>% st_drop_geometry() %>% as.data.frame()
# 
# 
# tmp[c("MaxNDVIDay","MaxIRGday", "SpringStartDay","SpringEndDay","IntegratedNDVI","MaxBrownDownDate",
#       "NDVI_scaled","IRG_scaled","SE_springDLCfit", "SpringLength","sumIRG","SpringScale")] <- sapply(tmp[c("MaxNDVIDay","MaxIRGday", "SpringStartDay","SpringEndDay","IntegratedNDVI","MaxBrownDownDate",
#                                                                                                             "NDVI_scaled","IRG_scaled","SE_springDLCfit", "SpringLength","sumIRG","SpringScale")], as.numeric)
# 
# tmp <- tmp[c("MaxNDVIDay","MaxIRGday", "SpringStartDay","SpringEndDay","IntegratedNDVI","MaxBrownDownDate",
#              "NDVI_scaled","IRG_scaled","SE_springDLCfit", "SpringLength","sumIRG","SpringScale")] %>% st_drop_geometry()
# 
# #covariance matrix
# round(cor(tmp),2)
# 
# #ifelse(round(cor(tmp),2) > .7,round(cor(tmp),2),"NA")#check the covariance of 