# Script for modelling

## Hypotheses

#Greenscape quality and fidelity affect surfing score
## absDFP ~ Duration*IYD + Order*IYD at 1, 2, 3 years

#Stops with greater fidelity are because of longterm stability
## IYD_km ~ sdiNDVI (random = id_yr, km)

#information gathering along route
## IYD ~ km, gam for nonlinear

# how have conditions changed along route
## hist of mean iNDVI

# how does mismatch change when using previously visited or not visited stop
## DaysDiff ~ Status

load("REnvs/ModellingEnvironment.rds")

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
library(mgcViz)
library(visreg)
library(tidymv)
remotes::install_github("stefanocoretta/tidygam@v0.1.0")
library(tidygam)
library(StatisticalModels)


#--------------------------#
# prep datasets ####

setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
load("Data_out/Data/Stopover/RDH_StopoverBundle_14t22_20230425.RData")
load("Data_out/Data/Stopover/FidelityMetrics_20230428.RData")
load("Data_out/Data/Stopover/RDH_AlternateStops_Bischof_20230428.RData")
load("Data_out/Data/Greenscape/RDH_14t22_greenscapes_20230425.RData")
load("Data_out/Data/Stopover/RDH_UsedAvailStops_14t22_20230428.RData")
load("C:/Users/lwilde2/Desktop/RDH Database/Processed Data/RDH_AllMigrations_Bischof_2014t2022_20230425.RData")

# readRDS("./REnvs/Analysis_04282023.RDS")

# extract fidelity from lists
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

x <- himeanfid1 %>% filter(AID %in% himeanfid2$AID | AID %in% himeanfid3$AID) %>% dplyr::select(AID)
x1 <- himeanfid2 %>% filter(AID %in% himeanfid3$AID) %>% select(AID) %>% dplyr::select(AID)


#greenscape data
RDH_14t22_greenscapes
#stops
onstop_yr_fin
#full data
data1
#fidelity
himeanfid1
himeanfid2
himeanfid3

fid1
fid2
fid3

#adaptability of fidelity?
#absDFP ~ greenscape * IYD
#DFP ~ greenscape * IYD

#origins of fidelity?
#landscape
#response to mismatch


ggplot(fid1) + theme(legend.position = "none") + coord_cartesian(xlim = c(0,220), ylim = c(0,100000)) + geom_point(aes(x = km_mark_1, y = iyd_1, group = AID, color = AID),  alpha = .2) + geom_smooth(aes(x = km_mark_1, y = iyd_1)) #geom_line(aes(x = km_mark_1, y = iyd_1, group = AID, color = AID)) +


# First hypothesis, accumulated iyd should hurt ability to surf - decay of information?

t1 <- onstop_yr_fin %>% left_join(fid1 %>% dplyr::select(iyd_1, stop.n_1, km_mark_1, id_yr) %>% rename(stop.n.c = stop.n_1), by = c("stop.n.c", "id_yr"))
t1 <- t1 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

t2 <- onstop_yr_fin %>% left_join(fid2 %>% dplyr::select(iyd_2, stop.n_2, km_mark_2, id_yr) %>% rename(stop.n.c = stop.n_2), by = c("stop.n.c", "id_yr"))
t2 <- t2 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

t3 <- onstop_yr_fin %>% left_join(fid3 %>% dplyr::select(iyd_3, stop.n_3, km_mark_3, id_yr) %>% rename(stop.n.c = stop.n_3), by = c("stop.n.c", "id_yr"))
t3 <- t3 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP)) 

#cumulative iyd of each point
t1 <- t1 %>% group_by(id_yr) %>% mutate(ciyd = cumsum(iyd_1))
t2 <- t2 %>% group_by(id_yr) %>% mutate(ciyd = cumsum(iyd_2))
t3 <- t3 %>% group_by(id_yr) %>% mutate(ciyd = cumsum(iyd_3))

#lmer of relationship
m1 <- glmer(absDFP ~ scale(ciyd) + (1|id_yr) + (1|AID) + (1|Year), family = poisson(link = "log"), data = t1)
summary(m1)

m2 <- glmer(absDFP ~ scale(ciyd) + (1|id_yr) + (1|AID) + (1|Year), family = poisson(link = "log"), data = t2)
summary(m2)

m3 <- glmer(absDFP ~ scale(ciyd) + (1|id_yr) + (1|AID) + (1|Year), family = poisson(link = "log"), data = t3)
summary(m3)

plot(DHARMa::simulateResiduals(m1))
plot(DHARMa::simulateResiduals(m2))
plot(DHARMa::simulateResiduals(m3))

#rely on trigamma for poisson, needs a log link function
r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m3)


newdata <- data.frame(absDFP = seq(min(na.omit(t1$absDFP)), max(na.omit(t1$absDFP)), length = 50), ciyd = seq(min(na.omit(t1$ciyd)), max(na.omit(t1$ciyd)), length = 50), id_yr = "108_2018", AID = "108", Year = 2016)

p1 <- PredictGLMER(m1, newdata, se.fit = T, seMultiplier = 1.96)
p2 <- PredictGLMER(m2, newdata, se.fit = T, seMultiplier = 1.96)
p3 <- PredictGLMER(m3, newdata, se.fit = T, seMultiplier = 1.96)

ggplot() + geom_ribbon(data = cbind(p1, newdata), aes(x = ciyd, ymin = yminus, ymax = yplus), fill = "#008B8B", alpha = .75) + geom_line(data = cbind(p1, newdata), aes(x = ciyd, y=y), color = "black") + geom_ribbon(data = cbind(p2, newdata), aes(x = ciyd, ymin = yminus, ymax = yplus), fill = "#008B8B", alpha = .3) + geom_line(data = cbind(p2, newdata), aes(x = ciyd, y=y), color = "black") + geom_ribbon(data = cbind(p3, newdata), aes(x = ciyd, ymin = yminus, ymax = yplus), fill = "#008B8B", alpha = .1) + geom_line(data = cbind(p3, newdata), aes(x = ciyd, y=y), color = "black") + theme_classic() + labs(y = "|Days from Peak Green-up|", x = bquote('cumulative Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"))

#ggsave(filename = "Figures/ResponseToCumulativeIYD.svg", width = 24, height = 24, units = "cm", dpi = 600)
#ggsave(filename = "Figures/ResponseToCumulativeIYD.jpg", width = 24, height = 24, units = "cm", dpi = 600)


# hypothesis two, dealing with different conditions
t1
green <- RDH_14t22_greenscapes %>% rename(id_yr = AY_ID, AID = AnimalID)
t1g <- t1 %>% left_join(green, by = c("id_yr")) %>% ungroup() %>% mutate(ciyds = scale(ciyd))

t1g <- t1g %>% mutate(GUDs = scale(greenUpDur), SPSs = scale(springScale), SVSs = scale(svSlope), GUDc = ifelse(GUDs <= 0, "short", "long"), SPSc = ifelse(SPSs <= 0,  "fast", "slow" ), SVSc = ifelse(SVSs <= 0,"rand", "ord")) %>% rename(Year = Year.x, AID = AID.x)


mod.IYDbyGUD <- mgcv::gamm(absDFP ~ s(ciyds, bs = "cs", k = 5, by = factor(t1g$GUDc)) +  SVSs, method = "REML", random = list(id_yr=~ciyds|id_yr, Year = ~1, AID=~1), family = poisson(link = "log"), data = t1g, niterPQL = 20) #SPSs +

#mod.IYDbySPS <- mgcv::gamm(absDFP ~ s(ciyds, bs = "cs", k = 5, by = factor(t1g$SPSc)) + GUDs + SVSs, method = "REML", random = list(id_yr=~ciyds|id_yr, Year = ~1, AID=~1), family = poisson(link = "log"), data = t1g, niterPQL = 20)

mod.IYDbySVS <- mgcv::gamm(absDFP ~ s(ciyds, bs = "cs", k = 5, by = factor(SVSc)) + GUDs, method = "REML", random = list(id_yr=~ciyds|id_yr, Year = ~1, AID=~1), family = poisson(link = "log"), data = t1g, niterPQL = 20) #SPSs +

summary(mod.IYDbyGUD$gam)
summary(mod.IYDbySVS$gam)

summary(mod.IYDbyGUD$lme)
summary(mod.IYDbySVS$lme)

par(mfrow = c(2,2))
plot.gam(mod.IYDbyGUD$gam)
plot.gam(mod.IYDbySVS$gam)

dev.off()

newdata <- data.frame(expand.grid(ciyds = seq(-3,3,length = 50), SVSc = unique(t1g$SVSc), GUDs = seq(-3,3,length = 50), id_yr = t1g$id_yr[4], AID = t1g$AID[4], Year = t1g$Year[4]))
predict(mod.IYDbySVS, newdata, type = "response")




# NOT RUN BELOW #






#------------------------#
# adaptability of fidelity ####

#prepare the dataset
# two tests
# do you have lower mismatch when closer to last year
# do you have lower if you the last stop was last years

# cumulative IYD compared to cumulative absDFP

# build dataset for each year

t1 <- onhigh_yr_fin %>% left_join(hifid1 %>% dplyr::select(iyd_1, stop.n_1, km_mark_1, id_yr) %>% rename(stop.n.c = stop.n_1), by = c("stop.n.c", "id_yr"))
t1 <- t1 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

t1 <- t1 %>% filter(!is.na(iyd_1))

t2 <- t1 %>% group_by(id_yr, Date) %>% sample_n(1) %>% ungroup() %>% group_by(id_yr) %>% summarise(cDFP = sum(DFP), cabsDFP = sum(absDFP), ciyd = sum(na.omit(iyd_1)), id_yr = unique(id_yr)) # here we need control for different lengths of the season

t2 <- t2 %>% left_join(Spring.summary %>% dplyr::select(id_yr, length), by = "id_yr")

t2 <- t2 %>% st_drop_geometry() %>% left_join(RDH_14t22_greenscapes %>% rename(id_yr = AY_ID, AID = AnimalID), by = "id_yr")

head(t2)
cor(t2[,c("length", "cDFP", "cabsDFP", "ciyd")])

t2s <- t2 %>% mutate(ciyd = scale(ciyd), greenUpDur = scale(greenUpDur), svSlope = scale(svSlope), springScale = scale(springScale)) #, cabsDFP = cabsDFP/length, cDFP = cDFP/length

# t2s <- t2s %>% mutate(GUDc = ifelse(greenUpDur < -.5, "short", ifelse(greenUpDur >= -.5 & greenUpDur <= .5, "mod", "long")), RATEc = ifelse(svSlope < -.5, "slow", ifelse(svSlope >= -.5 & svSlope <= .5, "mod", "fast")), ORDc = ifelse(springScale < -.5, "rand", ifelse(springScale >= -.5 & springScale <= .5, "mod", "Ord")))

#shouldnt need to correct for length, non-proportional to the data calculated




mod.adap_linear <- mgcv::gamm(cabsDFP ~ ciyd*greenUpDur + ciyd*svSlope + ciyd*springScale, random = list(id_yr=~1, Year=~1, AID=~1), method = "REML", family = poisson(link = "log"), data = t2s)
summary(mod.adap_linear$lme)
summary(mod.adap_linear$gam)


mod.absDFP_nl <- mgcv::gamm(cabsDFP ~ s(ciyd, bs = "cs", k = 5, by = GUDc) + svSlope + springScale, method = "REML", random = list(id_yr=~ciyd|id_yr, Year=~1, AID=~1), family = poisson(link = "log"), data = t2s)


mod.absDFP_nl <- mgcv::gamm(cabsDFP ~ s(ciyd) + greenUpDur + svSlope + springScale, method = "REML", random = list(id_yr=~ciyd|id_yr, Year=~1, AID=~1), family = poisson(link = "log"), data = t2s) #+ s(greenUpDur, svSlope , springScale), random = list(id_yr=~ciyd|id_yr, Year=~1, AID=~1)
summary(mod.absDFP_nl$gam)
summary(mod.absDFP_nl$lme)

plot.gam(mod.absDFP_nl$gam, shade = T, shade.col = "#008B8B", col = "black", rug = F, lwd = 2.5, ylab = "cumulate |DFP|", xlab = "cumulative Fidelity Index", residuals= T)

mod.DFP_nl <- mgcv::gamm(cDFP ~ s(ciyd, bs = "cs", k = 5) + (greenUpDur + svSlope + springScale), random = list(id_yr=~ciyd|id_yr, Year=~1, AID=~1), family = "gaussian", method = "REML", data = t2s)
summary(mod.DFP_nl$gam)
summary(mod.DFP_nl$lme)














#predict from models
absP <- predict_gam(mod.absDFP_nl, type = "terms")

absP %>% ggplot(aes(x = ciyd, y = fit)) + geom_smooth_ci()


newabs <- data.frame(ciyd = seq(min(t2$ciyd), max(t2$ciyd), by = 100), greenUpDur = seq(min(na.omit(t2$greenUpDur)), max(na.omit(t2$greenUpDur)), by = 100), svSlope = seq(min(na.omit(t2$svSlope)), max(na.omit(t2$svSlope)), by = 100), springScale = seq(min(na.omit(t2$springScale)), max(na.omit(t2$springScale)), by = 100), id_yr = unique(t2$id_yr)[3], Year = unique(t2$Year)[3], AID = unique(t2$AID)[3]) #expand.grid()
newabs$fit <- predict(mod.absDFP_nl$gam, newabs,re.form=F,se.fit = TRUE, terms = "s(ciyd)")$fit
newabs$se <- predict(mod.absDFP_nl$gam, newabs, re.form=F,se.fit = TRUE, terms = "s(ciyd)")$se.fit

ggplot(newabs) + geom_line(aes(x = ciyd, y = fit))

# 
# newabs$real <- exp(newabs$fit)
# newabs$upp <- exp(newabs$fit + 1.96*newabs$se)
# newabs$low <- exp(newabs$fit - 1.96*newabs$se)


mod.test <- mgcv::gamm(cabsDFP ~ s(greenUpDur, by = ciyd, bs = "cs", k = 5) + s(springScale, by = ciyd, bs = "cs", k = 5) + s(svSlope, by = ciyd, bs = "cs", k = 5) , random = list(id_yr=~1, Year=~1, AID=~1), method = "REML", family = poisson(link = "log"), data = t2s)











#--------------------------------#
# information along route ####

#all stops
onstop_yr_fin <- onstop_yr_fin %>% group_by(stop.n.c) %>% mutate(km = round(mean(km_mark), 0))
test <- fid1 %>% rename(stop.n.c = stop.n_1) %>% left_join(onstop_yr_fin %>% dplyr::select(-AID), by = c("id_yr", "stop.n.c")) %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))
test[which(test$iyd_1 == 0),"iyd_1"] <- 1

#high use only
testhi <- hifid1 %>% rename(stop.n.c = stop.n_1) %>% left_join(onstop_yr_fin %>% dplyr::select(-AID), by = c("id_yr", "stop.n.c")) %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))
testhi[which(testhi$iyd_1 == 0),"iyd_1"] <- 1

testhi$iyd_1 <- round(testhi$iyd_1, 0)
test$iyd_1 <- round(test$iyd_1, 0)

testhi_simple <- testhi %>% group_by(stop.n.c) %>% summarise(curr_y = unique(curr_Y), AID = unique(AID), id_yr = unique(id_yr), iyd_1 = unique(iyd_1), MaxIRGday = mean(MaxIRGday), absDFP = mean(absDFP), km = round(mean(km_mark_1)))


head(test)

#does fidelity get better as they go? all
mod.SampAlongRoute.All.gam <- gamm(iyd_1 ~ s(km, bs = "cs", k = 5), random=list(id_yr=~km|id_yr, Year=~1, AID=~1), method = "REML", na.action=na.exclude, data = test, family = poisson(link = "log")) 
summary(mod.SampAlongRoute.All.gam$gam)
summary(mod.SampAlongRoute.All.gam$lme)

newS <- data.frame(km = seq(min(na.omit(test$km)),200,length = 200), AID = unique(test$AID)[10], id_yr = unique(test$id_yr)[10])
newS$fit <- predict(mod.SampAlongRoute.All.gam$gam, newS,re.form=NA,se.fit = TRUE, type="link")$fit
newS$se <- predict(mod.SampAlongRoute.All.gam$gam, newS, re.form=NA,se.fit = TRUE, type="link")$se.fit

newS$real <- exp(newS$fit)
newS$upp <- exp(newS$fit + 1.96*newS$se)
newS$low <- exp(newS$fit - 1.96*newS$se)

ggplot(newS) + geom_line(aes(x = km, y = real), color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = upp), linetype = "dashed", color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = low), linetype = "dashed", color = "#008B8B", size = 1.2)



mod.SampAlongRoute.hi.gam <- gamm(iyd_1 ~ s(km, bs = "cs", k = 5), random=list(id_yr=~km|id_yr, Year=~1, AID=~1), method = "REML", na.action=na.exclude, data = testhi, family = poisson(link = "log")) #convergence issues, but doesnt work with Gamma either
summary(mod.SampAlongRoute.hi.gam$gam)
summary(mod.SampAlongRoute.hi.gam$lme)

newSH <- data.frame(km = seq(min(na.omit(testhi$km)),220,length = 200), absDFP = seq(min(na.omit(testhi$absDFP)),43,length = 200), AID = unique(testhi$AID)[10], id_yr = unique(testhi$id_yr)[10])
newSH$fit <- predict(mod.SampAlongRoute.hi.gam$gam, newSH,re.form=NA,se.fit = TRUE, type="link")$fit
newSH$se <- predict(mod.SampAlongRoute.hi.gam$gam, newSH, re.form=NA,se.fit = TRUE, type="link")$se.fit

newSH$real <- exp(newSH$fit)
newSH$upp <- exp(newSH$fit + 1.96*newSH$se)
newSH$low <- exp(newSH$fit - 1.96*newSH$se)

ggplot(newSH) + geom_line(aes(x = km, y = real), color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = upp), linetype = "dashed", color = "#008B8B", size = 1.2) + geom_line(aes(x = km, y = low), linetype = "dashed", color = "#008B8B", size = 1.2) + theme_classic()

ggsave("Figures/FidelityByDistance_rough.png")

#how does this compare to absDFP
data.daily <- data1 %>% group_by(id_yr, Date) %>% sample_n(1)

data.daily <- data.daily %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

mod.gws <- gamm(absDFP ~ s(km_mark, bs = "cs", k = 5), random=list(id_yr=~km_mark|id_yr, Year=~1, AID=~1), method = "REML", na.action=na.exclude, data = data.daily, family = poisson(link = "log"))
summary(mod.gws$gam)
summary(mod.gws$lme)

newGWS <- data.frame(km_mark = seq(min(na.omit(data.daily$km_mark)),220,length = 200), AID = unique(data.daily$AID)[10], id_yr = unique(data.daily$id_yr)[10], Year = unique(data.daily$Year)[3])
newGWS$fit <- predict(mod.gws$gam, newGWS,re.form=NA,se.fit = TRUE, type="link")$fit
newGWS$se <- predict(mod.gws$gam, newGWS, re.form=NA,se.fit = TRUE, type="link")$se.fit

newGWS$real <- exp(newGWS$fit)
newGWS$upp <- exp(newGWS$fit + 1.96*newGWS$se)
newGWS$low <- exp(newGWS$fit - 1.96*newGWS$se)

ggplot(newGWS) + geom_line(aes(x = km_mark, y = real), color = "#008B8B", size = 1.2) + geom_line(aes(x = km_mark, y = upp), linetype = "dashed", color = "#008B8B", size = 1.2) + geom_line(aes(x = km_mark, y = low), linetype = "dashed", color = "#008B8B", size = 1.2) + theme_classic()

ggsave("Figures/MismatchDistance_rough.png")

#-------------------------------#
# ####

# what could be driving fidelity

#IYD ~ km + mismatch + greenscape + 



mod.DriveGlobal.Hi <- gam::gam(iyd_1 ~ km_mark_1 + absDFP + IntegratedNDVI + (1|id_yr) + (1|), family = poisson(link = "log"), data = testhi, na.action = "na.fail")











#---------------------------------#
#compare used and alternate####

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



