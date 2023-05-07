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
library(patchwork)
library(gammit)

#--------------------------#
# prep datasets ####

setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
load("Data_out/Data/Stopover/RDH_StopoverBundle_14t22_20230425.RData")
load("Data_out/Data/Stopover/FidelityMetrics_20230503.RData")
#load("Data_out/Data/Stopover/RDH_AlternateStops_Bischof_20230428.RData")
load("Data_out/Data/Greenscape/RDH_14t22_greenscapes_20230425.RData")
#load("Data_out/Data/Stopover/RDH_UsedAvailStops_14t22_20230428.RData")
load("C:/Users/lwilde2/Desktop/RDH Database/Processed Data/RDH_AllMigrations_Bischof_2014t2022_20230425.RData")

# readRDS("./REnvs/Analysis_04282023.RDS")

# extract fidelity from lists
# meanfid1 <- IYD_mean_fidelity_list[[1]]
# meanfid2 <- IYD_mean_fidelity_list[[2]]
# meanfid3 <- IYD_mean_fidelity_list[[3]]

fid1 <- IYD_stop_fidelity_list[[1]]
fid2 <- IYD_stop_fidelity_list[[2]]
fid3 <- IYD_stop_fidelity_list[[3]]

# himeanfid1 <- IYD_HighUse_mean_fidelity_list[[1]]
# himeanfid2 <- IYD_HighUse_mean_fidelity_list[[2]]
# himeanfid3 <- IYD_HighUse_mean_fidelity_list[[3]]

hifid1 <- IYD_HighUse_stop_fidelity_list[[1]]
hifid2 <- IYD_HighUse_stop_fidelity_list[[2]]
hifid3 <- IYD_HighUse_stop_fidelity_list[[3]]

# x <- himeanfid1 %>% filter(AID %in% himeanfid2$AID | AID %in% himeanfid3$AID) %>% dplyr::select(AID)
# x1 <- himeanfid2 %>% filter(AID %in% himeanfid3$AID) %>% select(AID) %>% dplyr::select(AID)


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

t1 <- onhigh_yr_fin %>% left_join(hifid1 %>% dplyr::select(iyd_1, stop.n_1, km_mark_1, id_yr) %>% rename(stop.n.c = stop.n_1), by = c("stop.n.c", "id_yr")) %>% arrange(id_yr, POSIXct)
t1 <- t1 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

t2 <- onhigh_yr_fin %>% left_join(hifid2 %>% dplyr::select(iyd_2, stop.n_2, km_mark_2, id_yr) %>% rename(stop.n.c = stop.n_2), by = c("stop.n.c", "id_yr")) %>% arrange(id_yr, POSIXct)
t2 <- t2 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

t3 <- onhigh_yr_fin %>% left_join(hifid3 %>% dplyr::select(iyd_3, stop.n_3, km_mark_3, id_yr) %>% rename(stop.n.c = stop.n_3), by = c("stop.n.c", "id_yr")) %>% arrange(id_yr, POSIXct)
t3 <- t3 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP)) 

#cumulative iyd of each point
# t1 <- t1 %>% group_by(stop.n.c) %>% mutate(iydsimp = ifelse(duplicated(iyd_1), 0, iyd_1)) %>% ungroup() %>% group_by(id_yr) %>% mutate(ciyd = cumsum(iydsimp)) %>% ungroup() %>% group_by(id_yr) %>% slice(which(row_number() %% 3 == 1))
# t2 <- t2 %>% group_by(stop.n.c) %>% mutate(iydsimp = ifelse(duplicated(iyd_2), 0, iyd_2)) %>% ungroup() %>% group_by(id_yr) %>% mutate(ciyd = cumsum(iydsimp)) %>% group_by(id_yr) %>% slice(which(row_number() %% 3 == 1))
# t3 <- t3 %>% group_by(stop.n.c) %>% mutate(iydsimp = ifelse(duplicated(iyd_3), 0, iyd_3)) %>% ungroup() %>% group_by(id_yr) %>% mutate(ciyd = cumsum(iydsimp)) %>% group_by(id_yr) %>% slice(which(row_number() %% 3 == 1))


#confirmed no relatinoship between
ggplot(t1 %>% group_by(stop.n.c) %>% summarise(avDFP = mean(absDFP), TimeStopped = unique(TimeStopped))) + geom_point(aes(x = TimeStopped, y = avDFP)) + geom_smooth(aes(x = TimeStopped, y = avDFP), method = "loess")

#need to model this, with random effects

#averaging to each stop
t1_sum <- t1 %>% group_by(stop.n.c) %>% summarise(absDFP = round(mean(absDFP),0), DFP = round(mean(DFP),0), TimeStopped = unique(TimeStopped), id_yr = unique(id_yr), AID = unique(AID), Year = unique(Year), iyd = unique(iyd_1), km = mean(km_mark))
t1_sum <- t1_sum %>% group_by(id_yr) %>% mutate(ciyd = cumsum(iyd))

t2_sum <- t2 %>% group_by(stop.n.c) %>% summarise(absDFP = round(mean(absDFP),0), DFP = round(mean(DFP),0), TimeStopped = unique(TimeStopped), id_yr = unique(id_yr), AID = unique(AID), Year = unique(Year), iyd = unique(iyd_2), km = mean(km_mark))
t2_sum <- t2_sum %>% group_by(id_yr) %>% mutate(ciyd = cumsum(iyd))

t3_sum <- t3 %>% group_by(stop.n.c) %>% summarise(absDFP = round(mean(absDFP),0), DFP = round(mean(DFP),0), TimeStopped = unique(TimeStopped), id_yr = unique(id_yr), AID = unique(AID), Year = unique(Year), iyd = unique(iyd_3), km = mean(km_mark))
t3_sum <- t3_sum %>% group_by(id_yr) %>% mutate(ciyd = cumsum(iyd))

#lmer of relationship
#average per stop, and what about including timestopped?

m1 <- glmer(absDFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), family = poisson(link = "log"), data = t1_sum ) #%>% filter(iyd < 20000)
summary(m1)

m2 <- glmer(absDFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), family = poisson(link = "log"), data = t2_sum) #%>% filter(iyd < 20000)
summary(m2)

m3 <- glmer(absDFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), family = poisson(link = "log"), data = t3_sum) #%>% filter(iyd < 20000)
summary(m3)

plot(DHARMa::simulateResiduals(m1))
plot(DHARMa::simulateResiduals(m2))
plot(DHARMa::simulateResiduals(m3))

#rely on trigamma for poisson, needs a log link function
r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m3)


newdata <- data.frame(absDFP = seq(min(na.omit(t1_sum$absDFP)), max(na.omit(t1_sum$absDFP)), length = 50), iyd = seq(min(na.omit(t1_sum$iyd)), max(na.omit(t1_sum$iyd)), length = 50), id_yr = "108_2018", AID = "108", Year = 2016)

 p1 <- PredictGLMER(m1, newdata, se.fit = T, seMultiplier = 1.96)
 p2 <- PredictGLMER(m2, newdata, se.fit = T, seMultiplier = 1.96)
 p3 <- PredictGLMER(m3, newdata, se.fit = T, seMultiplier = 1.96)




plot1 <- ggplot() + geom_ribbon(data = cbind(p1, newdata), aes(x = iyd, ymin = exp(yminus), ymax = exp(yplus)), fill = "#008B8B", alpha = .75) + geom_line(data = cbind(p1, newdata), aes(x = iyd, y=exp(y)), color = "black") + theme_classic() + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + coord_cartesian(xlim = c(0,50000), ylim = c(20,40)) 

plot2 <- ggplot() + geom_ribbon(data = cbind(p2, newdata), aes(x = iyd, ymin = exp(yminus), ymax = exp(yplus)), fill = "#008B8B", alpha = .3) + geom_line(data = cbind(p2, newdata), aes(x = iyd, y=exp(y)), color = "black") + theme_classic() + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"))+ coord_cartesian(xlim = c(0,50000), ylim = c(20,40))

plot3<- ggplot() + geom_ribbon(data = cbind(p3, newdata), aes(x = iyd, ymin = exp(yminus), ymax = exp(yplus)), fill = "#008B8B", alpha = .1) + geom_line(data = cbind(p3, newdata), aes(x = iyd, y=exp(y)), color = "black") + theme_classic() + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"))+ coord_cartesian(xlim = c(0,50000), ylim = c(20,40))

plot1 + plot2 + plot3


#ggsave(filename = "Figures/ResponseToCumulativeIYD.svg", width = 24, height = 24, units = "cm", dpi = 600)
ggsave(filename = "Figures/ResponseToIYD_high_20230503.jpg", width = 34, height = 24, units = "cm", dpi = 600)



#now DFP
m1 <- glmer(DFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), data = t1_sum ) #%>% filter(iyd < 20000)
summary(m1)

m2 <- glmer(DFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), data = t2_sum) #%>% filter(iyd < 20000)
summary(m2)

m3 <- glmer(DFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), data = t3_sum) #%>% filter(iyd < 20000)
summary(m3)


r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m3)

# test relationship between timestopped and absDFP


newdata <- data.frame(DFP = seq(min(na.omit(t1_sum$DFP)), max(na.omit(t1_sum$DFP)), length = 50), iyd = seq(min(na.omit(t1_sum$iyd)), max(na.omit(t1_sum$iyd)), length = 50), id_yr = "108_2018", AID = "108", Year = 2016)

p1 <- PredictGLMER(m1, newdata, se.fit = T, seMultiplier = 1.96)
p2 <- PredictGLMER(m2, newdata, se.fit = T, seMultiplier = 1.96)
p3 <- PredictGLMER(m3, newdata, se.fit = T, seMultiplier = 1.96)




plot1 <- ggplot() + geom_ribbon(data = cbind(p1, newdata), aes(x = iyd, ymin = yminus, ymax = yplus), fill = "#008B8B", alpha = .75) + geom_line(data = cbind(p1, newdata), aes(x = iyd, y=y), color = "black") + theme_classic() + labs(y = "Days from Peak Green-up", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + coord_cartesian(xlim = c(0,50000), ylim = c(15,45)) 

plot2 <- ggplot() + geom_ribbon(data = cbind(p2, newdata), aes(x = iyd, ymin = (yminus), ymax = (yplus)), fill = "#008B8B", alpha = .3) + geom_line(data = cbind(p2, newdata), aes(x = iyd, y=(y)), color = "black") + theme_classic() + labs(y = "Days from Peak Green-up", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"))+ coord_cartesian(xlim = c(0,50000), ylim = c(15,45))


plot1 + plot2 + plot3

ggsave(filename = "Figures/ResponseToIYD_high_DFP_20230503.jpg", width = 34, height = 24, units = "cm", dpi = 600)


#---------------------------------#
# summary stats ####

summary(t1_sum$iyd/1000)
summary(t2_sum$iyd/1000)
summary(t3_sum$iyd/1000)



ggplot() + geom_density(data = t1_sum, aes(x = iyd/1000), bounds = c(0, 250), fill = "blue", alpha = .2) + geom_density(data = t2_sum, aes(x = iyd/1000), bounds = c(0, 250), fill = "red", alpha = .2) + geom_density(data = t3_sum, aes(x = iyd/1000), bounds = c(0, 250), fill = "green", alpha = .2) + coord_cartesian(xlim = c(-2,20)) + scale_x_continuous(limits = c(0,22))


#---------------------------------#
# start the gams ####

mod.absDFP1 <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped, data = t1_sum %>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))
mod.absDFP2 <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped, data = t2_sum%>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))
mod.absDFP3 <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped, data = t3_sum%>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))


newdata = data.frame(km = mean(t1_sum$km), iyd = seq(1,30000,length = 100), TimeStopped = mean(t1_sum$TimeStopped), AID = "108", Year = 2017, id_yr = "108_2017")
pred1 <- gammit::predict_gamm(mod.absDFP1$gam , newdata = newdata, se = T)
pred2 <- gammit::predict_gamm(mod.absDFP2$gam , newdata = newdata, se = T)
pred3 <- gammit::predict_gamm(mod.absDFP3$gam , newdata = newdata, se = T)

newdata1 <- cbind(newdata, pred1)
newdata1$Year <- 1

newdata2 <- cbind(newdata, pred2)
newdata2$Year <- 2

newdata3 <- cbind(newdata, pred3)
newdata3$Year <- 3

new <- rbind(newdata1, newdata2, newdata3)

ggplot(new) + geom_ribbon(aes(x = iyd/1000, ymax = exp(prediction) + exp(se), ymin = exp(prediction) - exp(se), group = factor(Year)), fill = "#008B8B", alpha = .33) + geom_line(aes(x = iyd/1000, y = exp(prediction), linetype = factor(Year)), size = .8) + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance (km)'~( Fidelity ^-1)), linetype = "") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + scale_linetype_manual(values = c("solid", "dashed", "dotted"), label = c("t-1", "t-2", "t-3")) + coord_cartesian(ylim = c(15,45), xlim = c(0,30)) + scale_y_continuous(limits= c(14.8,50.2), breaks = seq(15,50,5)) + scale_x_continuous(limits= c(0,30), breaks = seq(0,30,5))
 
ggsave(filename = "Figures/ResponseToIYD_gam_absDFP_trunc_20230504.jpg", width = 30, height = 30, units = "cm", dpi = 600)








mod.DFP1 <- mgcv::gamm(DFP ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped, data = t1_sum, method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1))
mod.DFP2 <- mgcv::gamm(DFP ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped, data = t2_sum, method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1))
mod.DFP3 <- mgcv::gamm(DFP ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped, data = t3_sum, method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1))


newdata = data.frame(km = mean(t1_sum$km), iyd = seq(1,50000,length = 100), TimeStopped = mean(t1_sum$TimeStopped), AID = "108", Year = 2017, id_yr = "108_2017")
pred1 <- gammit::predict_gamm(mod.DFP1$gam , newdata = newdata, se = T)
pred2 <- gammit::predict_gamm(mod.DFP2$gam , newdata = newdata, se = T)
pred3 <- gammit::predict_gamm(mod.DFP3$gam , newdata = newdata, se = T)

newdata1 <- cbind(newdata, pred1)
newdata1$Year <- 1

newdata2 <- cbind(newdata, pred2)
newdata2$Year <- 2

newdata3 <- cbind(newdata, pred3)
newdata3$Year <- 3

new <- rbind(newdata1, newdata2, newdata3)

ggplot(new) + geom_ribbon(aes(x = iyd, ymax = (prediction) + (se), ymin = (prediction) - (se), group = factor(Year)), fill = "#008B8B", alpha = .33) + geom_line(aes(x = iyd, y = (prediction), linetype = factor(Year)), size = .8) + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance (m)'~( Fidelity ^-1)), linetype = "") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + scale_linetype_manual(values = c("solid", "dashed", "dotted"), label = c("t-1", "t-2", "t-3")) + coord_cartesian(ylim = c(5,45), xlim = c(0,50000)) + scale_y_continuous(limits= c(4.8,50.2), breaks = seq(5,50,5)) + scale_x_continuous(limits= c(0,50000), breaks = seq(0,50000,10000))

ggsave(filename = "Figures/ResponseToIYD_gam_DFP_20230503.jpg", width = 24, height = 24, units = "cm", dpi = 600)

AIC(mod.absDFP)
AIC(mod.DFP)
AIC(mod.stop_km)

summary(mod.stop$lme)



E6 <- resid(mod.absDFP1$lme)
F6 <- fitted(mod.absDFP1$lme)

#validate model
par(mfrow = c(3,1), mar = c(5,5,2,2))
hist(E6,
     breaks = 25,
     main = "",
     xlab = "Residuals")
plot(x = F6,
     y = E6,
     xlab = "Fitted values",
     ylab = "Residuals")
e <- acf(E6); e
dev.off()
qqnorm(E6); qqline(E6, col="red")


E6 <- resid(mod.DFP1$lme)
F6 <- fitted(mod.DFP1$lme)

#validate model
par(mfrow = c(3,1), mar = c(5,5,2,2))
hist(E6,
     breaks = 25,
     main = "",
     xlab = "Residuals")
plot(x = F6,
     y = E6,
     xlab = "Fitted values",
     ylab = "Residuals")
e <- acf(E6); e
dev.off()
qqnorm(E6); qqline(E6, col="red")







#---------------------------#
# goodness of fit ####


mod.absDFP1.alt <- mgcv::gam(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t1_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML", family = poisson(link = "log"))
mod.absDFP2.alt <- mgcv::gam(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t2_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML", family = poisson(link = "log"))
mod.absDFP3.alt <- mgcv::gam(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped+ s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t3_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML",  family = poisson(link = "log"))

summary(mod.absDFP1.alt)$dev.expl
summary(mod.absDFP2.alt)$dev.expl
summary(mod.absDFP3.alt)$dev.expl

#without iyd
mod.absDFP1.altx <- mgcv::gam(round(absDFP,0) ~  km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t1_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML", family = poisson(link = "log"))
mod.absDFP2.altx <- mgcv::gam(round(absDFP,0) ~  km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t2_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML", family = poisson(link = "log"))
mod.absDFP3.altx <- mgcv::gam(round(absDFP,0) ~  km + TimeStopped+ s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t3_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML",  family = poisson(link = "log"))

summary(mod.absDFP1.altx)$dev.expl
summary(mod.absDFP2.altx)$dev.expl
summary(mod.absDFP3.altx)$dev.expl

mod.DFP1.alt <- mgcv::gam(DFP ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t1_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML")
mod.DFP2.alt <- mgcv::gam(DFP ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t2_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML")
mod.DFP3.alt <- mgcv::gam(DFP ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t3_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML")


summary(mod.DFP1.alt)$dev.expl
summary(mod.DFP2.alt)$dev.expl
summary(mod.DFP3.alt)$dev.expl




#-------------------------------#
# varied by greenscape ####

#chose two year since has most deviance explained
green <- RDH_14t22_greenscapes %>% rename(id_yr = AY_ID, AID = AnimalID)
t1g <- t1_sum %>% left_join(green, by = c("id_yr")) %>% ungroup() %>% mutate(iyds = scale(iyd))

#t2g <- t2g %>% mutate(GUDs = scale(greenUpDur), SPSs = scale(springScale), SVSs = scale(svSlope), GUDc = ifelse(GUDs <= 0, "short", "long"), SPSc = ifelse(SPSs <= 0,  "fast", "slow" ), SVSc = ifelse(SVSs <= 0,"rand", "ord")) 

t1g <- t1g %>% drop_na(iyd) %>% mutate(iydc = cut(iyd, breaks = quantile(iyd, seq(0, 1, by = 0.5)),include.lowest = TRUE, labels = FALSE)) %>% rename(Year = Year.x, AID = AID.x)

#
names(t2g)

# mod.absDFP.GUDc <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5, by = factor(GUDc)) + km + TimeStopped, data = t2g %>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))
# mod.absDFP.SVSc <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5, by = factor(SVSc))  + km + TimeStopped, data = t2g%>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))
# mod.absDFP.SPSc <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5, by = factor(SPSc))  + km + TimeStopped, data = t2g%>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))



mod.absDFP.GUD <- mgcv::gamm(round(absDFP,0) ~ s(greenUpDur, bs = "cs", k = 5, by = factor(iydc)), data = t1g %>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))
mod.absDFP.SPS <- mgcv::gamm(round(absDFP,0) ~ s(springScale, bs = "cs", k = 5, by = factor(iydc)), data = t1g%>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))
mod.absDFP.SVS <- mgcv::gamm(round(absDFP,0) ~ s(svSlope, bs = "cs", k = 5, by = factor(iydc)), data = t1g %>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))



newdata.gud <- data.frame(expand.grid(km = mean(t2g$km), iydc = factor(c(1,2)), TimeStopped = mean(t2g$TimeStopped), AID = "108", Year = 2017, id_yr = "108_2017", greenUpDur = seq(min(na.omit(t2g$greenUpDur)), max(na.omit(t2g$greenUpDur)), length = 50)))
newdata.svs <- data.frame(expand.grid(km = mean(t2g$km), iydc = factor(c(1,2)), TimeStopped = mean(t2g$TimeStopped), AID = "108", Year = 2017, id_yr = "108_2017", svSlope = seq(min(na.omit(t2g$svSlope)), max(na.omit(t2g$svSlope)), length = 50)))
newdata.sps = data.frame(expand.grid(km = mean(t2g$km), iydc = factor(c(1,2)), TimeStopped = mean(t2g$TimeStopped), AID = "108", Year = 2017, id_yr = "108_2017", springScale = seq(min(na.omit(t2g$springScale)), max(na.omit(t2g$springScale)), length = 50)))



predGUD <- gammit::predict_gamm(mod.absDFP.GUD$gam , newdata = newdata.gud, se = T)
predSVS <- gammit::predict_gamm(mod.absDFP.SVS$gam , newdata = newdata.svs, se = T)
predSPS <- gammit::predict_gamm(mod.absDFP.SPS$gam , newdata = newdata.sps, se = T)

newdatapredGUD <- cbind(newdata.gud, predGUD)
newdatapredSVS <- cbind(newdata.svs, predSVS)
newdatapredSPS <- cbind(newdata.sps, predSPS)




gud <- ggplot(newdatapredGUD) + geom_point(aes(x = greenUpDur, y = log(absDFP)), size = 3.5, alpha = .15, shape = 21, color = "black", fill = "grey30", data = t1g) + geom_line(aes(x = greenUpDur, y = (prediction), color = factor(iydc)), size = 1.8) + labs(y = "ln(|Days from Peak Green-up|)", title = "", x = "Green-up duration (d)", color = "Fidelity Quartile", fill = "Fidelity Quartile") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 17,color = "grey18"), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c( "#0000FF", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c( "#0000FF", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1)) + scale_x_continuous(limits = c(50,150), breaks = seq(50,150,25)); gud

svs <- ggplot(newdatapredSVS) + geom_point(aes(x = asin(svSlope^.5), y = log(absDFP)), fill = "grey30", size = 3.5, alpha = .15, shape = 21, color = "black", data = t1g) + geom_line(aes(x = asin(svSlope^.5), y = (prediction), color = factor(iydc)), size = 1.8) + labs(y = "", title = "", x = "arcsin(Green-up order)", color = "Fidelity Quartile", fill = "Fidelity Quartile") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_blank(), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c( "#0000FF", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c( "#0000FF", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1))  + scale_x_continuous(limits = c(.2,.71), breaks = seq(.2,.7,.1)); svs

sps <- ggplot(newdatapredSPS) + geom_point(aes(x = springScale, y = log(absDFP)), fill = "grey30", size = 3.5, alpha = .15, shape = 21, color = "black", data = t1g) + geom_line(aes(x = springScale, y = (prediction), color = factor(iydc)), size = 1.8) + labs(y = "", x =  bquote('Spring scale'~('green-up rate' ^-1)), color = bquote('Morrison Fidelity'~('IYD')), fill = bquote('Morrison Fidelity'~('IYD'))) + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_blank(),legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c( "#0000FF", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c( "#0000FF", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1))  + scale_x_continuous(limits = c(14.95, 22.05), breaks = seq(15,22,1)); sps




# sps <- ggplot(newdatapredSPSc) + geom_ribbon(aes(x = iyd/1000, ymax = exp(prediction) + exp(se), ymin = exp(prediction) - exp(se), group = factor(SPSc), fill = factor(SPSc)), alpha = .33) + geom_line(aes(x = iyd/1000, y = exp(prediction), fill = factor(SPSc)), size = .8) + labs(y = "", title = bquote('Green-up'~rate^-1), x = bquote('Inter-year Distance (km)'~( Fidelity ^-1)), fill = "") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"), title = element_text(size = 20,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"))+ scale_fill_manual(values = c("#00994C", "#994C00"), label = c("fast", "slow"))
# 
# svs <- ggplot(newdatapredSVSc) + geom_ribbon(aes(x = iyd/1000, ymax = exp(prediction) + exp(se), ymin = exp(prediction) - exp(se), group = factor(SVSc), fill = factor(SVSc)), alpha = .33) + geom_line(aes(x = iyd/1000, y = exp(prediction), fill = factor(SVSc)), size = .8) + labs(y = "", x = "", fill = "", title = "Green-up Order") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c("#994C00", "#00994C"), label = c("random", "ordered"))


gud + svs + sps

ggsave(filename = "Figures/GreenscapeMetricsByAvg_20230506.jpg", width = 72, height = 24, units = "cm", dpi = 600)



#--------------------------------#
# what does high Morrison mean ####

stoploc <- hifid1 %>% mutate(stop.n.c = stop.n_1, km_mark = km_mark_1, TimeStopped = TimeStopped_1, Year = curr_Y) %>% group_by(stop.n.c) %>% mutate(km = round(mean(km_mark),0)) %>% select(id_yr, stop.n.c, km, TimeStopped, Year, prev_Y_1) %>% group_by(stop.n.c) %>% slice(1) %>% st_drop_geometry()

head(stoploc)


topbott <- t1_sum %>% drop_na(iyd) %>% mutate(km = round(km,0), Year = as.numeric(Year), prev_Y_1 = Year - 1) 

#create a dataframe where we have iyd, this years stop, last years stop, then we can find the minimum distance between a given stop and all stops in the previous year thus getting at distance between stops

topbott_now <- topbott %>% left_join(stoploc %>% dplyr::select(-prev_Y_1) %>% mutate(Year = as.numeric(Year)), by = c("stop.n.c", "km", "Year", "id_yr"))
topbott_prev <- topbott %>% left_join(stoploc %>% dplyr::select(-Year) %>% mutate(prev_Y_1 = as.numeric(prev_Y_1)), by = c("stop.n.c", "km", "prev_Y_1", "id_yr"))


topbott_prev



#-------------------------#
# landscape test ####




load("./REnvs/ImageOfClimateSampling_20230428.RDS")



clim1 <- t1_sum %>% mutate(km = round(km, 0 ) ) %>% as.data.frame() %>% left_join(clim.var.data %>% group_by(km) %>% summarise(sd = mean(sd)), by = "km",relationship = "many-to-many")
clim2 <- t2_sum %>% mutate(km = round(km, 0 ) ) %>% as.data.frame() %>% left_join(clim.var.data %>% group_by(km) %>% summarise(sd = mean(sd)), by = "km",relationship = "many-to-many")
clim3 <- t3_sum %>% mutate(km = round(km, 0 ) ) %>% as.data.frame() %>% left_join(clim.var.data %>% group_by(km) %>% summarise(sd = mean(sd)), by = "km",relationship = "many-to-many")


#plot
clim1 <- clim1 %>% mutate(fiyd = as.factor(ifelse(iyd <= 5000, 1, 0)))
clim2 <- clim2 %>% mutate(fiyd = as.factor(ifelse(iyd <= 5000, 1, 0)))
clim3 <- clim3 %>% mutate(fiyd = as.factor(ifelse(iyd <= 5000, 1, 0)))

ggplot(clim1 %>% filter(!is.na(fiyd))) + geom_point(aes(x = sd, y = TimeStopped))

ggplot(clim1)  + geom_density(aes(x = sd, fill = fiyd), alpha = .5) + geom_vline(aes(xintercept = meSD, color = fiyd), size = 1,  data = clim1 %>% group_by(fiyd) %>% summarise(meSD = median(na.omit(sd))))
ggplot(clim2) + geom_density(aes(x = sd, fill = fiyd), alpha = .5) + geom_vline(aes(xintercept = meSD, color = fiyd), size = 1,  data = clim2 %>% group_by(fiyd) %>% summarise(meSD = median(na.omit(sd))))
ggplot(clim3) + geom_density(aes(x = sd, fill = fiyd), alpha = .5) + geom_vline(aes(xintercept = meSD, color = fiyd), size = 1,  data = clim3 %>% group_by(fiyd) %>% summarise(meSD = median(na.omit(sd))))


set.seed(2820)
sd1 <- clim1 %>% drop_na(sd) %>% group_by(fiyd) %>% sample_n(150)
sd2 <- clim2 %>% drop_na(sd) %>% group_by(fiyd) %>% sample_n(100)
sd3 <- clim3 %>% drop_na(sd) %>% group_by(fiyd) %>% sample_n(40)

t.test(sd1$sd ~ sd1$fiyd, paired = TRUE)
t.test(sd2$sd ~ sd2$fiyd, paired = TRUE)
t.test(sd3$sd ~ sd3$fiyd, paired = TRUE)



# hypothesis two, dealing with different conditions
t

green <- RDH_14t22_greenscapes %>% rename(id_yr = AY_ID, AID = AnimalID)
t2g <- t2 %>% left_join(green, by = c("id_yr")) %>% ungroup() %>% mutate(ciyds = scale(ciyd))

t2g <- t2g %>% mutate(GUDs = scale(greenUpDur), SPSs = scale(springScale), SVSs = scale(svSlope), GUDc = ifelse(GUDs <= 0, "short", "long"), SPSc = ifelse(SPSs <= 0,  "fast", "slow" ), SVSc = ifelse(SVSs <= 0,"rand", "ord")) %>% rename(Year = Year.x, AID = AID.x)


mod.IYDbyGUD <- mgcv::gamm(absDFP ~ s(ciyds, bs = "cs", k = 5, by = factor(GUDc)) +  SVSs, method = "REML", random = list(id_yr=~ciyds|id_yr, Year = ~ciyds|Year, AID=~1), family = poisson(link = "log"), data = t2g, niterPQL = 20) #SPSs +

mod.IYDbySPS <- mgcv::gamm(absDFP ~ s(ciyds, bs = "cs", k = 5, by = factor(SPSc)) + GUDs , method = "REML", random = list(id_yr=~ciyds|id_yr, Year = ~ciyds|Year, AID=~1), family = poisson(link = "log"), data = t2g, niterPQL = 20)

mod.IYDbySVS <- mgcv::gamm(absDFP ~ s(ciyds, bs = "cs", k = 5, by = factor(SVSc)) + GUDs, method = "REML", random = list(id_yr=~ciyds|id_yr, Year = ~ciyds|Year, AID=~1), family = poisson(link = "log"), data = t2g, niterPQL = 20) #SPSs +

summary(mod.IYDbyGUD$gam)
summary(mod.IYDbySPS$gam)
summary(mod.IYDbySVS$gam)

summary(mod.IYDbyGUD$lme)
summary(mod.IYDbySPS$lme) # top model
summary(mod.IYDbySVS$lme)

AICc(mod.IYDbyGUD$lme)
AICc(mod.IYDbySPS$lme)
AICc(mod.IYDbySVS$lme)

par(mfrow = c(3,2))
plot.gam(mod.IYDbyGUD$gam)
plot.gam(mod.IYDbySPS$gam)
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



