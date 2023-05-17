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

load("REnvs/EnvForRunningMostGams_20230509.RDS")

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
#remotes::install_github("stefanocoretta/tidygam@v0.1.0")
library(tidygam)
library(StatisticalModels)
library(patchwork)
library(gammit)

#--------------------------#
# prep datasets ####

setwd("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/")
load("Data_out/Data/Stopover/RDH_StopoverBundle_14t22_20230511.RData")
load("Data_out/Data/Stopover/FidelityMetrics_HighUse_20230513.RData")
#load("Data_out/Data/Stopover/RDH_AlternateStops_Bischof_20230428.RData")
load("Data_out/Data/Greenscape/RDH_14t22_greenscapes_20230425.RData")
#load("Data_out/Data/Stopover/RDH_UsedAvailStops_14t22_20230428.RData")
load("C:/Users/lwilde2/Desktop/RDH Database/Processed Data/RDH_AllMigrations_Bischof_2014t2022_20230425.RData")

load("Data_out/Data/Stopover/FidelityDataset_HighUse_20230513.RData")


#-----------------------#
# remove indiv with incomp migs ####

sumcheck <- data2 %>% group_by(id_yr) %>% summarise(start = min(km_mark), end = max(km_mark))
remove <- sumcheck %>% filter(start > 100 | end < 150) %>% dplyr::select(id_yr) %>% st_drop_geometry()

remove <- remove %>% mutate(AID = str_sub(id_yr, 0,3))

IYD_stop_fidelity <- IYD_stop_fidelity %>% filter(!AID.1 %in% remove$AID)


#------------------------------#
#histograms of stops ####

num.migs <- onhigh_yr_fin %>% group_by(Year) %>% summarise(count = length(unique(id_yr)))
num.km <- onhigh_yr_fin %>% mutate(km = round(km_mark, 0)) %>% group_by(Year, km) %>% summarise(count = length(unique(id_yr)))

data.stops.corridor <- num.km %>% st_drop_geometry() %>% left_join(num.migs %>% st_drop_geometry(), by = "Year") %>% mutate(prop = count.x / count.y)


ggplot(data.stops.corridor, aes(x = km, y = Year, group = Year)) + geom_density_ridges2(aes(fill = Year)) + scale_x_continuous(limits = c(0,240), expand = c(0,0)) #+ stat_smooth(geom = 'area', method = 'loess', alpha = 1/2, aes(fill = Year)) 


#---------------------------------#
# Calcs per stop and id_yr ####

IYD_stop_fidelity$meanDev <- rowMeans(IYD_stop_fidelity[,c("dev_1","dev_2","dev_3","dev_4","dev_5","dev_6","dev_7")], na.rm = TRUE)
IYD_stop_fidelity$meanIYD <- rowMeans(IYD_stop_fidelity[,c("iyd_1","iyd_2","iyd_3","iyd_4","iyd_5","iyd_6","iyd_7")], na.rm = TRUE)


summary(IYD_stop_fidelity$meanDev)
summary(IYD_stop_fidelity$meanIYD)

sumid_yr <- IYD_stop_fidelity %>% mutate(id_yr = paste(AID.1,"_",curr_Y.1,sep = "")) %>% group_by(id_yr) %>% summarise(OverMeanDev = mean(meanDev), OverMeanIYD = mean(meanIYD), OverSumDev = sum(meanDev))
sumaid <- IYD_stop_fidelity %>% mutate(id_yr = paste(AID.1,"_",curr_Y.1,sep = "")) %>% group_by(AID.1) %>% summarise(OverMeanDev = mean(meanDev,na.rm = TRUE), OverMeanIYD = mean(meanIYD,na.rm = TRUE), OverSumDev = sum(meanDev,na.rm = TRUE)); sumaid


fid1 <- IYD_stop_fidelity_list[[1]]
fid2 <- IYD_stop_fidelity_list[[2]]
fid3 <- IYD_stop_fidelity_list[[3]]
fid4 <- IYD_stop_fidelity_list[[4]]
fid5 <- IYD_stop_fidelity_list[[5]]
fid6 <- IYD_stop_fidelity_list[[6]]
fid7 <- IYD_stop_fidelity_list[[7]]

#time window
fid1 <- fid1 %>% filter((startJul_1 - startJul_prev_1) < 30)
fid2 <- fid2 %>% filter((startJul_2 - startJul_prev_2) < 30)
fid3 <- fid3 %>% filter((startJul_3 - startJul_prev_3) < 30)
fid4 <- fid4 %>% filter((startJul_4 - startJul_prev_4) < 30)
fid5 <- fid5 %>% filter((startJul_5 - startJul_prev_5) < 30)
fid6 <- fid6 %>% filter((startJul_6 - startJul_prev_6) < 30)
fid7 <- fid7 %>% filter((startJul_7 - startJul_prev_7) < 30)


#greenscape data
RDH_14t22_greenscapes
#stops
onstop_yr_fin
#full data
data1
#fidelity
fid1
fid2
fid3




fid1 <- fid1 %>% arrange(id_yr, km_mark_1) %>% group_by(id_yr) %>% mutate(StopNumCurr = row_number()) %>% ungroup()
fid2 <- fid2 %>% arrange(id_yr, km_mark_2) %>% group_by(id_yr) %>% mutate(StopNumCurr = row_number()) %>% ungroup()
fid3 <- fid3 %>% arrange(id_yr, km_mark_3) %>% group_by(id_yr) %>% mutate(StopNumCurr = row_number()) %>% ungroup()



#remove problem deer with incomplete migrations

fid1 <- fid1 %>% filter(!AID.1 %in% remove$AID)
fid2 <- fid2 %>% filter(!AID.2 %in% remove$AID)
fid3 <- fid3 %>% filter(!AID.3 %in% remove$AID)
fid4 <- fid4 %>% filter(!AID.4 %in% remove$AID)
fid5 <- fid5 %>% filter(!AID.5 %in% remove$AID)
fid6 <- fid6 %>% filter(!AID.6 %in% remove$AID)

length(unique(fid1$AID.1)); length(unique(fid1$id_yr))
length(unique(fid2$AID.2)); length(unique(fid2$id_yr))
length(unique(fid3$AID.3)); length(unique(fid3$id_yr))
length(unique(fid4$AID.4)); length(unique(fid4$id_yr))
length(unique(fid5$AID.5)); length(unique(fid5$id_yr))
length(unique(fid6$AID.6)); length(unique(fid6$id_yr))
length(unique(fid7$AID.7)); length(unique(fid7$id_yr))

#adaptability of fidelity?
#absDFP ~ greenscape * IYD
#DFP ~ greenscape * IYD

#origins of fidelity?
#landscape
#response to mismatch

length(Reduce(unique, list(unique(fid1$AID.1), unique(fid2$AID.2), unique(fid3$AID.3), unique(fid4$AID.4), unique(fid5$AID.5))))
length(Reduce(unique, list(unique(fid1$id_yr), unique(fid2$id_yr), unique(fid3$id_yr), unique(fid4$id_yr), unique(fid5$id_yr))))

#-------------------------------#
# consistency in fidelity ####

#to km

#in individuals

#process data
# 
# all <- do.call(rbind, iyddat)
# 
# all$ydiff <- as.numeric(all$curr_Y) - as.numeric(all$prev_Y)
# 
# all.sum <- all %>% group_by(id_yr, ydiff) %>% summarise(stop.n = n(), mean.iyd = mean(iyd))
# 
# ggplot(all.sum) + geom_line(aes(x = ydiff, y = mean.iyd/1000, color = id_yr), linewidth = 1.4, alpha = .5) + theme_classic() + theme(legend.position = "none") + labs(x = "Years Between Migration Events", y = bquote('Mean Inter-year Distance (km) '~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + coord_cartesian(xlim = c(1,7), ylim = c(0,80)) + scale_x_continuous(expand = c(0,0), limits = c(1,7), breaks = seq(1,7,1))  + scale_y_continuous(expand = c(0,0), limits = c(0,80), breaks = seq(0,80,10))
# 
# head(all)
# 
# mod.km.eff <- mgcv::gamm(round(iyd,0) ~ s(km_mark, by = factor(ydiff), bs = "cs", k = 5), method = "GCV", random = list(id_yr=~1, curr_Y = ~1, AID=~1), family = poisson(link = "log"), data = all %>% filter(km_mark < 240))
# mod.km.eff2 <- mgcv::gamm(round(iyd,0) ~ s(km_mark, bs = "cs", k = 5), method = "GCV", random = list(id_yr=~1, curr_Y = ~1, AID=~1), weights = 1/ydiff, family = poisson(link = "log"), data = all %>% filter(km_mark < 240))
# 
# plot.gam(mod.km.eff2$gam)
# 
# newdata = data.frame(expand.grid(km_mark = seq(1,240,by = 1), ydiff = seq(1,7,by = 1), AID = "108", curr_Y = 2017, id_yr = "108_2017"))
# newdata$pred <- gammit::predict_gamm(mod.km.eff2$gam , newdata = newdata, se = T)$prediction
# newdata$se <- gammit::predict_gamm(mod.km.eff2$gam , newdata = newdata, se = T)$se
# 
# ggplot(data = newdata) + geom_line(aes(x = km_mark, y = pred), linewidth = 1.2) + geom_ribbon(aes(x = km_mark, ymin = pred - se, ymax = pred + se), alpha = .2) + theme_classic()
# 
# 



# ggplot(fid1) + theme(legend.position = "none") + coord_cartesian(xlim = c(0,220), ylim = c(0,100000)) + geom_point(aes(x = km_mark_1, y = iyd_1, group = AID, color = AID),  alpha = .2) + geom_smooth(aes(x = km_mark_1, y = iyd_1)) #geom_line(aes(x = km_mark_1, y = iyd_1, group = AID, color = AID)) +


# First hypothesis, accumulated dev should hurt ability to surf - decay of information?

sumfid1 <- fid1 %>% group_by(id_yr) %>% summarise(Duration = (max(endJul_1) - min(startJul_1)), SumDev = sum(dev_1)/Duration)
fid1 <- fid1 %>% group_by(id_yr) %>% mutate(cumulDev = cumsum(dev_1))

sumfid2 <- fid2 %>% group_by(id_yr) %>% summarise(Duration = (max(endJul_2) - min(startJul_2)), SumDev = sum(dev_2)/Duration)
fid2 <- fid2 %>% group_by(id_yr) %>% mutate(cumulDev = cumsum(dev_2))

sumfid3 <- fid3 %>% group_by(id_yr) %>% summarise(Duration = (max(endJul_3) - min(startJul_3)), SumDev = sum(dev_3)/Duration)
fid3 <- fid3 %>% group_by(id_yr) %>% mutate(cumulDev = cumsum(dev_3))


t1 <- onhigh_yr_fin %>% left_join(fid1 %>% dplyr::select(iyd_1, dev_1, cumulDev, stop.n_1, km_mark_1) %>% rename(stop.n.c = stop.n_1), by = c("stop.n.c", "id_yr")) %>% arrange(id_yr, POSIXct)
t1 <- t1 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

t2 <- onhigh_yr_fin %>% left_join(fid2 %>% dplyr::select(iyd_2, dev_2, cumulDev, stop.n_2, km_mark_2) %>% rename(stop.n.c = stop.n_2), by = c("stop.n.c", "id_yr")) %>% arrange(id_yr, POSIXct)
t2 <- t2 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP))

t3 <- onhigh_yr_fin %>% left_join(fid3 %>% dplyr::select(iyd_3, dev_3, cumulDev, stop.n_3, km_mark_3) %>% rename(stop.n.c = stop.n_3), by = c("stop.n.c", "id_yr")) %>% arrange(id_yr, POSIXct)
t3 <- t3 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP)) 



#need to model this, with random effects

#averaging to each stop
t1_sum <- t1 %>% group_by(stop.n.c) %>% summarise(absDFP = round(mean(absDFP),0), DFP = round(mean(DFP),0), TimeStopped = unique(TimeStopped), id_yr = unique(id_yr), AID = unique(AID), Year = unique(Year), iyd = unique(iyd_1), dev = unique(dev_1), cumdev = unique(cumulDev), km = mean(km_mark))


t2_sum <- t2 %>% group_by(stop.n.c) %>% summarise(absDFP = round(mean(absDFP),0), DFP = round(mean(DFP),0), TimeStopped = unique(TimeStopped), id_yr = unique(id_yr), AID = unique(AID), Year = unique(Year), iyd = unique(iyd_2), dev = unique(dev_2), cumdev = unique(cumulDev), km = mean(km_mark))


t3_sum <- t3 %>% group_by(stop.n.c) %>% summarise(absDFP = round(mean(absDFP),0), DFP = round(mean(DFP),0), TimeStopped = unique(TimeStopped), id_yr = unique(id_yr), AID = unique(AID), Year = unique(Year), iyd = unique(iyd_3), dev = unique(dev_3), cumdev = unique(cumulDev), km = mean(km_mark))


compensation <- data2 %>% mutate(DFP = JDate - MaxIRGday, absDFP = abs(DFP)) %>% arrange(id_yr, POSIXct) %>% group_by(id_yr) %>% summarise(fDFP = first(DFP), lDFP = last(DFP), start =first(JDate), end = last(JDate), Comp = (last(DFP) - first(DFP)), CompAbs = (last(absDFP) - first(absDFP))) %>% st_drop_geometry()



t1_sum <- t1_sum %>% left_join(compensation, by = "id_yr")
t2_sum <- t2_sum %>% left_join(compensation, by = "id_yr")
t3_sum <- t3_sum %>% left_join(compensation, by = "id_yr")

t1_sum_timing <- t1_sum %>% filter(Year != 2015) %>% group_by(Year) %>% mutate(timing.c = cut(start, breaks = quantile(start, c(0,.25,.75,1.00)),include.lowest = TRUE, labels = FALSE))
t2_sum_timing <- t2_sum %>% filter(Year != 2015) %>% group_by(Year) %>% mutate(timing.c = cut(start, breaks = quantile(start, c(0,.25,.75,1.00)),include.lowest = TRUE, labels = FALSE))
t3_sum_timing <- t3_sum %>% filter(!Year %in% c(2015)) %>% group_by(Year) %>% mutate(timing.c = cut(start, breaks = quantile(start, c(0,.25,.75,1.00)),include.lowest = TRUE, labels = FALSE))


#-----------------------------------#
# ability to compensate ####

t1_comp <- t1_sum_timing %>% st_drop_geometry() 
t2_comp <- t2_sum_timing %>% ungroup() %>% st_drop_geometry() %>% rename(iyd2 = iyd, dev2 = dev) %>% select(iyd2, dev2, stop.n.c, TimeStopped)
t3_comp <- t3_sum_timing %>% ungroup() %>% st_drop_geometry() %>% rename(iyd3 = iyd, dev3 = dev) %>% select(iyd3, dev3, stop.n.c, TimeStopped)

t_all_comp <- t1_comp %>% left_join(t2_comp, by = "stop.n.c") %>% left_join(t3_comp, by = "stop.n.c")

t_all_comp$meanIYD <- rowMeans(t_all_comp[c("iyd", "iyd2", "iyd3")], na.rm = T)
t_all_comp$meanDev <- rowMeans(t_all_comp[c("dev", "dev2", "dev3")], na.rm = T)

#t_all_comp <- t_all_comp %>% drop_na(meanIYD, meanDev)

t_all_comp_sum <- t_all_comp %>% ungroup() %>% group_by(id_yr) %>% dplyr::summarise(AID = unique(AID), Year = unique(Year), sumIYD = sum(meanIYD, na.rm = T), meanIYD = mean(meanIYD), dev = mean(dev), meanDev = mean(meanDev), cumdev = max(cumdev)/(sum(TimeStopped)/24), Comp= unique(Comp), CompAbs = unique(CompAbs), timing.c = unique(timing.c)) %>% mutate(fComp = ifelse(timing.c > 1, ifelse(Comp < 0,1,0), ifelse(Comp < 0,0,1))); t_all_comp_sum



compdata <- read.csv("C:/Users/lwilde2/Documents/PhD_AdaptiveFidelity/Data_out/Data/Compensation/AOrtega_RDHDeerCompensationStatus.csv")

t_all_comp_sumx <- t_all_comp_sum %>% left_join(compdata, by = "id_yr") %>% drop_na(CompFact)

table(t_all_comp_sumx$CompFact)

modbinom <- mgcv::gamm(CompFact ~ s(meanDev, by = factor(timing.c), bs = "cs", k = 10), list(Year = ~1, AID = ~1), family = binomial(link = "cloglog"), data = t_all_comp_sumx)
summary(modbinom$gam)


# moddegree <- mgcv::gamm(Comp ~ s(meanDev, by = factor(timing.c), bs = "cs", k = 5), list(Year = ~1, AID = ~1), data = t_all_comp %>% filter(km > 50 & km < 210))
# summary(moddegree$gam)

ggplot(t_all_comp) + geom_point(aes(x = km, y = dev, size = log(absDFP), color = factor(timing.c))) + geom_line(aes(x = km, y = dev, group = factor(id_yr), color = factor(timing.c))) + coord_cartesian(xlim = c(50,200))



newdata = data.frame(expand.grid(meanDev = seq(min(t_all_comp_sum$meanDev, na.rm = T), max(t_all_comp_sum$meanDev, na.rm = T), length = 100), timing.c = c(1,2,3), AID = 108, Year = 2019))
predcumdev <- gammit::predict_gamm(modbinom$gam , newdata = newdata, se = T)

clog.inv <- function(values){
  e1<- (-1)*(exp(values))
  1-exp(e1)
}

predcumdev$prediction <- clog.inv(predcumdev$prediction)
predcumdev$se <- clog.inv(predcumdev$se)/3

ggplot(cbind(newdata, predcumdev)) + geom_line(aes(x = meanDev, y = prediction, color = factor(timing.c), group = factor(timing.c)), size = 1.2) + geom_ribbon(aes(x = meanDev, ymin = (prediction - se), ymax = (prediction + se), group = factor(timing.c)), color = "grey60", alpha = .2) + coord_cartesian(xlim = c(0,4000), ylim = c(0,1)) + scale_x_continuous(expand = c(0,0), limits = c(0,4000), breaks = seq(0,8000,500)) + scale_y_continuous(expand = c(0,0), limits = c(-.5,1.5), breaks = seq(0,1,.2)) + theme_classic()



modIYDmismatch <- mgcv::gamm(absDFP ~ s(meanIYD, by = factor(timing.c), bs = "cs", k = 5), list(id_yr = ~1, Year = ~1, AID = ~1), data=t_all_comp, family = poisson(link = "log"))
summary(modIYDmismatch$gam)

newdata = data.frame(expand.grid(meanIYD = seq(min(t_all_comp$meanIYD, na.rm = T), max(t_all_comp$meanIYD, na.rm = T), length = 100), timing.c = c(1,2,3), AID = 108, Year = 2019, id_yr = "108_2019"))
predIYDmis <- gammit::predict_gamm(modIYDmismatch$gam , newdata = newdata, se = T)

ggplot(cbind(newdata, predIYDmis)) + geom_line(aes(x = meanIYD, y = exp(prediction), color = factor(timing.c), group = factor(timing.c)), size = 1.2) + geom_ribbon(aes(x = meanIYD, ymin = exp(prediction - se), ymax = exp(prediction + se), group = factor(timing.c)), color = "grey60", alpha = .2) + coord_cartesian(xlim = c(0,4000), ylim = c(-25,150)) + scale_x_continuous(expand = c(0,0), limits = c(0,4100), breaks = seq(0,4000,500)) + scale_y_continuous(expand = c(0,0), limits = c(-40,160), breaks = seq(-25,150,25)) + theme_classic()


# 
# modbinomIYD <- mgcv::gamm(Comp ~ s(sumDev, by = factor(timing.c), bs = "cs", k = 5), list(Year = ~1, AID = ~1), data = t_all_comp_sum)
# summary(modbinomIYD$gam)
# 
# newdata = data.frame(expand.grid(sumDev = seq(min(t_all_comp_sum$sumDev, na.rm = T), max(t_all_comp_sum$sumDev, na.rm = T), length = 100), timing.c = c(1,2,3), AID = 108, Year = 2019))
# predcumdev <- gammit::predict_gamm(modbinomIYD$gam , newdata = newdata, se = T)
# 
# clog.inv <- function(values){
#   e1<- (-1)*(exp(values))
#   1-exp(e1)
# }
# 
# predcumdev$prediction <- predcumdev$prediction
# predcumdev$se <- predcumdev$se

ggplot(cbind(newdata, predcumdev)) + geom_line(aes(x = meanDev, y = prediction, color = factor(timing.c), group = factor(timing.c)), size = 1.2) + geom_ribbon(aes(x = meanDev, ymin = (prediction - se), ymax = (prediction + se), group = factor(timing.c), fill = factor(timing.c)), alpha = .2) + coord_cartesian(xlim = c(0,3000), ylim = c(0,1.05)) + scale_x_continuous(expand = c(0,0), limits = c(0,9100), breaks = seq(0,3000,500)) + scale_y_continuous(expand = c(0,0), limits = c(-.5,15), breaks = seq(0,1,.2)) + theme_classic() + theme(axis.title.x = element_text(size = 40,color = "grey18"), axis.title.y = element_text(size = 40,color = "grey18"),axis.text.x = element_text(size = 30,color = "grey18"),axis.text.y = element_text(size = 30,color = "grey18"), legend.text = element_text(size = 30,color = "grey18"), legend.title = element_text(size = 40,color = "grey18"),legend.key.size = unit(1, 'cm'), legend.key.height = unit(1, 'cm'), legend.key.width = unit(2, 'cm')) + labs(x = "Mean Route Deviation (m)", y = "Probability of Compensation", fill = "Migration Start Date", color = "Migration Start Date") + scale_color_manual(values = c("#4DBBD5FF","#3C5488FF","black"), labels = c("early", "mid", "late")) + scale_fill_manual(values = c("#4DBBD5FF","#3C5488FF","black"), labels = c("early", "mid", "late")) 


ggsave(filename = "Figures/ProbabilityOfCompensation_20230517.jpg", width = 40, height = 24, units = "cm", dpi = 600)

comptest

mod1 <- mgcv::gamm(DFP ~ s(dev, by = factor(timing.c), bs = "cs", k = 5), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), data = t1_sum_timing)


k1 <- t1_sum_timing %>% group_by(id_yr) %>% summarise(Comp = unique(Comp), cumdev = max(cumdev), timing.c = unique(timing.c), Year = factor(unique(Year)), AID = factor(unique(AID)), fComp = ifelse(Comp >= 0, 0, 1))

ggplot(k1) + geom_point(aes(x = cumdev, y = fComp))

modbinom <- mgcv::gamm(fComp ~ s(cumdev, by = factor(timing.c), bs = "cs", k = 5), list(Year = ~1, AID = ~1), family = binomial(link = "logit"), data = k1)
summary(modbinom$gam)

newdata = data.frame(expand.grid(cumdev = seq(min(k1$cumdev, na.rm = T), max(k1$cumdev, na.rm = T), length = 100), timing.c = c(1,2,3), AID = 108, Year = 2019))
predcumdev <- gammit::predict_gamm(modbinom$gam , newdata = newdata, se = T)

clog.inv <- function(values){
  e1<- (-1)*(exp(values))
  1-exp(e1)
}

predcumdev$prediction <- clog.inv(predcumdev$prediction)
predcumdev$se <- clog.inv(predcumdev$se)

ggplot(cbind(newdata, predcumdev)) + geom_line(aes(x = cumdev, y = prediction, color = factor(timing.c), group = factor(timing.c)), size = 1.2) + geom_ribbon(aes(x = cumdev, ymin = (prediction - se), ymax = (prediction + se), group = factor(timing.c)), color = "grey60", alpha = .2) + scale_x_continuous(expand = c(0,0), limits = c(0,8000), breaks = seq(0,8000,2000)) + scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1,.2)) + theme_classic()



# mod1 <- mgcv::gamm(Comp ~ s(cumdev, by = factor(timing.c), bs = "cs", k = 5), method = "REML", random = list(Year = ~1, AID=~1), data = k1)
# plot(mod1$gam)
# summary(mod1$gam)
# 
# newdata = data.frame(expand.grid(cumdev = seq(min(k1$cumdev, na.rm = T), max(k1$cumdev, na.rm = T), length = 100), timing.c = c(1,2,3), AID = "108", Year = 2017))
# predcumdev <- gammit::predict_gamm(mod1$gam , newdata = newdata, se = T)
# 
# ggplot(cbind(newdata, predcumdev)) + geom_line(aes(x = cumdev, y = prediction, color = factor(timing.c), group = factor(timing.c)), size = 1.2) + geom_ribbon(aes(x = cumdev, ymin = (prediction - se), ymax = (prediction + se), group = factor(timing.c)), color = "grey60", alpha = .2) + theme_classic() + labs(x = "Cumulative Deviations (corrected)", y = "Degree of Compensation (Days from Peak Green-up)") + geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1.2) + scale_x_continuous(expand = c(0,0), limits = c(0,8000), breaks = seq(0,8000,2000)) + scale_y_continuous(expand = c(0,0), limits = c(5,100), breaks = c(0,20,40,60,80,100))


hist(t_all_comp$meanIYD)

mod.adap.1 <- mgcv::gamm(absDFP ~ s(dev, bs = "cs", k = 5), random = list(id_yr = ~1, Year = ~1, AID=~1), family = poisson(link = "log"), data = t1_sum_timing) #by = factor(timing.c), 
mod.adap.2 <- mgcv::gamm(absDFP ~ s(meanIYD, bs = "cs", by = factor(timing.c), k = 5), random = list(id_yr = ~1, Year = ~1, AID=~1), family = poisson(link = "log"), data = t_all_comp %>% filter(meanIYD< 40000))

summary(mod.adap.1$gam)

newdata1 = data.frame(expand.grid(dev = seq(min(t1_sum_timing$dev, na.rm = T), max(t1_sum_timing$dev, na.rm = T), length = 100), timing.c = c(1,2,3), AID = "108", Year = 2017))
preddev <- gammit::predict_gamm(mod.adap.1$gam , newdata = newdata1, se = T)

ggplot(cbind(newdata1, preddev)) + geom_line(aes(x = dev/1000, y = exp(prediction), group = factor(timing.c)), size = 1.5) + geom_ribbon(aes(x = dev/1000, ymin = exp(prediction - se), ymax = exp(prediction + se), group = factor(timing.c)), fill = "#008B8B", alpha = 0.23) + labs(y = "|Days from Peak Green-up|", x = "Deviation From Remembered Route (km)", linetype = "") + theme_classic()  + theme(axis.title.x = element_text(size = 24,color = "grey18"), axis.title.y = element_text(size = 24,color = "grey18"),axis.text.x = element_text(size = 16,color = "grey18"),axis.text.y = element_text(size = 16,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + coord_cartesian(ylim = c(20,50), xlim = c(0,8)) + scale_y_continuous(expand = c(0,0), limits= c(19,51), breaks = seq(20,50,5)) + scale_x_continuous(expand = c(0,0), limits= c(0,8), breaks = seq(0,8,1)) # + scale_linetype_manual(values = c("solid", "dashed", "dotted"), label = c("Early", "Middle", "Late")) + coord_cartesian(ylim = c(20,140), xlim = c(0,8)) + scale_y_continuous(expand = c(0,0), limits= c(10,140), breaks = seq(20,140,20)) + scale_x_continuous(expand = c(0,0), limits= c(0,8), breaks = seq(0,8,1))

ggsave(filename = "Figures/CatchingWave_Deviation_20230513.jpg", width = 24, height = 24, units = "cm", dpi = 600)

newdata2 = data.frame(expand.grid(meanIYD = seq(0, 40000, length = 100), timing.c = c(1,2,3), AID = "108", Year = 2017))
preddev2 <- gammit::predict_gamm(mod.adap.2$gam , newdata = newdata2, se = T)

ggplot(cbind(newdata2, preddev2)) + geom_line(aes(x = meanIYD/1000, y = exp(prediction), linetype = factor(timing.c), group = factor(timing.c)), size = 1.5) + geom_ribbon(aes(x = meanIYD/1000, ymin = exp(prediction - se), ymax = exp(prediction + se), group = factor(timing.c)), fill = "#008B8B", alpha = .53) + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance (km)'~( Fidelity ^-1)), linetype = "") + theme_classic()  + theme(axis.title.x = element_text(size = 40,color = "grey18"), axis.title.y = element_text(size = 40,color = "grey18"),axis.text.x = element_text(size = 30,color = "grey18"),axis.text.y = element_text(size = 30,color = "grey18"),legend.text = element_text(size = 30,color = "grey18"), legend.title = element_text(size = 40,color = "grey18"),legend.key.size = unit(1, 'cm'), legend.key.height = unit(1, 'cm'), legend.key.width = unit(2, 'cm')) + scale_linetype_manual(values = c("solid", "dashed", "dotted"), label = c("Early", "Middle", "Late")) + coord_cartesian(ylim = c(5,61), xlim = c(0,120)) + scale_y_continuous(expand = c(0,0), limits= c(5,60.5), breaks = seq(5,60,5)) + scale_x_continuous(limits= c(0,120), breaks = seq(0,120,30))

ggsave(filename = "Figures/CatchingWave_IYD_20230514.jpg", width = 24, height = 24, units = "cm", dpi = 600)

# 
# #lmer of relationship
# #average per stop, and what about including timestopped?
# 
# m1 <- glmer(absDFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), family = poisson(link = "log"), data = t1_sum ) #%>% filter(iyd < 20000)
# summary(m1)
# 
# m2 <- glmer(absDFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), family = poisson(link = "log"), data = t2_sum) #%>% filter(iyd < 20000)
# summary(m2)
# 
# m3 <- glmer(absDFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), family = poisson(link = "log"), data = t3_sum) #%>% filter(iyd < 20000)
# summary(m3)
# 
# plot(DHARMa::simulateResiduals(m1))
# plot(DHARMa::simulateResiduals(m2))
# plot(DHARMa::simulateResiduals(m3))
# 
# #rely on trigamma for poisson, needs a log link function
# r.squaredGLMM(m1)
# r.squaredGLMM(m2)
# r.squaredGLMM(m3)
# 
# 
# newdata <- data.frame(absDFP = seq(min(na.omit(t1_sum$absDFP)), max(na.omit(t1_sum$absDFP)), length = 50), iyd = seq(min(na.omit(t1_sum$iyd)), max(na.omit(t1_sum$iyd)), length = 50), id_yr = "108_2018", AID = "108", Year = 2016)
# 
#  p1 <- PredictGLMER(m1, newdata, se.fit = T, seMultiplier = 1.96)
#  p2 <- PredictGLMER(m2, newdata, se.fit = T, seMultiplier = 1.96)
#  p3 <- PredictGLMER(m3, newdata, se.fit = T, seMultiplier = 1.96)
# 
# 
# 
# 
# plot1 <- ggplot() + geom_ribbon(data = cbind(p1, newdata), aes(x = iyd, ymin = exp(yminus), ymax = exp(yplus)), fill = "#008B8B", alpha = .75) + geom_line(data = cbind(p1, newdata), aes(x = iyd, y=exp(y)), color = "black") + theme_classic() + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + coord_cartesian(xlim = c(0,50000), ylim = c(20,40)) 
# 
# plot2 <- ggplot() + geom_ribbon(data = cbind(p2, newdata), aes(x = iyd, ymin = exp(yminus), ymax = exp(yplus)), fill = "#008B8B", alpha = .3) + geom_line(data = cbind(p2, newdata), aes(x = iyd, y=exp(y)), color = "black") + theme_classic() + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"))+ coord_cartesian(xlim = c(0,50000), ylim = c(20,40))
# 
# plot3<- ggplot() + geom_ribbon(data = cbind(p3, newdata), aes(x = iyd, ymin = exp(yminus), ymax = exp(yplus)), fill = "#008B8B", alpha = .1) + geom_line(data = cbind(p3, newdata), aes(x = iyd, y=exp(y)), color = "black") + theme_classic() + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"))+ coord_cartesian(xlim = c(0,50000), ylim = c(20,40))
# 
# plot1 + plot2 + plot3
# 
# 
# #ggsave(filename = "Figures/ResponseToCumulativeIYD.svg", width = 24, height = 24, units = "cm", dpi = 600)
# ggsave(filename = "Figures/ResponseToIYD_high_20230503.jpg", width = 34, height = 24, units = "cm", dpi = 600)
# 
# 
# 
# #now DFP
# m1 <- glmer(DFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), data = t1_sum ) #%>% filter(iyd < 20000)
# summary(m1)
# 
# m2 <- glmer(DFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), data = t2_sum) #%>% filter(iyd < 20000)
# summary(m2)
# 
# m3 <- glmer(DFP ~ scale(iyd) + (1|id_yr) + (1|AID) + (1|Year), data = t3_sum) #%>% filter(iyd < 20000)
# summary(m3)
# 
# 
# r.squaredGLMM(m1)
# r.squaredGLMM(m2)
# r.squaredGLMM(m3)
# 
# # test relationship between timestopped and absDFP
# 
# 
# newdata <- data.frame(DFP = seq(min(na.omit(t1_sum$DFP)), max(na.omit(t1_sum$DFP)), length = 50), iyd = seq(min(na.omit(t1_sum$iyd)), max(na.omit(t1_sum$iyd)), length = 50), id_yr = "108_2018", AID = "108", Year = 2016)
# 
# p1 <- PredictGLMER(m1, newdata, se.fit = T, seMultiplier = 1.96)
# p2 <- PredictGLMER(m2, newdata, se.fit = T, seMultiplier = 1.96)
# p3 <- PredictGLMER(m3, newdata, se.fit = T, seMultiplier = 1.96)
# 
# 
# 
# 
# plot1 <- ggplot() + geom_ribbon(data = cbind(p1, newdata), aes(x = iyd, ymin = yminus, ymax = yplus), fill = "#008B8B", alpha = .75) + geom_line(data = cbind(p1, newdata), aes(x = iyd, y=y), color = "black") + theme_classic() + labs(y = "Days from Peak Green-up", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + coord_cartesian(xlim = c(0,50000), ylim = c(15,45)) 
# 
# plot2 <- ggplot() + geom_ribbon(data = cbind(p2, newdata), aes(x = iyd, ymin = (yminus), ymax = (yplus)), fill = "#008B8B", alpha = .3) + geom_line(data = cbind(p2, newdata), aes(x = iyd, y=(y)), color = "black") + theme_classic() + labs(y = "Days from Peak Green-up", x = bquote('Inter-year Distance'~(Fidelity ^-1))) + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18"))+ coord_cartesian(xlim = c(0,50000), ylim = c(15,45))
# 
# 
# plot1 + plot2 + plot3
# 
# ggsave(filename = "Figures/ResponseToIYD_high_DFP_20230503.jpg", width = 34, height = 24, units = "cm", dpi = 600)
# 
# 
# #---------------------------------#
# # summary stats ####
# 
# summary(t1_sum$iyd/1000)
# summary(t2_sum$iyd/1000)
# summary(t3_sum$iyd/1000)
# 
# 
# 
# ggplot() + geom_density(data = t1_sum, aes(x = iyd/1000), bounds = c(0, 250), fill = "blue", alpha = .2) + geom_density(data = t2_sum, aes(x = iyd/1000), bounds = c(0, 250), fill = "red", alpha = .2) + geom_density(data = t3_sum, aes(x = iyd/1000), bounds = c(0, 250), fill = "green", alpha = .2) + coord_cartesian(xlim = c(-2,20)) + scale_x_continuous(limits = c(0,22))


#---------------------------------#
# start the gams ####

mod.absDFP1 <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5), data = t1_sum %>% filter(iyd < 15000), method = "GCV", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log")) # + km + TimeStopped
mod.absDFP2 <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) , data = t2_sum%>% filter(iyd < 15000), method = "GCV", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))
mod.absDFP3 <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) , data = t3_sum%>% filter(iyd < 15000), method = "GCV", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))


newdata = data.frame(km = mean(t1_sum$km), iyd = seq(1,15000,length = 100), TimeStopped = mean(t1_sum$TimeStopped), AID = "108", Year = 2017, id_yr = "108_2017")
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

ggplot(new) + geom_ribbon(aes(x = iyd/1000, ymax = exp(prediction) + exp(se), ymin = exp(prediction) - exp(se), group = factor(Year)), fill = "#008B8B", alpha = .33) + geom_line(aes(x = iyd/1000, y = exp(prediction), linetype = factor(Year)), size = .8) + labs(y = "|Days from Peak Green-up|", x = bquote('Inter-year Distance (km)'~( Fidelity ^-1)), linetype = "") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + scale_linetype_manual(values = c("solid", "dashed", "dotted"), label = c("t-1", "t-2", "t-3")) + coord_cartesian(ylim = c(15,45), xlim = c(0,15)) + scale_y_continuous(limits= c(14.8,50.2), breaks = seq(15,50,5)) + scale_x_continuous(limits= c(0,15), breaks = seq(0,15,3))
 
ggsave(filename = "Figures/ResponseToIYD_gam_absDFP_20230511.jpg", width = 30, height = 30, units = "cm", dpi = 600)








mod.DFP1 <- mgcv::gamm(DFP ~ s(iyd, bs = "cs", k = 5)  + km + TimeStopped, data = t1_sum, method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1))
mod.DFP2 <- mgcv::gamm(DFP ~ s(iyd, bs = "cs", k = 5)  + km + TimeStopped, data = t2_sum, method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1))
mod.DFP3 <- mgcv::gamm(DFP ~ s(iyd, bs = "cs", k = 5)  + km + TimeStopped, data = t3_sum, method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1))


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

ggplot(new) + geom_ribbon(aes(x = iyd, ymax = (prediction) + (se), ymin = (prediction) - (se), group = factor(Year)), fill = "#008B8B", alpha = .33) + geom_line(aes(x = iyd, y = (prediction), linetype = factor(Year)), size = .8) + labs(y = "Days from Peak Green-up", x = bquote('Inter-year Distance (m)'~( Fidelity ^-1)), linetype = "") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + scale_linetype_manual(values = c("solid", "dashed", "dotted"), label = c("t-1", "t-2", "t-3")) + coord_cartesian(ylim = c(5,45), xlim = c(0,50000)) + scale_y_continuous(limits= c(4.8,50.2), breaks = seq(5,50,5)) + scale_x_continuous(limits= c(0,50000), breaks = seq(0,50000,10000))

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


mod.absDFP1.alt <- mgcv::gam(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t1_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML", family = poisson(link = "log"))
mod.absDFP2.alt <- mgcv::gam(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t2_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML", family = poisson(link = "log"))
mod.absDFP3.alt <- mgcv::gam(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5) + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t3_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML",  family = poisson(link = "log"))

summary(mod.absDFP1.alt)$dev.expl
summary(mod.absDFP2.alt)$dev.expl
summary(mod.absDFP3.alt)$dev.expl

#without iyd
mod.absDFP1.altx <- mgcv::gam(round(absDFP,0) ~  s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t1_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML", family = poisson(link = "log"))
mod.absDFP2.altx <- mgcv::gam(round(absDFP,0) ~   s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t2_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML", family = poisson(link = "log"))
mod.absDFP3.altx <- mgcv::gam(round(absDFP,0) ~   s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t3_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML",  family = poisson(link = "log"))

summary(mod.absDFP1.altx)$dev.expl
summary(mod.absDFP2.altx)$dev.expl
summary(mod.absDFP3.altx)$dev.expl

# mod.DFP1.alt <- mgcv::gam(DFP ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t1_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML")
# mod.DFP2.alt <- mgcv::gam(DFP ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t2_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML")
# mod.DFP3.alt <- mgcv::gam(DFP ~ s(iyd, bs = "cs", k = 5) + km + TimeStopped + s(fid, bs = "re") + s(fy, bs = "re") + s(faid, bs = "re"), data = t3_sum %>% mutate(fid = factor(id_yr), fy = factor(Year), faid = factor(AID)), method = "REML")
# 
# 
# summary(mod.DFP1.alt)$dev.expl
# summary(mod.DFP2.alt)$dev.expl
# summary(mod.DFP3.alt)$dev.expl



#------------------------#
# benefit of deviation ####

devtest <- t1_sum %>% filter(km < 241) #t1 %>% filter(km_mark_1 <= 240) %>% mutate(km = round(km_mark_1, 0)) %>% group_by(id_yr, km) %>% summarise(iyd = mean(iyd_1), mismatch = mean(absDFP) , mis.sign = mean(DFP), AID = unique(AID), Year = unique(Year))
# devtest <- devtest %>% filter(iyd)

devtest_sum <- devtest %>% group_by(id_yr) %>% summarise(meanIYD = mean(iyd, na.rm = T))

devtest <- devtest %>% left_join(devtest_sum %>% st_drop_geometry(), by = "id_yr")

devtest <- devtest %>% drop_na(iyd, dev, meanIYD) %>% ungroup() %>% mutate(fidc = cut(meanIYD, breaks = quantile(meanIYD, c(.25,.5,.75,1.00)),include.lowest = TRUE, labels = FALSE), devc = cut(dev, breaks = quantile(iyd, c(.25,.5,.75,1.00)),include.lowest = TRUE, labels = FALSE))





mod.flex <- mgcv::gamm(absDFP ~ s(km, by = factor(fidc), bs = "cs", k = 5), method = "GCV", data = devtest, family = quasipoisson, random = list(id_yr=~1, Year = ~1, AID=~1))


mod.sign.flex <- mgcv::gamm(DFP ~ s(km, by = factor(fidc), bs = "cs", k = 10), method = "REML",data = devtest, random = list(id_yr=~1, Year = ~1, AID=~1))

mod.dev <- mgcv::gamm(absDFP ~ s(km, by = factor(devc), bs = "cs", k = 5), method = "GCV", data = devtest, family = quasipoisson, random = list(id_yr=~1, Year = ~1, AID=~1))
mod.sign.dev <- mgcv::gamm(DFP ~ s(km, by = factor(devc), bs = "cs", k = 5), method = "REML",data = devtest, random = list(id_yr=~1, Year = ~1, AID=~1))


summary(mod.sign.dev$gam)

#plot.gam(mod.dev$gam)

newdata <- data.frame(expand.grid(km = seq(0,240,by = 1), fidc = c(1,2,3), devc = c(1,2,3), Year = 2017, id_yr = "255_2018", AID = "255"))

predflex <- gammit::predict_gamm(mod.flex$gam , newdata = newdata, se = T)
predflexvsign <- gammit::predict_gamm(mod.sign.flex$gam , newdata = newdata, se = T)

preddev <- gammit::predict_gamm(mod.dev$gam , newdata = newdata, se = T)
preddevvsign <- gammit::predict_gamm(mod.sign.dev$gam , newdata = newdata, se = T)



#ggplot(cbind(newdata, predflex)) + geom_ribbon(aes(x = km, ymin = exp(prediction - se), ymax = exp(prediction + se), group = factor(fidc)), fill = "darkcyan", alpha = .33) + geom_line(aes(x = km, y = exp(prediction), group = factor(fidc), linetype = factor(fidc)), linewidth = 1.4, color = "darkcyan") + scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c("High", "Moderate", "Low")) + labs(x = "Route Distance (km)", y = "|Days from Peak Green-up|", linetype = "Fidelity") + theme_classic()

p1 <- ggplot(cbind(newdata, predflexvsign)) + geom_ribbon(aes(x = km, ymin = (prediction - se), ymax = (prediction + se), group = factor(fidc)), fill = "darkcyan", alpha = .33) + geom_line(aes(x = km, y = (prediction), group = factor(fidc), linetype = factor(fidc)), linewidth = 1.4, color = "darkcyan") + scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c("High", "Moderate", "Low")) + labs(x = "Distance from Winter Range (km)", y = "Days from Peak Green-up", linetype = "Fidelity Score", title = "") + theme_classic() + theme(axis.title.x = element_text(size = 40,color = "grey18"), axis.title.y = element_text(size = 40,color = "grey18"),axis.text.x = element_text(size = 30,color = "grey18"),axis.text.y = element_text(size = 30,color = "grey18"), legend.text = element_text(size = 30,color = "grey18"), legend.title = element_text(size = 40,color = "grey18"),legend.key.size = unit(1, 'cm'), legend.key.height = unit(1, 'cm'), legend.key.width = unit(2, 'cm')) + coord_cartesian(ylim = c(-20,63), xlim = c(0,250))+ scale_y_continuous(expand = c(0,0), limits = c(-30,65.2), breaks = seq(-20,60,10)) + scale_x_continuous(expand = c(0,0), limits = c(0,260), breaks = seq(0,240,40)) + geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1.2); p1

ggsave(filename = "Figures/OptimalStops_20230516.jpg", width = 35, height = 24, units = "cm", dpi = 600)

ggplot(cbind(newdata, predflex)) + geom_ribbon(aes(x = km, ymin = exp(prediction - se), ymax = exp(prediction + se), group = factor(fidc)), fill = "darkcyan", alpha = .33) + geom_line(aes(x = km, y = exp(prediction), group = factor(fidc), linetype = factor(fidc)), linewidth = 1.4, color = "darkcyan") + scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c("High", "Moderate", "Low")) + labs(x = "Route Distance (km)", y = "Days from Peak Green-up", linetype = "Fidelity", title = "Stops") + theme_classic()




#ggplot(cbind(newdata, preddev)) + geom_ribbon(aes(x = km, ymin = exp(prediction - se), ymax = exp(prediction + se), group = factor(devc)), fill = "darkcyan", alpha = .33) + geom_line(aes(x = km, y = exp(prediction), group = factor(devc), linetype = factor(devc)), linewidth = 1.4, color = "darkcyan") + scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c("High", "Moderate", "Low")) + labs(x = "Route Distance (km)", y = "|Days from Peak Green-up|", linetype = "Fidelity") + theme_classic()

p2 <- ggplot(cbind(newdata, preddevvsign)) + geom_ribbon(aes(x = km, ymin = (prediction - se), ymax = (prediction + se), group = factor(devc)), fill = "darkcyan", alpha = .33) + geom_line(aes(x = km, y = (prediction), group = factor(devc), linetype = factor(devc)), linewidth = 1.4, color = "darkcyan") + scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c("High", "Moderate", "Low")) + labs(x = "Route Distance (km)", y = "Days from Peak Green-up", linetype = "Fidelity", title = "Route") + theme_classic() + theme(axis.title.x = element_text(size = 40,color = "grey18"), axis.title.y = element_text(size = 40,color = "grey18"),axis.text.x = element_text(size = 30,color = "grey18"),axis.text.y = element_text(size = 30,color = "grey18"),legend.text = element_text(size = 30,color = "grey18"), legend.title = element_text(size = 40,color = "grey18"),legend.key.size = unit(2, 'cm'), legend.key.height = unit(1, 'cm'), legend.key.width = unit(2, 'cm')) + coord_cartesian(ylim = c(-20,63), xlim = c(0,240))+ scale_y_continuous(expand = c(0,0), limits = c(-30,65.2), breaks = seq(-20,60,10)) + scale_x_continuous(expand = c(0,0), limits = c(0,255), breaks = seq(0,240,40)) + geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1.2)



ggsave(filename = "Figures/OptimalDeviations_20230513.jpg", width = 40, height = 30, units = "cm", dpi = 600)



p1+p2

mod.alt <- mgcv::gam(round(mismatch,0) ~ s(km, by = factor(fidc), bs = "cs", k = 5) + s(factor(devtest$id_yr), bs= "re") + s(factor(devtest$AID), bs = "re") + s(factor(devtest$Year), bs = "re"), method= "REML", data = devtest, family = quasipoisson())
summary(mod.alt)$dev.expl #random = list(id_yr=~1, Year = ~1, AID=~1)

mod.alts <- mgcv::gam(mis.sign ~ s(km, by = factor(fidc), bs = "cs", k = 5) + s(factor(devtest$id_yr), bs= "re") + s(factor(devtest$AID), bs = "re") + s(factor(devtest$Year), bs = "re"), method= "REML", data = devtest)
summary(mod.alts)$dev.expl

#-------------------------------#
# varied by greenscape ####

trans.arcsine <- function(x){
  asin(sign(x) * sqrt(abs(x)))
}


#chose two year since has most deviance explained
green <- RDH_14t22_greenscapes %>% rename(id_yr = AY_ID, AID = AnimalID)
head(green)
m1 <- glmer(IRG_120 ~ greenUpDur + (1+greenUpDur|Year), data = green)
m2 <- glmer(log(DFP_120) ~ greenUpDur + (1+greenUpDur|Year), data = green)
m3 <- glmer(log(DFP_120) ~ trans.arcsine(svSlope) + (1+svSlope|Year), data = green, na.action = "na.omit")
m4 <- glmer(log(DFP_120) ~ springScale + (1+springScale|Year), data = green)


summary(m1)
summary(m2)
summary(m3)
summary(m4)

newdatagud <- data.frame(expand.grid(Year = seq(2014, 2022, by = 1), greenUpDur = seq(min(green$greenUpDur),max(green$greenUpDur),length = 100)))
newdatasps <- data.frame(expand.grid(Year = seq(2014, 2022, by = 1), springScale = seq(min(RDH_14t22_greenscapes$springScale),max(RDH_14t22_greenscapes$springScale),length = 100)))
newdatasvs <- data.frame(expand.grid(Year = seq(2014, 2022, by = 1), svSlope = seq(min(RDH_14t22_greenscapes$svSlope),max(RDH_14t22_greenscapes$svSlope),length = 100)))
newdatasvs$svSlope <- trans.arcsine(newdatasvs$svSlope)
                                 

newdatagud$pred <-  predict(m1, newdatagud, re.form = ~(1+greenUpDur|Year))
newdatagud$predd <-  predict(m2, newdatagud, re.form = ~(1+greenUpDur|Year))
newdatasps$pred <- predict(m4, newdatasps, re.form = ~(1+springScale|Year))
newdatasvs$pred <- predict(m3, newdatasvs, re.form = ~(1+svSlope|Year))

p1 <- ggplot(green) + geom_point(aes(y = IRG_120, x = greenUpDur, group = factor(Year), color = factor(Year))) + geom_line(aes(x = greenUpDur, y = pred, group = factor(Year), color = factor(Year)), linewidth = 1.2, newdatagud) + geom_line(aes(x = greenUpDur, y = mean), color = "black", linewidth = 1.2, newdatagud %>% group_by(greenUpDur) %>% summarise(mean = mean(pred))) + theme_classic() + theme(legend.position = "none");p1

p2 <- ggplot(green) + geom_point(aes(y = log(DFP_120), x = greenUpDur, group = factor(Year), color = factor(Year))) + geom_line(aes(x = greenUpDur, y = predd, group = factor(Year), color = factor(Year)), linewidth = 1.2, newdatagud) + geom_line(aes(x = greenUpDur, y = mean), color = "black", linewidth = 1.2, newdatagud %>% group_by(greenUpDur) %>% summarise(mean = mean(predd))) + theme_classic() + theme(legend.position = "none")+ scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1)) + scale_x_continuous(limits = c(50,150), breaks = seq(50,150,25)); p2

p3 <- ggplot(green) + geom_point(aes(y = log(DFP_120), x = trans.arcsine(svSlope), group = factor(Year), color = factor(Year))) + geom_line(aes(x = trans.arcsine(svSlope), y = pred, group = factor(Year), color = factor(Year)), linewidth = 1.2, newdatasvs) + geom_line(aes(x = trans.arcsine(svSlope), y = mean), color = "black", linewidth = 1.2, data =  newdatasvs %>% group_by(svSlope) %>% summarise(mean = mean(pred))) + theme_classic() + theme(legend.position = "none") + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1))  + scale_x_continuous(limits = c(-.5,1.3), breaks = seq(-.5,1.3,.3)); p3

p4 <- ggplot(green) + geom_point(aes(y = log(DFP_120), x = springScale, group = factor(Year), color = factor(Year))) + geom_line(aes(x = springScale, y = pred, group = factor(Year), color = factor(Year)), linewidth = 1.2, newdatasps) + geom_line(aes(x = springScale, y = mean), color = "black", linewidth = 1.2, newdatasps %>% group_by(springScale) %>% summarise(mean = mean(pred))) + theme_classic() + theme(legend.position = "none") + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1))  + scale_x_continuous(limits = c(9.95, 32.05), breaks = seq(10,32,1)); p4

p1 + p2 + p4 + p3




#how does fidelity mediate
green <- RDH_14t22_greenscapes %>% rename(id_yr = AY_ID, AID = AnimalID)
head(green)

t1_mig <- t1_sum %>% group_by(id_yr) %>% summarise(fid = mean(na.omit(iyd)), TimeOnStop = sum(TimeStopped), TimeStopAvg = mean(TimeStopped)) %>% drop_na(fid)
t2_mig <- t2_sum %>% group_by(id_yr) %>% summarise(fid = mean(na.omit(iyd)), TimeOnStop = sum(TimeStopped), TimeStopAvg = mean(TimeStopped)) %>% drop_na(fid)
t3_mig <- t3_sum %>% group_by(id_yr) %>% summarise(fid = mean(na.omit(iyd)), TimeOnStop = sum(TimeStopped), TimeStopAvg = mean(TimeStopped)) %>% drop_na(fid)

t1g <- t1_mig %>% left_join(green, by = c("id_yr")) %>% ungroup() %>% mutate(iyds = scale(fid))
t2g <- t2_mig %>% left_join(green, by = c("id_yr")) %>% ungroup() %>% mutate(iyds = scale(fid))
t3g <- t3_mig %>% left_join(green, by = c("id_yr")) %>% ungroup() %>% mutate(iyds = scale(fid))

#t2g <- t2g %>% mutate(GUDs = scale(greenUpDur), SPSs = scale(springScale), SVSs = scale(svSlope), GUDc = ifelse(GUDs <= 0, "short", "long"), SPSc = ifelse(SPSs <= 0,  "fast", "slow" ), SVSc = ifelse(SVSs <= 0,"rand", "ord")) 

t1gc <- t1g %>% st_drop_geometry() %>% mutate(fidc = cut(fid, breaks = quantile(fid, c(0,.5,1.00)),include.lowest = TRUE, labels = FALSE)) %>% select(-c(SVmax, midpoint))
t2gc <- t2g %>% st_drop_geometry() %>% mutate(fidc = cut(fid, breaks = quantile(fid, c(0,.5,1.00)),include.lowest = TRUE, labels = FALSE)) %>% select(-c(SVmax, midpoint))
t3gc <- t3g %>% st_drop_geometry() %>% mutate(fidc = cut(fid, breaks = quantile(fid, c(0,.5,1.00)),include.lowest = TRUE, labels = FALSE)) %>% select(-c(SVmax, midpoint))

#


# mod.absDFP.GUDc <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5, by = factor(GUDc)) + km + TimeStopped, data = t2g %>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))
# mod.absDFP.SVSc <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5, by = factor(SVSc))  + km + TimeStopped, data = t2g%>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))
# mod.absDFP.SPSc <- mgcv::gamm(round(absDFP,0) ~ s(iyd, bs = "cs", k = 5, by = factor(SPSc))  + km + TimeStopped, data = t2g%>% filter(iyd < 30000), method = "REML", random = list(id_yr=~1, Year = ~1, AID=~1), family = poisson(link = "log"))

t1gc$Year <- as.factor(t1gc$Year)

t1gc <- t1gc %>% mutate(tasv = trans.arcsine(svSlope))

mod.IRG.GUD <- mgcv::gamm(IRG_120 ~ s(greenUpDur, bs = "cs", k = 5, by = factor(fidc)) + greenUpDur, data = t1gc, method = "REML", random = list(Year =  ~1+greenUpDur))
mod.DFP.GUD <- mgcv::gamm(log(DFP_120) ~ s(greenUpDur, bs = "cs", k = 5, by = factor(fidc)) + greenUpDur, data = t1gc, method = "REML", random = list(Year = ~1+greenUpDur))
mod.DFP.SPS <- mgcv::gamm(log(DFP_120) ~ s(springScale, bs = "cs", k = 5, by = factor(fidc)) + springScale, data = t1gc , method = "REML", random = list(Year = ~1+springScale))
mod.DFP.SVS <- mgcv::gamm(log(DFP_120) ~ s(tasv, bs = "cs", k = 5, by = factor(fidc)) + tasv, data = t1gc, method = "REML", random = list(Year = ~1+tasv))

summary(mod.IRG.GUD$gam)
summary(mod.DFP.GUD$gam)
summary(mod.DFP.SPS$gam)
summary(mod.DFP.SVS$gam)

summary(mod.IRG.GUD$lme)
summary(mod.DFP.GUD$lme)
summary(mod.DFP.SPS$lme)
summary(mod.DFP.SVS$lme)

years <- data.frame(years = unique(t1g$Year))
years <- na.omit(years)

newdata.gud <- data.frame(expand.grid(fidc = factor(c(1,2)), AID = "108", Year = years$years, id_yr = "108_2017", greenUpDur = seq(min(na.omit(t1g$greenUpDur)), max(na.omit(t1g$greenUpDur)), length = 50)))
newdata.svs <- data.frame(expand.grid(fidc = factor(c(1,2)), AID = "108", Year = years$years, id_yr = "108_2017", tasv = seq(min(na.omit(t1gc$tasv)), max(na.omit(t1gc$tasv)), length = 50)))
newdata.sps = data.frame(expand.grid(fidc = factor(c(1,2)), AID = "108", Year = years$years, id_yr = "108_2017", springScale = seq(min(na.omit(t1g$springScale)), max(na.omit(t1g$springScale)), length = 50)))

extract_ranef(mod.IRG.GUD$gam)

predIRGGUD <- gammit::predict_gamm(mod.IRG.GUD$gam , newdata = newdata.gud, se = T); predIRGGUD
predGUD <- gammit::predict_gamm(mod.DFP.GUD$gam , newdata = newdata.gud, random = T,se = T)
predSVS <- gammit::predict_gamm(mod.DFP.SVS$gam , newdata = newdata.svs,  se = T)
predSPS <- gammit::predict_gamm(mod.DFP.SPS$gam , newdata = newdata.sps, random = T,se = T)

newdatapredGUD <- cbind(newdata.gud, predGUD)
newdatapredSVS <- cbind(newdata.svs, predSVS)
newdatapredSPS <- cbind(newdata.sps, predSPS)

irggud <- ggplot(cbind(newdata.gud, predIRGGUD)) + geom_point(aes(x = greenUpDur, y = (IRG_120)), size = 3.5, alpha = .15, shape = 21, color = "black", fill = "grey30", data = t1g) + geom_line(aes(x = greenUpDur, y = (prediction), color = factor(fidc)), size = 1.8) + labs(y = "Instan. Rate of Green-up", title = "", x = "Green-up duration (d)", color = "Fidelity Quartile", fill = "Fidelity Quartile") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 17,color = "grey18"), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(.5,.9), breaks = seq(.5,.9,.1)) + scale_x_continuous(limits = c(50,150), breaks = seq(50,150,25)); irggud

gud <- ggplot(cbind(newdata.gud, predGUD)) + geom_point(aes(x = greenUpDur, y = log(DFP_120)), size = 3.5, alpha = .15, shape = 21, color = "black", fill = "grey30", data = t1g) + geom_line(aes(x = greenUpDur, y = (prediction), color = factor(fidc)), size = 1.8) + labs(y = "ln(Days from Peak Green-up)", title = "", x = "Green-up duration (d)", color = "Fidelity Quartile", fill = "Fidelity Quartile") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 17,color = "grey18"), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1)) + scale_x_continuous(limits = c(50,150), breaks = seq(50,150,25)); gud

svs <- ggplot(cbind(newdata.svs, predSVS)) + geom_point(aes(x = tasv, y = log(DFP_120)), fill = "grey30", size = 3.5, alpha = .15, shape = 21, color = "black", data = t1gc) + geom_line(aes(x = tasv, y = (prediction), color = factor(fidc)), size = 1.8) + labs(y = "ln(Days from Peak Green-up)", title = "", x = "arcsin(Green-up order)", color = "Fidelity Quartile", fill = "Fidelity Quartile") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 17,color = "grey18"), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18"))  + scale_fill_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1))  + scale_x_continuous(limits = c(.2,.71), breaks = seq(.2,.7,.1)); svs

sps <- ggplot(cbind(newdata.sps, predSPS)) + geom_point(aes(x = springScale, y = log(DFP_120)), fill = "grey30", size = 3.5, alpha = .15, shape = 21, color = "black", data = t1g) + geom_line(aes(x = springScale, y = (prediction), color = factor(fidc)), size = 1.8) + labs(y = "ln(Days from Peak Green-up)", x =  bquote('Spring scale'~('green-up rate' ^-1)), color = bquote('Morrison Fidelity'~('IYD')), fill = bquote('Morrison Fidelity'~('IYD'))) + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 17,color = "grey18"), title = element_text(size = 20,color = "grey18"), legend.position = "none") + scale_fill_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1))  + scale_x_continuous(limits = c(14.95, 22.05), breaks = seq(15,22,1)); sps


irggud + gud + sps + svs +plot_layout(nrow=2)

ggsave(filename = "Figures/GreenscapeMetricsByAvg_20230511.jpg", width = 48, height = 36, units = "cm", dpi = 600)




# recreate the yearly plots, but with gam
# 
# mod.IRG.GUDY <- mgcv::gam(IRG_120 ~ greenUpDur + s(greenUpDur, by = Year, bs = "re"), data = green, method = "REML")
# mod.DFP.GUDY <- mgcv::gam(log(DFP_120) ~ greenUpDur + s(greenUpDur, Year, bs = "re"), data = green, method = "REML")
# mod.DFP.SPSY <- mgcv::gam(log(DFP_120) ~ springScale + s(springScale, Year, bs = "re"), data = green, method = "REML")
# mod.DFP.SVSY <- mgcv::gam(log(DFP_120) ~ trans.arcsine(green$svSlope) + s(trans.arcsine(green$svSlope), Year, bs = "re"), data = green, method = "REML")
# 
# 
# years <- data.frame(years = unique(green$Year))
# years <- na.omit(years)
# 
# newdata.gud <- data.frame(expand.grid(AID = "108", Year = years$years, id_yr = "108_2017", greenUpDur = seq(min(na.omit(green$greenUpDur)), max(na.omit(green$greenUpDur)), length = 50)))
# newdata.svs <- data.frame(expand.grid(AID = "108", Year = years$years, id_yr = "108_2017", svSlope = seq(min(na.omit(green$svSlope)), max(na.omit(green$svSlope)), length = 50)))
# newdata.sps = data.frame(expand.grid(AID = "108", Year = years$years, id_yr = "108_2017", springScale = seq(min(na.omit(green$springScale)), max(na.omit(green$springScale)), length = 50)))
# 
# predIRGGUD <- gammit::predict_gamm(mod.IRG.gudY ,  newdata = newdata.gud, random = T, se = T); predIRGGUD
# predGUD <- gammit::predict_gamm(model = mod.DFP.GUDY , newdata = newdata.gud, random = T,se = T)
# predSVS <- gammit::predict_gamm(mod.DFP.SVSY , newdata = newdata.svs, random = T,  se = T)
# predSPS <- gammit::predict_gamm(mod.DFP.SPSY , newdata = newdata.sps, random = T, random = T,se = T)
# 
# newdatapredGUD <- cbind(newdata.gud, predGUD)
# newdatapredSVS <- cbind(newdata.svs, predSVS)
# newdatapredSPS <- cbind(newdata.sps, predSPS)
# 
# ggplot(cbind(newdata.gud, predIRGGUD)) + geom_point(aes(x = greenUpDur, y = (IRG_120)), size = 3.5, alpha = .15, shape = 21, color = "black", fill = "grey30", data = green) + geom_line(aes(x = greenUpDur, y = (prediction), group = factor(Year), color = factor(Year)), position = position_jitter(), size = 1.8) 
# 
# ggplot(cbind(newdata.gud, predGUD)) + geom_point(aes(x = greenUpDur, y = log(DFP_120)), size = 3.5, alpha = .15, shape = 21, color = "black", fill = "grey30", data = green) + geom_line(aes(x = greenUpDur, y = (prediction), group = factor(Year), color = factor(Year)), position = position_jitter(), size = 1.8) 
# 
# # 
# # + labs(y = "ln(|Days from Peak Green-up|)", title = "", x = "Green-up duration (d)", color = "Fidelity Quartile", fill = "Fidelity Quartile") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 17,color = "grey18"), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(.5,.9), breaks = seq(.5,.9,.1)) + scale_x_continuous(limits = c(50,150), breaks = seq(50,150,25)); irggud
# # 
# # gud <- ggplot(cbind(newdata.gud, predGUD)) + geom_point(aes(x = greenUpDur, y = log(DFP_120)), size = 3.5, alpha = .15, shape = 21, color = "black", fill = "grey30", data = t1g) + geom_line(aes(x = greenUpDur, y = (prediction), color = factor(fidc)), size = 1.8) + labs(y = "ln(|Days from Peak Green-up|)", title = "", x = "Green-up duration (d)", color = "Fidelity Quartile", fill = "Fidelity Quartile") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_text(size = 17,color = "grey18"), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1)) + scale_x_continuous(limits = c(50,150), breaks = seq(50,150,25)); gud
# 
# svs <- ggplot(cbind(newdata.svs, predSVS)) + geom_point(aes(x = trans.arcsine(svSlope), y = log(DFP_120)), fill = "grey30", size = 3.5, alpha = .15, shape = 21, color = "black", data = t1g) + geom_line(aes(x = trans.arcsine(svSlope), y = (prediction), color = factor(fidc)), size = 1.8) + labs(y = "", title = "", x = "arcsin(Green-up order)", color = "Fidelity Quartile", fill = "Fidelity Quartile") + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_blank(), legend.position = "none",legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c("#2F4F4F","#0000FF", "#87CEFA"), label = c("High", "Mod", "Low")) + scale_color_manual(values = c("#2F4F4F", "#0000FF", "#87CEFA"), label = c("High", "Mod", "Low")) + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1))  + scale_x_continuous(limits = c(.2,.71), breaks = seq(.2,.7,.1)); svs
# 
# sps <- ggplot(cbind(newdata.sps, predSPS)) + geom_point(aes(x = springScale, y = log(DFP_120)), fill = "grey30", size = 3.5, alpha = .15, shape = 21, color = "black", data = t1g) + geom_line(aes(x = springScale, y = (prediction), color = factor(fidc)), size = 1.8) + labs(y = "", x =  bquote('Spring scale'~('green-up rate' ^-1)), color = bquote('Morrison Fidelity'~('IYD')), fill = bquote('Morrison Fidelity'~('IYD'))) + theme_classic()  + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 17,color = "grey18"),axis.text.y = element_blank(),legend.text = element_text(size = 16,color = "grey18"), title = element_text(size = 20,color = "grey18")) + scale_fill_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_color_manual(values = c("#2F4F4F", "#87CEFA"), label = c("Above Average", "Below Average")) + scale_y_continuous(limits = c(2.875, 3.925), breaks = seq(2.9,3.9,.1))  + scale_x_continuous(limits = c(14.95, 22.05), breaks = seq(15,22,1)); sps


#--------------------------------#
# what does high Morrison mean ####

stoploc <- fid1 %>% mutate(stop.n.c = stop.n_1, km_mark = km_mark_1, TimeStopped = TimeStopped_1, Year = curr_Y) %>% group_by(stop.n.c) %>% mutate(km = round(mean(km_mark),0)) %>% select(id_yr, stop.n.c, km, TimeStopped, Year, prev_Y_1) %>% group_by(stop.n.c) %>% slice(1) %>% st_drop_geometry()

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



#-----------------------------------#
# how does IYD change over space ####


t12 <- t1_sum %>% left_join(t2_sum %>% mutate(id_yr = str_sub(stop.n.c,0,8)) %>% st_drop_geometry() %>% rename(iyd2 = iyd) %>% select(c(stop.n.c, iyd2, id_yr)), by = c("stop.n.c","id_yr"))

t123 <- t12 %>% left_join(t3_sum %>% mutate(id_yr = str_sub(stop.n.c,0,8))%>% st_drop_geometry() %>% rename(iyd3 = iyd) %>% select(stop.n.c, iyd3, id_yr), by = c("stop.n.c","id_yr"))

t123 <- t123 %>% st_drop_geometry() 


t123$meanfid <- rowMeans(t123[,c("iyd", "iyd2", "iyd3")], na.rm = TRUE)

t123 <- t123 %>% mutate(fY = factor(Year), fA = factor(AID), fI = factor(id_yr))


mod.km_fid <- mgcv::gam(round(meanfid,0) ~ s(km, bs = "cs", k = 5) + s(km, fY, bs = "re")  + s(fA, bs = "re")  + s(fI, bs = "re"), data = t123 %>% filter(km < 240), family = poisson(link = "log"))

summary(mod.km_fid)

plot.gam(mod.km_fid)

newdata <- data.frame(km = seq(0,240,1), fY = unique(t123$fY)[1], fA = unique(t123$fA)[1], fI = unique(t123$fI)[1])

predict <- predict_gamm(mod.km_fid, newdata, se = T); predict

intercept <- predict_gamm(mod.km_fid, newdata = data.frame(expand.grid(km = seq(0,240,1), fY = unique(t123$fY)[1], fA = unique(t123$fA), fI = unique(t123$fI)[1])), re_form = "s(fA)", se = F)

ggplot(cbind(newdata, predict)) + geom_ribbon(aes(x = km, ymin = exp(prediction-se), ymax = exp(prediction+se)), fill = "grey45", alpha = .2)+ geom_line(aes(x = km, y = exp(prediction)), color = "#008B8B", linewidth = 1.3) + theme_classic() + labs(x = "Route Distance (km)", y = "Fidelity to Previous Stops") + theme(axis.title.x = element_text(size = 20,color = "grey18"), axis.title.y = element_text(size = 20,color = "grey18"),axis.text.x = element_text(size = 13,color = "grey18"),axis.text.y = element_text(size = 13,color = "grey18"),legend.text = element_text(size = 16,color = "grey18")) + coord_cartesian(xlim = c(0,245) ,ylim = c(0,31000)) +  scale_x_continuous(expand = c(0,0), limits = c(-20,252), breaks = seq(0,240,40)) + scale_y_continuous(expand = c(0,0), limits = c(0,40000), breaks = seq(0,30000,5000))

ggsave(filename = "Figures/IYDAlongRoute_20230513.jpg", width = 34, height = 24, units = "cm", dpi = 600)


ggplot(cbind(newdata, intercept)) + geom_point(aes(x = 1, y = prediction, group = fA, color = fA))

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



# mod.DriveGlobal.Hi <- gam::gam(iyd_1 ~ km_mark_1 + absDFP + IntegratedNDVI + (1|id_yr) + (1|), family = poisson(link = "log"), data = testhi, na.action = "na.fail")




#---------------------------------#
#compare used and alternate####


