#Code for analyzing days from peak green-up on seasonal ranges and at the start and end of spring migration
#Written by Anna Ortega, Wyoming Cooperative Fish and Wildlife Research Unit
#December 28, 2021

#Load required libraries
library(cowplot)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggridges)
library(Rmisc)
library(rstatix)
library(stringr)

#Clean working environment
rm(list=ls())

#### IMPORT CLEAN GPS COLLAR DATA ####
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data/SpringMigrations")
load("RD2H_SpringMigrations_2011_thru_2020_WithTimingDistanceIRG_28DEC2021.RData")
head(r)

#Import mean DFP for winter ranges and summer ranges
#First, import mean date of peak IRG on winter range
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data/WinterRange/IRG")
winter<-read.csv("IRGWinterSeasonalRanges.csv",header=TRUE,sep=",")
colnames(winter)<-c("id_yr","MaxIRG_WinterRange")
winter$seasonalrange<-"Migration start"

mig.winter<-merge(r,winter,by=c("id_yr"))
mig.winter<-mig.winter[order(mig.winter$id_yr,mig.winter$date),]
mig.winter<-mig.winter[duplicated(mig.winter$id_yr)==FALSE,]  #gives first row of dataframe for each month
head(mig.winter)

#Determine DFP at the start of spring migration
mig.winter$DFP_Start<-mig.winter$SpringMSJul-mig.winter$MaxIRG_WinterRange

#Determine mean date of peak IRG on winter range for each year of the study
sum.1<-summarySE(mig.winter,measurevar = c("MaxIRG_WinterRange"),groupvars = c("year"))
sum.1$mean<-round(mean(sum.1$MaxIRG_WinterRange))
sum.1$dev<-round(sum.1$MaxIRG_WinterRange-sum.1$mean)
sum.1<-sum.1[c("year","dev")]

w<-merge(mig.winter,sum.1,by=c("year"))
w$abs.DFP.start<-abs(w$DFP_Start)
w<-w[c("year","id_yr","dev","abs.DFP.start")]
colnames(w)<-c("year","id_yr","dev.winter","abs.DFP.start")

mig.winter<-mig.winter[c(1,12,15,28)]
colnames(mig.winter)<-c("id_yr","DOY_Start","timing","DFP_Start")
mig.winter$Cat<-"Migration start"
range(mig.winter$DOY_Start)

#Animals start migration anywhere from 1-150 days apart
table(mig.winter$DOY)

#Second, import mean date of peak IRG on summer range
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data/SummerRange/IRG")
summer<-read.csv("IRGSummerSeasonalRanges.csv",header=TRUE,sep=",")
colnames(summer)<-c("id_yr","MaxIRG_SummerRange")
summer$seasonalrange<-"Migration end"

mig.summer<-merge(r,summer,by=c("id_yr"))
mig.summer<-mig.summer[order(mig.summer$id_yr,mig.summer$date),]
mig.summer<-mig.summer[duplicated(mig.summer$id_yr)==FALSE,]  #gives first row of dataframe for each month

#Determine DFP at the end of spring migration
mig.summer$DFP_End<-mig.summer$SpringMEJul-mig.summer$MaxIRG_SummerRange

#Determine mean date of peak IRG on summer range for each year of the study
sum.2<-summarySE(mig.summer,measurevar = c("MaxIRG_SummerRange"),groupvars = c("year"))
sum.2$mean<-round(mean(sum.2$MaxIRG_SummerRange))
sum.2$dev<-round(sum.2$MaxIRG_SummerRange-sum.2$mean)
sum.2<-sum.2[c("year","dev")]

s<-merge(mig.summer,sum.2,by=c("year"))
s$abs.DFP.end<-abs(s$DFP_End)
s<-s[c("year","id_yr","dev","abs.DFP.end")]
colnames(s)<-c("year","id_yr","dev.summer","abs.DFP.end")

t<-merge(w,s,by=c("id_yr","year"))

m1<-lm(abs.DFP.start~dev.winter,data=t)
summary(m1)

m2<-lm(abs.DFP.end~dev.winter,data=t)
summary(m2)

mig.summer<-mig.summer[c(1,13,15,28)]
colnames(mig.summer)<-c("id_yr","DOY_End","timing","DFP_End")
mig.summer$Cat<-"Migration end"
range(mig.winter$DOY_Start)

#Combine DFP and DOY at start and end of spring migration
#Paired sample design
colnames(mig.winter)<-c("id_yr","DOY","timing","DFP","Cat")
colnames(mig.summer)<-c("id_yr","DOY","timing","DFP","Cat")
head(mig.winter)
head(mig.summer)

all<-rbind(mig.winter,mig.summer)
head(all)
all$year<- as.numeric(as.character(as.factor(str_split_fixed(all$id_yr, "_", 2)[,2])))
all<-all[c("id_yr","year","timing","Cat","DOY","DFP")]

#Export data
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Metadata")
save(all,file="RD2H_DFP_DOY_AtStartAndEndOfSpringMigration_20FEB2023.RData")

#Briefly look at variation in date of peak IRG on summer range among years
summer$year<- as.factor(str_split_fixed(summer$id_yr, "_", 2)[,2])
winter$year<- as.factor(str_split_fixed(winter$id_yr, "_", 2)[,2])
head(mig.winter)

# setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Metadata")
# write.csv(mig.winter,"MismatchAtStartofSpringMigration.csv")
# write.csv(mig.summer,"MismatchAtEndofSpringMigration.csv")

#Ensure that levels of departure and start/end are correctly ordered
all$Cat<-factor(all$Cat, levels=c("Migration start","Migration end"))
all$Level<-ifelse(all$Cat=="Migration start","Start","End")
all$Cat<-factor(all$Cat, levels=c("Migration start","Migration end"))
all$Level<-factor(all$Level, levels=c("Start","End"))

all$timing<-factor(all$timing, levels=c("early","mid","late"))
levels(all$Cat)
levels(all$timing)

# #Determine absolute days from peak
# all$abs<-abs(all$DFP)
# mig.winter$abs<-abs(mig.winter$DFP)
# mig.summer$abs<-abs(mig.summer$DFP)

# #Analyze summary of DFP at the start and end of spring migration
# library(Rmisc)
# 
# sum.1 = summarySE(mig.winter,measurevar="DFP", groupvars=c("timing"))
# sum.1$true.ci<-sum.1$se*1.96
# sum.1
# 
# #Analyze summary of start and end of spring migration
# sum.2 = summarySE(mig.summer,measurevar="DFP", groupvars=c("timing"))
# sum.2$true.ci<-sum.2$se*1.96
# sum.2

#Add id and year column to each dataframe
mig.summer$id<- str_split_fixed(mig.summer$id_yr, "_", 2)[,1]
mig.summer$year<- str_split_fixed(mig.summer$id_yr, "_", 2)[,2]

mig.winter$id<- str_split_fixed(mig.winter$id_yr, "_", 2)[,1]
mig.winter$year<- str_split_fixed(mig.winter$id_yr, "_", 2)[,2]
head(mig.winter)
colnames(mig.winter)<-c("id_yr","DOY_Start","timing","DFP_Start","Cat","id","year")
colnames(mig.summer)<-c("id_yr","DOY_End","timing","DFP_End","Cat","id","year")

sw<-merge(mig.winter,mig.summer,by=c("id_yr"))
names(sw)

sw<-sw[c(1:4,8,10)]
colnames(sw)<-c("id_yr","DOY_Start","timing","DFP_Start","DOY_End","DFP_End")
head(sw)

#Are DFP at the start and end correlated?
m1<-lm(sw$DOY_Start~sw$DOY_End)
summary(m1)

head(sw)
sum.1<-summarySE(sw,measurevar = c("DFP_Start"),groupvars = c("timing"))
sum.1

#### TEST FOR DIFFERENCES IN MIGRATION TIMING AND DFP ####
library(ggpubr)

early<-subset(mig.summer,mig.summer$timing=="early")
m1<-mean(early$DOY_End)
m2<-median(early$DOY_End)

mid<-subset(mig.summer,mig.summer$timing=="mid")
m1<-mean(mid$DOY_End)
m2<-median(mid$DOY_End)

late<-subset(mig.summer,mig.summer$timing=="late")
m1<-mean(late$DOY_End)
m2<-median(late$DOY_End)

ggplot(data=late, aes(x = DOY_End, y = 1)) +
  geom_density_ridges(alpha = 0.60, scale = 500 ,size=0.5) + theme_ridges()+
  geom_vline(xintercept=m1)


#Conduct one-way multivariate analysis of variance to understand if DFP and DOY at the start
#and end of spring migration differ among early, mid, and late migrants
mod1 <- manova(cbind(DFP_Start, DFP_End,DOY_Start, DOY_End) ~ timing, data = sw)
summary(mod1)
#summary.aov(mod1)
summary(mod1)[[4]]

#Conduct Tukey HSD to further identify differences among early, mid, and late migrants
mod2<-aov(DFP_Start~timing, data = sw)
summary(mod2)
tukey.test <- TukeyHSD(mod2)
tukey.test$timing

mod3<-aov(DFP_End~timing, data = sw)
summary(mod3)
tukey.test <- TukeyHSD(mod3)
tukey.test$timing

mod4<-aov(DOY_Start~timing, data = sw)
summary(mod4)
tukey.test <- TukeyHSD(mod4)
tukey.test$timing

mod5<-aov(DOY_End~timing, data = sw)
summary(mod5)
tukey.test <- TukeyHSD(mod5)
tukey.test$timing

mig.winter$Cat<-factor(mig.winter$Cat, levels=c("Migration start","Migration end"))
mig.winter$timing<-factor(mig.winter$timing, levels=c("early","mid","late"))

mig.summer$Cat<-factor(mig.summer$Cat, levels=c("Migration start","Migration end"))
mig.summer$timing<-factor(mig.summer$timing, levels=c("early","mid","late"))

#### PLOT MIGRATING TIMING AND DFP AT START AND END OF SPRING MIGRATION ####
#Create density ridgeline plot showing DOY at the start of spring migration

#Color scheme:
#"#333399","#339966","#CC6600"

#If plotting for presentation
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/DissertationDefense/Graphs/Chapter1")
jpeg("MigrationTiming_EarlyMidLateEnd.jpeg", width = 4, height = 4, units = 'in', res = 350)

a1<-ggplot(data=mig.summer, aes(x = DOY_End, y = timing, colour=timing, fill=timing)) +
  geom_density_ridges(alpha = 0.60, scale = 500 ,size=0.5) + theme_ridges() +
  scale_fill_manual(values = c("#333399","#339966","#CC6600"),
                    name = "Depature from winter range",
                    labels=c("Early","Mid","Late")) +
  scale_color_cyclical(values = c("#333399","#339966","#CC6600"))+
  xlim(30,230)+
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5)))
  #ggtitle("a1")
a1

a1

d1 <- a1+xlab("Day of year") + ylab("Probability density")

e1<- d1 + theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.line.x = element_line(colour = "black"),
                             axis.line.y = element_line(colour = "black"))
e1
f1<- e1+ theme(axis.text.x=element_text(size=14,colour="black"),
               axis.text.y=element_blank(),
               axis.ticks.length.x = unit(.1, "cm"),
               axis.ticks.y = element_blank(),
               axis.title.x=element_text(size=14,colour="black",margin=margin(5,5,5,5)),
               axis.title.y=element_text(size=14,colour="black",margin=margin(5,5,5,5)),
               plot.title = element_text(size=14,face="bold",vjust = - 7,hjust=-0.09))
f1
g1<- f1 + theme (plot.background=element_blank(),
                 plot.margin = unit(c(.1,.1,.1,.1),"cm"),
                 legend.position="none",
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"))

g1

# h1<-g1+annotate('text', x= 35,y=490,
#                                 label="Start",
#                                 size=3)
# h1

dev.off()

#Create density ridgeline plot showing DOY at the end of spring migration
a2<-ggplot(data=mig.summer, aes(x = DOY_End, y = timing, colour=timing, fill=timing)) +
  geom_density_ridges(alpha = 0.60, scale = 500 ,size=0.5) + theme_ridges() +
  scale_fill_manual(values = c("#9966CC","#339966","#FF9933"),
                    name = "Depature from winter range",
                    labels=c("Early","Mid","Late")) +
  scale_color_cyclical(values = c("#9966CC","#339966","#FF9933"))+
  xlim(30,230)+
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5)))+
  ggtitle("a2")
a2

d2 <- a2+xlab("Day of year") + ylab("Probability density")
d2

e2<- d2 + theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.line.x = element_line(colour = "black"),
                             axis.line.y = element_line(colour = "black"))
e2
f2<- e2+ theme(axis.text.x=element_text(size=9, colour="black"),
               axis.text.y=element_blank(),
               axis.ticks.length.x = unit(.1, "cm"),
               axis.ticks.y = element_blank(),
               axis.title.x=element_text(size=9,margin=margin(5,5,5,5)),
               axis.title.y=element_text(size=9,margin=margin(5,5,5,5)),
               plot.title = element_text(size=9,face="bold",vjust = - 7,hjust=-0.09))

f2
g2<- f2 + theme (plot.background=element_blank(),
                 plot.margin = unit(c(.1,.1,.1,.1),"cm"),
                 legend.position="none",
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"))

g2

h2<-g2+annotate('text', x= 35,y=490,
                label="End",
                size=3)
h2

#Create density ridgeline plot showing DFP at the start of spring migration
a3<-ggplot(data=mig.winter, aes(x = DFP_Start, y = timing, colour=timing, fill=timing)) +
  geom_density_ridges(alpha = 0.60, scale = 500 ,size=0.5) + theme_ridges() +
  geom_vline(xintercept = 0, size=1)+
  scale_fill_manual(values = c("#9966CC","#339966","#FF9933"),
                    name = "Depature from winter range",
                    labels=c("Early","Mid","Late")) +
  scale_color_cyclical(values = c("#9966CC","#339966","#FF9933"))+
  xlim(-100,90)+
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5)))+
  ggtitle("b1")

a3

d3 <- a3+xlab("Days from peak IRG") + ylab("Density")
d3

e3<- d3 + theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.line.x = element_line(colour = "black"),
                             axis.line.y = element_line(colour = "black"))
e3
f3<- e3+ theme(axis.text.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.length.x = unit(.1, "cm"),
               axis.ticks.y = element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               plot.title = element_text(size=9,face="bold",vjust = - 7,hjust=-0.09))

g3<- f3 + theme (plot.background=element_blank(),
                 plot.margin = unit(c(.1,.1,.1,.1),"cm"),
                 legend.position="none",
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"))

g3

h3<-g3+annotate('text', x= -95,y=490,
                label="Start",
                size=3)
h3

#Create density ridgeline plot showing DOY at the end of spring migration
a4<-ggplot(data=mig.summer, aes(x = DFP_End, y = timing, colour=timing, fill=timing,group=timing)) +
  geom_density_ridges(alpha = 0.60, scale = 500 ,size=0.5) + theme_ridges() +
  geom_vline(xintercept = 0, size=1)+
  scale_fill_manual(values = c("#9966CC","#339966","#FF9933"),
                    name = "Depature from winter range",
                    labels=c("Early","Mid","Late")) +
  scale_color_cyclical(values = c("#9966CC","#339966","#FF9933"))+
  xlim(-100,90)+
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5)))+
  ggtitle("b2")

a4

d4 <- a4+xlab("Days from peak IRG") + ylab("Density")
d4

e4<- d4 + theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.line.x = element_line(colour = "black"),
                             axis.line.y = element_line(colour = "black"))
e4
f4<- e4+ theme(axis.text.x=element_text(size=9, colour="black"),
               axis.text.y=element_blank(),
               axis.ticks.length.x = unit(.1, "cm"),
               axis.ticks.y = element_blank(),
               axis.title.x=element_text(size=9,margin=margin(5,5,5,5)),
               axis.title.y=element_blank(),
               plot.title = element_text(size=9,face="bold",vjust = - 7,hjust=-0.09))

f4
g4<- f4 + theme (plot.background=element_blank(),
                   plot.margin = unit(c(.1,.1,.1,.1),"cm"),
                   legend.background = element_rect(),
                   legend.title=element_blank(),
                   legend.text=element_text(size=9),
                   legend.position="bottom",
                   legend.margin=margin(0,0,0,0),
                   legend.box.margin=margin(-10,-10,-1,-10),
                   axis.line.x = element_line(colour = "black"),
                   axis.line.y = element_line(colour = "black"))

g4

h4<-g4+annotate('text', x= -95,y=490,
                label="End",
                size=3)
h4

#Get legend from last plot
legend<-get_legend(h4)

#Remove legend from last plot
h4 <- h4 + theme (plot.background=element_blank(),
                      legend.position = "none",
                      legend.title = element_blank(),
                      legend.text = element_text(size=12),
                      plot.margin = unit(c(.1,.1,.1,.1),"cm"),
                      axis.line.x = element_line(colour = "black", size = 0.5),
                      axis.line.y = element_line(colour = "black", size = 0.5))


h4

#Plot multipanel
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Figures")
jpeg("Ortega_et_al_2021_FigureS2.jpeg", width = 5, height = 5, units = 'in', res = 350)

l<-grid.arrange(h1,h3,h2,h4,legend, ncol=2, nrow = 3, 
             layout_matrix = rbind(c(1,2),c(3,4),c(5,5)),
             widths = c(2.73,2.47), heights = c(2.10,2.5,0.25))
dev.off()

#Save as PDF for InfoGraphics Lab
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/DataForInfoGraphics/ExtendedDataFigures")
ggsave("Ortega_et_al_2021_FigureS2.PDF", l,width = 5, height = 5, units = 'in')

#Clean working environment
rm(list=ls()[! ls() %in% c("all")])

#### CONDUCT PAIRED T-TEST TO TEST FOR BEHAVIORAL COMPENSATION ####

#Subset early, mid, and late migrants
early<-subset(all,all$timing=="early")
mid<-subset(all,all$timing=="mid")
late<-subset(all,all$timing=="late")

#Test for normal distribution
d1 <- with(early, 
          DFP[Cat == "Migration start"] - DFP[Cat == "Migration end"])
d2 <- with(mid, 
          DFP[Cat == "Migration start"] - DFP[Cat == "Migration end"])
d3 <- with(late, 
          DFP[Cat == "Migration start"] - DFP[Cat == "Migration end"])

#Shapiro-Wilk normality test for the differences
shapiro.test(d1)
shapiro.test(d2)
shapiro.test(d3) #Only late migrants violate homogeneity

#Compute t-test (wilcoxon for late migrants)
res.early <- t.test(DFP ~ Cat, data = early, paired = TRUE, alternative = "two.sided")
res.early

res.mid <- t.test(DFP ~ Cat, data = mid, paired = TRUE, alternative = "two.sided")
res.mid

res.late <- wilcox.test(DFP ~ Cat, data = late, paired = TRUE, alternative = "two.sided")
res.late

#To obtain effect sizes:
Zstat<-qnorm(res.late$p.value/2)
Zstat

tibble.d<-tibble(late)

stat.test <- tibble.d  %>%
  rstatix::wilcox_test(DFP ~ Cat, paired = TRUE) %>%
  add_significance()
stat.test

tibble.d  %>%
  wilcox_effsize(DFP ~ Cat, paired = TRUE)

#### PLOT DFP AT START AND END OF SPRING MIGRATION IN PAIRED DESIGN ####
#Early migrants
a.early<-ggplot(data=early,aes(x=Level,y=DFP))+
  geom_boxplot(aes(colour=Level),outlier.shape = NA)+
  geom_line(aes(group=id_yr),colour="#000000",alpha=0.3)+
  geom_point(aes(colour=Level),size=1, shape=16)+
  geom_hline(yintercept = 0, linetype="dashed", size=0.75)+
  ylim(-70,30)+
  ggtitle("a")
a.early

b.early<- a.early + scale_color_manual(values=c("#9966CC","#9966CC"))
b.early

c.early<- b.early + theme_bw() + theme(panel.border = element_blank(),
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(),
                                             axis.line.x = element_line(colour = "black", size = 0.5),
                                             axis.line.y = element_line(colour = "black", size = 0.5))
c.early

d.early<- c.early+ theme(axis.text.x=element_text(size=9,colour="black"),
                         axis.ticks.length.x = unit(.1, "cm"),
                         axis.text.y=element_text(size=9,colour="black"),
                         axis.title.x=element_blank(),
                         axis.title.y=element_text(size=9,margin=margin(5,5,5,5)),
                         plot.title = element_text(size=9,face="bold",vjust = 0,hjust=-.06))
d.early

e.early<- d.early + theme (plot.background=element_blank(),
                                 plot.margin = unit(c(.1,.1,.05,.1),"cm"),
                                 legend.position = "none",
                                 axis.line.x = element_line(colour = "black", size = 0.5),
                                 axis.line.y = element_line(colour = "black", size = 0.5))

e.early

f.early<-e.early+xlab(" ") + ylab("Days from peak IRG")

f.early

#Mid-migrants
range(mid$DFP)
a.mid<-ggplot(data=mid,aes(x=Level,y=DFP))+
  geom_boxplot(aes(colour=Level),outlier.shape = NA)+
  geom_line(aes(group=id_yr),colour="#000000",alpha=0.3)+
  geom_point(aes(colour=Level),size=1, shape=16)+
  geom_hline(yintercept = 0, linetype="dashed", size=0.75)+
  ylim(-42,43)+
  ggtitle("b")
a.mid

b.mid<- a.mid + scale_color_manual(values=c("#339966","#339966"))
b.mid

c.mid<- b.mid + theme_bw() + theme(panel.border = element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       axis.line.x = element_line(colour = "black", size = 0.5),
                                       axis.line.y = element_line(colour = "black", size = 0.5))
c.mid

d.mid<- c.mid+ theme(axis.text.x=element_text(size=9,colour="black"),
                     axis.ticks.length.x = unit(.1, "cm"),
                     axis.text.y=element_text(size=9,colour="black"),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     plot.title = element_text(size=9,face="bold",vjust = 0,hjust=-.06))
d.mid

e.mid<- d.mid + theme (plot.background=element_blank(),
                           plot.margin = unit(c(.1,.1,.05,.1),"cm"),
                           legend.position = "none",
                           axis.line.x = element_line(colour = "black", size = 0.5),
                           axis.line.y = element_line(colour = "black", size = 0.5))

e.mid

f.mid<-e.mid+xlab(" ") + ylab(" ")

f.mid

#Late migrants
range(late$DFP)
a.late<-ggplot(data=late,aes(x=Level,y=DFP))+
  geom_boxplot(aes(colour=Level),outlier.shape = NA)+
  geom_line(aes(group=id_yr),colour="#000000",alpha=0.3)+
  geom_point(aes(colour=Level),size=1, shape=16)+
  geom_hline(yintercept = 0, linetype="dashed", size=0.75)+
  ylim(-26,73)+
  ggtitle("c")

a.late

b.late<- a.late + scale_color_manual(values=c("#FF9933","#FF9933"))
b.late

c.late<- b.late + theme_bw() + theme(panel.border = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.line.x = element_line(colour = "black", size = 0.5),
                                   axis.line.y = element_line(colour = "black", size = 0.5))
c.late

d.late<- c.late+ theme(axis.text.x=element_text(size=9,colour="black"),
                       axis.ticks.length.x = unit(.1, "cm"),
                       axis.text.y=element_text(size=9,colour="black"),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       plot.title = element_text(size=9,face="bold",vjust = 0,hjust=-.06))
d.late

e.late<- d.late + theme (plot.background=element_blank(),
                       plot.margin = unit(c(.1,.1,.05,.1),"cm"),
                       legend.position = "none",
                       axis.line.x = element_line(colour = "black", size = 0.5),
                       axis.line.y = element_line(colour = "black", size = 0.5))

e.late

f.late<-e.late+xlab(" ") + ylab(" ")

f.late

#Make multipanel
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Figures")
jpeg("Ortega_et_al_2021_FigureS2.jpeg", width = 6, height = 2.25, units = 'in', res = 350)

#Create blank plot with desired format
l<-grid.arrange(f.early,f.mid,f.late, ncol=3, nrow = 1, 
             layout_matrix = rbind(c(1,2,3)),
             bottom = textGrob("Position along migratory route",hjust=0.25,vjust = 0.5,gp=gpar(fontsize=9)),
             widths = c(2.27,2,2), heights = c(2.25))
dev.off()

setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/DataForInfoGraphics/ExtendedDataFigures")
ggsave("Ortega_et_al_2021_FigureS2.PDF", l,width = 6, height = 2.25, units = 'in')
