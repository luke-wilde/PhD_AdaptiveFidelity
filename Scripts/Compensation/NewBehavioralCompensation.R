#Code for identifying partial compensators, full compensators, surfers, and non-compensators
#Written by Anna Ortega, Wyoming Cooperative Fish and Wildlife Research Unit
#May 5, 2022

#Clean working environment
rm(list=ls())

#Import data with degree of mismatch at start and end of spring migration
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Metadata")
mis.start<-read.csv("MismatchAtStartofSpringMigration.csv",header=TRUE,sep=",")
head(mis.start)
mis.start<-mis.start[c(2,5)]
colnames(mis.start)<-c("id_yr","DFP_start")
mis.end<-read.csv("MismatchAtEndofSpringMigration.csv",header=TRUE,sep=",")
mis.end<-mis.end[c(2,5)]
colnames(mis.end)<-c("id_yr","DFP_end")
head(mis.end)

range(mis.start$DFP_start)
d<-merge(mis.start,mis.end,by=c("id_yr"))
head(d)

d$abs.DFP.start<-abs(d$DFP_start)
d$abs.DFP.end<-abs(d$DFP_end)
d$abs.diff.DFP<-d$abs.DFP.end-d$abs.DFP.start
head(d)

#Classify individuals as full compensators, partial compensators, perfect surfers or non-compensators
d$Comp<-NA

#For individuals that left behind peak IRG
d$Comp<-ifelse(d$DFP_start > 7 & d$abs.diff.DFP>0,"Behind non-compensator",d$Comp)
d$Comp<-ifelse(d$DFP_start > 7 & d$abs.diff.DFP<0, "Behind full compensator",d$Comp)

#For individuals that left during peak IRG
d$Comp<-ifelse(d$DFP_start >= -7 & d$DFP_start <=7& d$abs.diff.DFP>0,"During non-compensator",d$Comp)
d$Comp<-ifelse(d$DFP_start  >= -7 & d$DFP_start <=7 & d$DFP_end>=-7 & d$DFP_end<=7,"Perfect surfer",d$Comp)

#For individuals that left ahead peak IRG
d$Comp<-ifelse(d$DFP_start < -7 & d$abs.diff.DFP>0,"Ahead non-compensator",d$Comp) 
d$Comp<-ifelse(d$DFP_start < -7 & d$abs.diff.DFP<0, "Ahead full compensator",d$Comp)

#Further refine classification of full compensators and non-compensators into partial compensators
b.n<-subset(d,d$Comp=="Behind non-compensator")
b.n$Comp<-ifelse(b.n$abs.diff.DFP<=7,"Behind partial compensator",b.n$Comp)
table(b.n$Comp)

b.f<-subset(d,d$Comp=="Behind full compensator")
b.f$Comp<-ifelse(b.f$abs.diff.DFP>=-7,"Behind partial compensator",b.f$Comp)
table(b.f$Comp)

a.n<-subset(d,d$Comp=="Ahead non-compensator")
a.n$Comp<-ifelse(a.n$abs.diff.DFP<=7,"Ahead partial compensator",a.n$Comp)
table(a.n$Comp)

a.f<-subset(d,d$Comp=="Ahead full compensator")
a.f$Comp<-ifelse(a.f$abs.diff.DFP>=-7,"Ahead partial compensator",a.f$Comp)
table(a.f$Comp)

d<-d[d$Comp!="Behind non-compensator",]
d<-d[d$Comp!="Behind full compensator",]
d<-d[d$Comp!="Ahead non-compensator",]
d<-d[d$Comp!="Ahead full compensator",]

d<-rbind(d,b.n,b.f,a.n,a.f) #Should have 152 unique animal-years
table(d$Comp)
head(d)

d$compensation<-NA
d$compensation<-ifelse(d$Comp=="Ahead full compensator","Full compensator",d$compensation)
d$compensation<-ifelse(d$Comp=="Ahead non-compensator","Non-compensator",d$compensation)
d$compensation<-ifelse(d$Comp=="Ahead partial compensator","Partial compensator",d$compensation)
d$compensation<-ifelse(d$Comp=="Behind full compensator","Full compensator",d$compensation)
d$compensation<-ifelse(d$Comp=="Behind non-compensator","Non-compensator",d$compensation)
d$compensation<-ifelse(d$Comp=="Behind partial compensator","Partial compensator",d$compensation)
d$compensation<-ifelse(d$Comp=="During non-compensator","Non-compensator",d$compensation)
d$compensation<-ifelse(d$Comp=="Perfect surfer","Perfect surfer",d$compensation)

table(d$compensation)
table(d$Comp)
head(d)

#Double check that all compensators became closer to peak IRG
#at the end of spring migration compared with the start of spring migration
full.ahead<-subset(d,d$Comp=="Ahead full compensator")
range(full.ahead$abs.diff.DFP) #Should be negative

full.behind<-subset(d,d$Comp=="Behind full compensator")
range(full.behind$abs.diff.DFP) #Should be negative

non.ahead<-subset(d,d$Comp=="Ahead non-compensator")
range(non.ahead$abs.diff.DFP) #Should be positive

non.behind<-subset(d,d$Comp=="Behind non-compensator")
range(non.behind$abs.diff.DFP) #Should be positive

#Classify individuals as those that left ahead, behind, or during peak IRG
d$Level<-NA
d$Level<-ifelse(d$Comp=="Ahead full compensator","Ahead",d$Level)
d$Level<-ifelse(d$Comp=="Ahead partial compensator","Ahead",d$Level)
d$Level<-ifelse(d$Comp=="Behind partial compensator","Behind",d$Level)
d$Level<-ifelse(d$Comp=="Behind full compensator","Behind",d$Level)
d$Level<-ifelse(d$Comp=="Behind non-compensator","Behind",d$Level)
d$Level<-ifelse(d$Comp=="Perfect surfer","During",d$Level)
d$Level<-ifelse(d$Comp=="During non-compensator","During",d$Level)
d$Level<-ifelse(d$Comp=="Ahead non-compensator","Ahead",d$Level)

table(d$Level)

#Classify individuals as those that partially compensated, fully compensated, did not compensate, or surfed
d$Cat<-NA
d$Cat<-ifelse(d$Comp=="Ahead full compensator","Full compensator",d$Cat)
d$Cat<-ifelse(d$Comp=="Ahead partial compensator","Partial compensator",d$Cat)
d$Cat<-ifelse(d$Comp=="Behind partial compensator","Partial compensator",d$Cat)
d$Cat<-ifelse(d$Comp=="Behind full compensator","Full compensator",d$Cat)
d$Cat<-ifelse(d$Comp=="Behind non-compensator","Non-compensator",d$Cat)
d$Cat<-ifelse(d$Comp=="Perfect surfer","Perfect surfer",d$Cat)
d$Cat<-ifelse(d$Comp=="During non-compensator","Non-compensator",d$Cat)
d$Cat<-ifelse(d$Comp=="Ahead non-compensator","Non-compensator",d$Cat)
table(d$Cat)

table(is.na(d$Comp))
table(is.na(d$Cat))
table(is.na(d$Level))

table(d$Level)
table(d$Comp)

#Export data
head(d)
id<-d[c("id_yr","Cat","Level")]
colnames(id)<-c("id_yr","compensation","departure.rel.greenup")
head(id)
table(id$compensation)
table(id$departure.rel.greenup)

setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Metadata")
write.csv(id,"IDList_BehavioralCompensation_20FEB2023.csv")

library(Rmisc)
sum.1=summarySE(d,measurevar = c("abs.DFP.start"),groupvars = c("Cat"))
sum.1$true.ci<-sum.1$se*1.96
sum.1

sum.1$abs.DFP.start<-round(sum.1$abs.DFP.start,2)
sum.1$true.ci<-round(sum.1$true.ci,2)
sum.1<-sum.1[c("Cat","abs.DFP.start","true.ci")]
sum.1

sum.2=summarySE(d,measurevar = c("abs.DFP.end"),groupvars = c("Cat"))
sum.2$true.ci<-sum.2$se*1.96
sum.2

sum.2$abs.DFP.end<-round(sum.2$abs.DFP.end,2)
sum.2$true.ci<-round(sum.2$true.ci,2)
sum.2<-sum.2[c("Cat","abs.DFP.end","true.ci")]
sum.2

# #Export data for Source Data File in Nature Communications
# names(d)
# d<-d[c("id_yr","abs.DFP.start","abs.DFP.end","compensation","Level")]
# colnames(d)<-c("id_yr","abs.DFP.start","abs.DFP.end","compensation","departure.greenup")
# d<-d[order(d$id_yr),]
# table(d$compensation)
# table(d$departure.greenup)
# 
# setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Submission/NatureCommunicationsTransfer/SourceData")
# write.csv(d,"Compesation_AbsolutedDFP.csv")

#Conduct ordinal logistic regression
#to evaluate effect of mismatch on likelihood of being a full compensator
d$Cat<-factor(d$Cat, levels=c("Non-compensator","Perfect surfer","Partial compensator","Full compensator"))

library(MASS)
d$Cat<-as.factor(d$Cat)
model= polr(Cat ~ abs.DFP.start, data = d, Hess = TRUE)
summary(model)
(ctable <- coef(summary(model)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = p))
(ci <- confint(model))
exp(coef(model))
exp(cbind(OR = coef(model), ci))

# plot(Effect(focal.predictors = "abs.start",model),
#      ylab="Probability",xlab="Absolute Days-From-Peak at start of spring migration",
#      main=" ",cex.lab=9)

#Plot predicted probabilies
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Figures")
jpeg("Ortega_et_al_2021_FigureS7.jpeg", width = 6, height = 2, units = 'in', res = 350)

library(ggeffects)
library(ggplot2)
mydf<-ggpredict(model, terms = "abs.DFP.start[all]")
mydf$response.level
head(mydf)

data_text <- data.frame(
  response.level   = c("Full compensator","Partial compensator","Perfect surfer","Non-compensator"),
  label = c("a","b","c","d"))
data_text$response.level<-factor(data_text$response.level, levels=c("Full compensator","Partial compensator","Perfect surfer","Non-compensator"))

mydf$response.level<-factor(mydf$response.level, levels=c("Full compensator","Partial compensator","Perfect surfer","Non-compensator"))
library(egg)
a1<-ggplot(mydf, aes(x = x, y = predicted, colour = response.level)) +
  facet_grid(. ~ response.level,scales="free")+
  geom_line(colour="#000000")+
  geom_ribbon(data = mydf, aes(x=x,y = NULL, ymin = conf.low, ymax = conf.high,
                                     color = NULL, fill = NULL),alpha=0.25)+
  geom_text(data=data_text, aes(x=3, y=1, label=label), fontface='bold', size=3,colour="#000000")

a1
b1 <- a1 + labs(y = "Probability",
                x = "Absolute days from peak IRG at start of spring migration")
b1

c1<-b1 +
  theme_bw()+
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour="black", fill="#FFFFFF"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
c1

d1<- c1+ theme(axis.text.x=element_text(size=9,colour="black"),
               axis.text.y=element_text(size=9,colour="black"),
               axis.ticks.length=unit(.1, "cm"),
               axis.title.x=element_text(size=9,margin=margin(5,5,5,5)),
               axis.title.y=element_text(size=9,margin=margin(5,5,5,5)),
               plot.title = element_text(size=9,face="bold"))

d1
e1<- d1 + theme (plot.margin = unit(c(.1,.1,.1,.1),"cm"),
                 legend.background = element_rect(),
                 legend.title=element_blank(),
                 legend.text=element_text(size=9),
                 legend.position="bottom",
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(-10,-10,-1,-10),
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"))
e1

dev.off()

setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/DataForInfoGraphics/ExtendedDataFigures")
ggsave("Ortega_et_al_2021_FigureS7.PDF", e1,width = 6, height = 2, units = 'in')

#Clean working environment
rm(a.f,a.n,a1,b.f,b.n,b1,c1,ctable,d1,e1,full.ahead,full.behind,mis.end,mis.start,
   model,mydf,non.ahead,non.behind,sum.1,sum.2,ci,p)

#Compare cumulative IRG among compensation categories ####
#Import cumulative IRG (exposure to spring green-up among individuals)
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data")
f<-read.csv("CumIRG.csv",header=TRUE,sep=",")
head(f)

all<-merge(d,f,by=c("id_yr"))
all$Cat<-factor(all$Cat, levels=c("Full compensator","Partial compensator","Perfect surfer","Non-compensator"))
all$Level<-factor(all$Level, levels=c("Ahead","During","Behind"))
head(all)

# library(Rmisc)
sum.1 = summarySE(all,measurevar="cum.IRG.day", groupvars=c("Level","Cat"))
sum.1$true.ci<-sum.1$se*1.96
sum.1

m2<-aov(cum.IRG.day~Comp,data=all)
summary(m2)
 
tukey.test <- TukeyHSD(m2)
tukey.test

# library(ggplot2)
dodge <- position_dodge(width=0.5)  
# 
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Figures")
jpeg("Ortega_et_al_2021_FFFFFigureS8.jpeg", width = 6, height = 2.75, units = 'in', res = 350)
library(ggplot2)

a1<-ggplot(sum.1, aes(x=Level,y=cum.IRG.day))+
  geom_errorbar(aes(ymin=cum.IRG.day-true.ci,
                    ymax=cum.IRG.day+true.ci),
                width=0, size=0.5,colour="#999999",position=dodge)+
  geom_point(position=dodge,size=1)+
  ylim(0,4.1)

a1
#


f1 <- a1 + labs(y = "Exposure to spring green-up",
                x = "Start of spring migration relative to green wave")
f1

g1<-f1 +  facet_grid(. ~ Cat,scales="free") +
  theme_bw()+
  theme(strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9),
        strip.background = element_rect(colour="black", fill="#FFFFFF"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g1

h1<- g1+ theme(axis.text.x=element_text(size=9,colour="black",angle=45,hjust=1),
               axis.text.y=element_text(size=9,colour="black"),
               axis.ticks.length=unit(.1, "cm"),
               axis.title.x=element_text(size=9,margin=margin(15,15,15,15)),
               axis.title.y=element_text(size=9,margin=margin(15,5,15,1)),
               plot.title = element_text(size=9,face="bold"))

h1
i1<- h1 + theme (plot.margin = unit(c(.5,.5,.1,.5),"cm"),
                 legend.background = element_rect(),
                 legend.title=element_blank(),
                 legend.text=element_text(size=9),
                 legend.position="bottom",
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(-10,-10,-1,-10),
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"))
i1

dev.off()

#Clean working environment
rm(a1,all,dodge,f,f1,g1,h1,i1,m2,sum.1,tukey.test)

#Determine effect of age on probability of deer becoming
#full compensators versus partial compensators, perfect surfers, or non-compensators
#Import age data
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data/Age")
age<-read.csv("RD2H_AgeAtStartOfSpringMigration.csv", header=TRUE, sep=",")
head(age)
str(age)

#Floor age values
age$AgeF<-floor(age$Age)
table(is.na(age$AgeF)) #Should be all FALSE

#Merge age data with dataframe containing degree of compensation
w<-merge(d,age,by=c("id_yr"))
head(w)
u<-unique(w$id_yr) #51 unique individuals with age data
table(w$Age)
range(w$AgeF)

#Conduct multiple linear regression to see effect of age on odds of compensation
head(w)
range(w$AgeF)
head(w)

#Plot histogram of ages
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Figures")
jpeg("Ortega_et_al_2021_AgeDistribution_ForReviewers.jpeg", width = 6, height = 6, units = 'in', res = 350)

a1<-ggplot(w, aes(x=AgeF)) + 
  geom_histogram(binwidth=1,color="darkblue",fill="lightblue")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13))
a1

d1 <- a1+xlab("Age") + ylab("Count")

e1<- d1 + theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.line.x = element_line(colour = "black"),
                             axis.line.y = element_line(colour = "black"))
e1
f1<- e1+ theme(axis.text.x=element_text(size=9,colour="black"),
               axis.text.y=element_text(size=9,colour="black"),
               axis.ticks.length.x = unit(.1, "cm"),
               axis.ticks.length.y = unit(.1, "cm"),
               axis.title.x=element_text(size=9,margin=margin(5,5,5,5)),
               axis.title.y=element_text(size=9,margin=margin(5,5,5,5)),
               plot.title = element_text(size=9,face="bold",vjust = - 7,hjust=-0.09))
f1
g1<- f1 + theme (plot.background=element_blank(),
                 plot.margin = unit(c(.1,.1,.1,.1),"cm"),
                 legend.position="none",
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"))

g1

dev.off()

w$Cat<-factor(w$Cat, levels=c("Non-compensator","Perfect surfer","Partial compensator","Full compensator"))

library(MASS)
w$Cat<-as.factor(w$Cat)
str(w)

#Just to see the alternative and make sure model is working
#w$AgeF<-ifelse(w$Cat=="Full compensator",sample(7:9,1),w$AgeF)
#head(w)

model= polr(Cat ~ AgeF , data = w, Hess = TRUE)
summary(model)
(ctable <- coef(summary(model)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = p))
(ci <- confint(model))
exp(coef(model))
exp(cbind(OR = coef(model), ci))

library(effects)
Effect(focal.predictors = "AgeF",model)
plot(Effect(focal.predictors = "AgeF",model))

#No need to plot figure

#Clean working environment
rm(age,ctable,model,w,ci,p,u)

#Now, compare movement rate and stopover use among categories of compensation ####
#Import movement rate and stopover use
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data/MovementRateStopoverUse")
rate<-read.csv("MovementRate.csv")
stop<-read.csv("StopoverUse.csv")

#Merge compensation dataframe with movement rate
d.rate<-merge(d,rate,by=c("id_yr"))
d.stop<-merge(d,stop,by=c("id_yr"))

d.rate$Cat<-factor(d.rate$Cat, levels=c("Full compensator","Partial compensator","Perfect surfer","Non-compensator"))
d.rate$Level<-factor(d.rate$Level, levels=c("Ahead","During","Behind"))

d.stop$Cat<-factor(d.stop$Cat, levels=c("Full compensator","Partial compensator","Perfect surfer","Non-compensator"))
d.stop$Level<-factor(d.stop$Level, levels=c("Ahead","During","Behind"))

library(Rmisc)
sum.1 = summarySE(d.rate,measurevar="rate.km.day", groupvars=c("Cat","Level"))
sum.1$true.ci<-sum.1$se*1.96
sum.1

sum.2 = summarySE(d.stop,measurevar="stopover.day", groupvars=c("Cat","Level"))
sum.2$true.ci<-sum.2$se*1.96
sum.2

table(d.rate$compensation)
sub<-subset(d.stop,d.stop$compensation=="Non-compensator")
m1<-mean(sub$stopover.day)
m2<-median(sub$stopover.day)

ggplot(data=sub, aes(x = stopover.day, y = 1)) +
  geom_density_ridges(alpha = 0.60, scale = 500 ,size=0.5) + theme_ridges()+
  geom_vline(xintercept=m2)

m1<-aov(rate.km.day~Comp,data=d.rate)
summary(m1)
summary(m1)[[1]][1,5]

m2<-aov(stopover.day~Comp,data=d.stop)
summary(m2)

tukey.test <- TukeyHSD(m1)
tukey.test
print(tukey.test[[1]][3,4],digits=20)

all<-merge(sum.1,sum.2,by=c("Cat","Level"))
head(all)

names(all)
all<-all[c(1:2,4,6,10,12)]
colnames(all)<-c("Cat","Level","rate","se.rate","stopover.use","se.stop")
head(all)

range(all$rate)
range(all$stopover.use)

#Create graph showing movement rate and stopover use among compensation category
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Figures")
jpeg("Ortega_et_al_2021_FigureS6.jpeg", width = 6, height = 2.50, units = 'in', res = 350)
library(ggplot2)

head(all)
range(all$rate) #2-8
range(all$stopover.use) #6-40

#Transpose stopover use based on range of movement rate and stopover use
#Need to create graph with two y-axes of different scales
sum.2$variable<-sum.2$stopover.day*(8/40)
sum.2$CI<-(sum.2$se*1.96)*(8/40)
sum.2$type<-"Stopover use"
names(sum.2)

sum.2<-sum.2[c(1:2,9:11)]
head(sum.2)

sum.1$variable<-sum.1$rate
sum.1$CI<-sum.1$se*1.96
sum.1$type<-"Movement rate"

sum.1<-sum.1[c(1:2,9:11)]
head(sum.1)

all<-rbind(sum.1,sum.2)
head(all)

dodge <- position_dodge(width=0.5)  

data_text <- data.frame(
  Cat   = c("Full compensator","Partial compensator","Perfect surfer","Non-compensator"),
  type   = c("Movement rate","Stopover use"),
  label = c("a","b","c","d"))
data_text$Cat<-factor(data_text$Cat, levels=c("Full compensator","Partial compensator","Perfect surfer","Non-compensator"))

a1<-ggplot(all, aes(x=Level,y=variable,colour=type,group=type))+ 
    geom_errorbar(aes(ymin=variable-CI, 
                      ymax=variable+CI), 
                      size=0.5,width=0,colour="#999999",position=dodge)+
  geom_point(position=dodge,size=1)

a1

# b1<- a1 + geom_errorbar(aes(ymin=stopover.use*7/34-(se.stop*1.96)*7/34, 
#                             ymax=stopover.use*7/34+(se.stop*1.96)*7/34), 
#                             width=.1, size=1,colour="#999999", position=dodge)+
#           geom_point(aes(y = stopover.use*7/34,colour="Stopover use"),size=8, position=dodge)
# b1

d1 <- a1 + scale_y_continuous(sec.axis = sec_axis(~.*(40/8), name = "Days allocated to stopovers"))
d1

e1<-d1 + scale_colour_manual(values=c("#0000FF","#CC0000"))

e1

f1 <- e1 + labs(y = "Rate of movement (km/day)",
              x = "Departure relative to green wave",
              colour = "Movement parameter")
f1

# # New facet label names for supp variable
# graph.labs <- c("A", "B", "C")
# names(graph.labs) <- c("Full compensator", "Drifter","Surfer")
#g1<-f1 +  facet_grid(. ~ Cat,scales="free",labeller = labeller(Cat=graph.labs)) +
  
g1<-f1 +  facet_grid(. ~ Cat,scales="free") +
  theme_bw()+
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour="black", fill="#FFFFFF"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(data=data_text, aes(x=c(0.6,0.65,0.5,0.7), y=10, label=label), fontface='bold', size=3,colour="#000000")

g1


h1<- g1+ theme(axis.text.x=element_text(size=9,colour="black",angle=45,hjust=1),
               axis.text.y=element_text(size=9,colour="black"),
               axis.ticks.length=unit(.1, "cm"),
               axis.title.x=element_text(size=9,margin=margin(5,5,5,5)),
               axis.title.y=element_text(size=9,margin=margin(5,5,5,5),colour="#0000FF"),
               axis.title.y.right = element_text(vjust=2,colour="#CC0000"),
               plot.title = element_text(size=9,face="bold"))

h1
i1<- h1 + theme (plot.margin = unit(c(.1,.1,.1,.1),"cm"),
                 legend.background = element_rect(),
                 legend.title=element_blank(),
                 legend.text=element_text(size=9),
                 legend.position="bottom",
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(-10,-10,-1,-10),
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"))
i1

dev.off()

setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/DataForInfoGraphics/ExtendedDataFigures")
ggsave("Ortega_et_al_2021_FigureS6.PDF", i1,width = 6, height = 2.50, units = 'in')

#Now, export movement rate and stopover use for full and partial compensators in 2017
r.17<-subset(d.rate,d.rate$year=="2017")
s.17<-subset(d.stop,d.stop$year=="2017")

all.17<-merge(r.17,s.17,by=c("id_yr"),all=TRUE)
table(is.na(all.17$stopover.day)) #Do not have stopover use data for 10 individuals

#Reduce database to essential attributes only
names(all.17)
all.17<-all.17[c(1,7:9,14,27)]
colnames(all.17)<-c("id_yr","Comp","Level","Cat","rate.km.day","stopover.use")
head(all.17)

#Subset only full or partial compensators
all.17<-all.17[all.17$Cat %in% c("Full compensator"),]
head(all.17)
table(all.17$Cat)

#Remove unnecessary columns
names(all.17)
all.17<-all.17[c(1,5:6)]
head(all.17)

#Export movement rate and stopover use for full and partial compensatiors in 2017
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/VisualizingTheGreenWave")
write.csv(all.17,"AllCompensators_MovementRate_StopoverUse_2017.csv")

