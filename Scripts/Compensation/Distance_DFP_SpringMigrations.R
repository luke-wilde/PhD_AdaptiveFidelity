#Code for analyzing distance vs. days from peak green-up
#Written by Anna Ortega, Wyoming Cooperative Fish and Wildlife Research Unit
#March 5, 2023

#Load required libraries
library(dplyr)
library(ggplot)
library(raster)
library(stringr)

#Clean working environment
rm(list=ls())

#### IMPORT CLEAN GPS COLLAR DATA ####
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data/SpringMigrations")
load("AllLongDistance_SpringMigrations_2011_thru_2020_CleanForAnalyses.RData")

#Remove additional individuals from dataset for listed reasons
final.gps.data<-final.gps.data[final.gps.data$id_yr!="30_2013",] #Only one long-distance migrant in 2013
final.gps.data<-final.gps.data[final.gps.data$id_yr!="127_2015",] #Only one long-distance migrant in 2015

final.gps.data<-final.gps.data[final.gps.data$id_yr!="303_2018",] #Returned to winter range on 3/17/2018
final.gps.data<-final.gps.data[final.gps.data$id_yr!="303_2019",] #Returned to winter range on 2/22/2019

final.gps.data<-final.gps.data[final.gps.data$id_yr!="225_2016",] #First captured on Table Mtn stopover without prior GPS data
final.gps.data<-final.gps.data[final.gps.data$id_yr!="257_2016",] #First captured on Table Mtn stopover without prior GPS data
final.gps.data<-final.gps.data[final.gps.data$id_yr!="258_2016",] #First captured on Table Mtn stopover without prior GPS data
final.gps.data<-final.gps.data[final.gps.data$id_yr!="259_2016",] #First captured on Table Mtn stopover without prior GPS data
final.gps.data<-final.gps.data[final.gps.data$id_yr!="317_2017",] #First captured on Table Mtn stopover without prior GPS data
final.gps.data<-final.gps.data[final.gps.data$id_yr!="319_2017",] #First captured on Table Mtn stopover without prior GPS data
final.gps.data<-final.gps.data[final.gps.data$id_yr!="320_2017",] #First captured on Table Mtn stopover without prior GPS data
final.gps.data<-final.gps.data[final.gps.data$id_yr!="321_2017",] #First captured on Table Mtn stopover without prior GPS data

final.gps.data<-final.gps.data[final.gps.data$id_yr!="332_2018",] #Recaptured during migration
final.gps.data<-final.gps.data[final.gps.data$id_yr!="254_2019",] #Recaptured during migration
final.gps.data<-final.gps.data[final.gps.data$id_yr!="279_2019",] #Recaptured during migration

final.gps.data<-final.gps.data[final.gps.data$id_yr!="456_2020",] #Male 9-month-old that was captured in March 2020

id.list.1<-unique(final.gps.data$id)
id.list.1<-unique(final.gps.data$id_yr)

#Clean dataframe to export as csv file (for later use)
d<-final.gps.data[c(1:2,5:6)]
d<-d[c("id_yr","date","x","y")]
d<-d[order(d$id_yr,d$date),]
head(d)

#Write csv for final dataset
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Movebank")
write.csv(d,"RD2H_SpringMigrations_LongDistance_2011_thru_2020.csv")

# #Analyze sample sizes
# phase1<-subset(final.gps.data,final.gps.data$year<=2013)
# u1<-unique(phase1$id) #There should be 8 unique individuals
# u2<-unique(phase1$id_yr) #There should be 12 animal-years
# 
# phase2<-subset(final.gps.data,final.gps.data$year>=2014)
# u3<-unique(phase2$id) #There should be 64 unique individuals
# u4<-unique(phase2$id_yr) #There should be 140 animal-years

#Make data spatially aware
proj <- "+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs+towgs84=0,0,0"   #name a new proj
coordinates(final.gps.data) <- c("x","y")
proj4string(final.gps.data) <- CRS(proj)

#Rename dataframe
gps.data<-as.data.frame(final.gps.data)
rm(final.gps.data,d)

#Determine number of unique animal-years: 152 unique animal-years
id.list1<-unique(gps.data$id_yr)
id.list2<-unique(gps.data$id)

#Import migration timing data
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/Migrations")
Mig.Data<-read.csv("AllSpringFallMigrationTiming_2011_thru_2021.csv", header=TRUE, sep=",")
Mig.Data$id_yr<-paste(Mig.Data$AID, Mig.Data$Year, sep="_")
head(Mig.Data)

#Ensure that migration start and end dates are in POSIXct format
Mig.Data$SpringMS <- as.POSIXct(Mig.Data$SpringMS, format="%m/%d/%Y %H:%M")
Mig.Data$SpringME <- as.POSIXct(Mig.Data$SpringME, format="%m/%d/%Y %H:%M")
names(Mig.Data)

#Remove unnecessary columns and make sure there are no missing values for migration start and end dates
Mig.Data<-Mig.Data[c(6,9,19)]
Mig.Data<-Mig.Data[!is.na(Mig.Data[,1]),]
Mig.Data<-Mig.Data[!is.na(Mig.Data[,2]),]
head(Mig.Data)
names(Mig.Data)

#Change start and end dates to Julian day
library(lubridate)
Mig.Data$SpringMSJul<-yday(Mig.Data$SpringMS)
Mig.Data$SpringMEJul<-yday(Mig.Data$SpringME)

#Merge migration timing data with GPS collar data
d <- merge(gps.data, Mig.Data, by=c("id_yr"))

#Make sure that the new dataset has the same number of unique animal years as the GPS dataset
u1<-as.data.frame(unique(gps.data$id_yr))
u2<-as.data.frame(unique(d$id_yr))

#Order dataframe by id_yr and jul
d$jul<-yday(d$date)
d<-d[order(d$id_yr,d$date),]

#Determine mean departure date for each year
t<-d[duplicated(d$id_yr)==FALSE,]  #gives first row of dataframe for each month
range(t$SpringMSJul)

library(Rmisc)
sum = summarySE(t,measurevar="SpringMSJul",groupvars = c("year"))
sum$true.ci<-sum$se*1.96
sum

#Categorize all individuals into early, mid, and late migrants
#2011
long.11<-subset(d,d$year=="2011")
table(long.11$year)
long.11$timing<-"mid"
d.long<-long.11[duplicated(long.11$id_yr)==FALSE,]  #gives first row of dataframe for each month
summary(d.long$SpringMSJul)
long.11$timing<-ifelse(long.11$SpringMSJul<=91,"early",long.11$timing)
long.11$timing<-ifelse(long.11$SpringMSJul>=92,"late",long.11$timing)
table(long.11$timing)
long.11$c.SpringMSJul<-long.11$SpringMSJul-median(long.11$SpringMSJul)
range(long.11$c.SpringMSJul)

#2012
long.12<-subset(d,d$year=="2012")
table(long.12$year)
long.12$timing<-"mid"
d.long<-long.12[duplicated(long.12$id_yr)==FALSE,]  #gives first row of dataframe for each month
summary(d.long$SpringMSJul)
long.12$timing<-ifelse(long.12$SpringMSJul<=78,"early",long.12$timing)
long.12$timing<-ifelse(long.12$SpringMSJul>=109,"late",long.12$timing)
table(long.12$timing)
long.12$c.SpringMSJul<-long.12$SpringMSJul-median(long.12$SpringMSJul)

#2014
long.14<-subset(d,d$year=="2014")
table(long.14$year)
long.14$timing<-"mid"
d.long<-long.14[duplicated(long.14$id_yr)==FALSE,]  #gives first row of dataframe for each month
summary(d.long$SpringMSJul)
long.14$timing<-ifelse(long.14$SpringMSJul<=112,"early",long.14$timing)
long.14$timing<-ifelse(long.14$SpringMSJul>=119,"late",long.14$timing)
table(long.14$timing)
long.14$c.SpringMSJul<-long.14$SpringMSJul-median(long.14$SpringMSJul)

#2016
long.16<-subset(d,d$year=="2016")
table(long.16$year)
long.16$timing<-"mid"
d.long<-long.16[duplicated(long.16$id_yr)==FALSE,]  #gives first row of dataframe for each month
summary(d.long$SpringMSJul)
long.16$timing<-ifelse(long.16$SpringMSJul<=96,"early",long.16$timing)
long.16$timing<-ifelse(long.16$SpringMSJul>=136,"late",long.16$timing)
table(long.16$timing)
long.16$c.SpringMSJul<-long.16$SpringMSJul-median(long.16$SpringMSJul)

#2017
long.17<-subset(d,d$year=="2017")
table(long.17$year)
long.17$timing<-"mid"
d.long<-long.17[duplicated(long.17$id_yr)==FALSE,]  #gives first row of dataframe for each month
summary(d.long$SpringMSJul)
long.17$timing<-ifelse(long.17$SpringMSJul<=86,"early",long.17$timing)
long.17$timing<-ifelse(long.17$SpringMSJul>=123,"late",long.17$timing)
table(long.17$timing)
long.17$c.SpringMSJul<-long.17$SpringMSJul-median(long.17$SpringMSJul)

#2018
long.18<-subset(d,d$year=="2018")
table(long.18$year)
long.18$timing<-"mid"
d.long<-long.18[duplicated(long.18$id_yr)==FALSE,]  #gives first row of dataframe for each month
summary(d.long$SpringMSJul)
long.18$timing<-ifelse(long.18$SpringMSJul<=82,"early",long.18$timing)
long.18$timing<-ifelse(long.18$SpringMSJul>=113,"late",long.18$timing)
table(long.18$timing)
long.18$c.SpringMSJul<-long.18$SpringMSJul-median(long.18$SpringMSJul)

#2019
long.19<-subset(d,d$year=="2019")
table(long.19$year)
long.19$timing<-"mid"
d.long<-long.19[duplicated(long.19$id_yr)==FALSE,]  #gives first row of dataframe for each month
summary(d.long$SpringMSJul)
long.19$timing<-ifelse(long.19$SpringMSJul<=78,"early",long.19$timing)
long.19$timing<-ifelse(long.19$SpringMSJul>=120,"late",long.19$timing)
table(long.19$timing)
long.19$c.SpringMSJul<-long.19$SpringMSJul-median(long.19$SpringMSJul)

#2020
long.20<-subset(d,d$year=="2020")
table(long.20$year)
long.20$timing<-"mid"
d.long<-long.20[duplicated(long.20$id_yr)==FALSE,]  #gives first row of dataframe for each month
summary(d.long$SpringMSJul)
long.20$timing<-ifelse(long.20$SpringMSJul<=82,"early",long.20$timing)
long.20$timing<-ifelse(long.20$SpringMSJul>=111,"late",long.20$timing)
table(long.20$timing)
long.20$c.SpringMSJul<-long.20$SpringMSJul-median(long.20$SpringMSJul)

#Merge all long-distance migrants
all.long<-rbind(long.11,long.12,long.14,long.16,long.17,long.18,long.19,long.20)
names(all.long)

# #Determine what individuals had more than two years of GPS data
# head(all.long)
# dup<-all.long[duplicated(all.long$id_yr)==FALSE,]
# names(dup)
# dup<-dup[c(1,7,15)]
# head(dup)
# 
# dup$temp<-1
# 
# h=data.frame(id=id.list2,
#              num.years=NA)
# head(h)
# 
# for(i in 1:length(id.list2)){
#   a<-subset(dup,dup$id==id.list2[[i]])
#   count<-sum(a$temp)
#   h[i,"num.years"]<-count
# }
# 
# k<-subset(h,h$num.years>=2)
# id.list3<-unique(k$id)
# 
# dup<-dup[(dup$id %in% k$id),]
# dup<-dup[order(dup$id_yr,dup$timing),]

# setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Metadata")
# write.csv(dup,"AnimalsWithMultipleYearsOfData.csv")

range(all.long$c.SpringMSJul)
table(all.long$SpringMSJul)
range(all.long$SpringMSJul) #Earliest migration on Feb 18 and latest migration on June 1

#Rename dataframe
f<-as.data.frame(all.long)

#Clean working environment and rename dataframe
rm(list=ls()[! ls() %in% c("f")])

#### EXTRACT GREEN-WAVE SURFING PARAMETERS ####
source("D:/R-code/Functions/ndviExtract.R")

metrics <- c("IRGVals","fittedVals","maxIRGdate") #IRGVals are fitted IRG and fittedVals are fitted NDVI
"temp" %in% names(f)   #verify this is false
for(w in 1:length(metrics)){
  f$temp <- ndviEXTRACT(XYdata=f, NDVImetric=metrics[w], NDVIfolder="D:/R-code/NDVI_Data_2020",
                        maxcpus=11, xname="x", yname="y", datesname="date", scaleIRG=TRUE,
                        xyCRS=CRS("+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs+towgs84=0,0,0"))
  names(f)[names(f)=="temp"] <- metrics[w]
}

range(f$IRGVals, na.rm=TRUE)
range(f$fittedVals, na.rm=TRUE)
range(f$maxIRGdate, na.rm=TRUE)

#Replace negative NDVI values with zero
f$fittedVals<-ifelse(f$fittedVals<0,0,f$fittedVals)
range(f$IRGVals, na.rm=TRUE) #Should be 0-1
range(f$fittedVals, na.rm=TRUE) #Should be 0-1
range(f$maxIRGdate, na.rm=TRUE)

#Rename dataframe
d<-f
rm(f)
d<-d[order(d$id_yr,d$date),]

head(d)

#### USE NSDS TO CALCULATE DISTANCE FROM WINTER RANGE FOR EACH LOCATION ####
#Create nsd object
nsd<-function(x, y){
  ns <- length(x)
  msd <- NULL
  for(i in 1:ns){
    msd <- c(msd, sqrt((x[i]-x[1])^2+(y[i]-y[1])^2) )
  }
  return(msd)
}

d$NSD<-rep(NA, nrow(d))
full.df<-d[0, ]

ids<-unique(d$id_yr)

#NSD analysis
for(i in 1:length(ids)){
  di<-subset(d, id_yr==ids[i])
  di$NSD<-nsd(di$x, di$y)
  full.df<-rbind(full.df, di)
}
head(full.df)

#Convert NSD values from meters to kilometers
full.df$NSD.km<-full.df$NSD*0.001
range(full.df$NSD.km) #0-293
head(full.df)

#Temporarily subset an individual to make sure that distance looks correct
temp<-subset(full.df,full.df$id_yr=="255_2020")
head(temp)
tail(temp)

#Determine the range of kilometers for each NSD value
#If an individual is between 0 and 1 km...then they are in the first km of migration
#If an individual is between 52 and 53...then they are in the 53rd km of migration
f<-full.df %>%
  group_by(gr=cut(NSD.km, breaks= seq(0, 293, by = 1)) )
head(f)

#Split the column with km range and create column with km marker
f$gr1 <- str_split_fixed(f$gr, ",", 2)[,2]
f$gr2 <- str_split_fixed(f$gr1, "]", 2)[,1]
f$km.mark<-as.numeric(f$gr2)
f$km.mark[is.na(f$km.mark)] <-0
names(f)
f<-f[c(1:19,25)]
table(is.na(f$km.mark))
range(f$km.mark)

#Ensure data is in correct format
f$date<-as.POSIXct(f$date, format="%Y-%m-%d %H:%M:%S")

#Now determine total length of migration for each individual
f<-f[order(f$id_yr,f$date),]

ids<-unique(f$id_yr)
df<-data.frame()

#Grabs the last known location of each deer
for (i in 1:length(ids)){
  temp<- subset(f, id_yr==ids[i])
  last <- tail(temp, n=1)
  df<- rbind(df, last)
}

names(df)
df<-df[c(1,15,20)]
colnames(df)<-c("id_yr","timing","distance")

library(Rmisc)
sum = summarySE(df,measurevar="distance",groupvars = "timing")
sum$true.ci<-sum$se*1.96
sum

range(df$distance) #Range of distance is 134-293 km

#Merge distance dataset with GPS dataset
r<-merge(f,df,by=c("id_yr"))
head(r)
range(r$distance)

#Clean up working environment
rm(list=ls()[! ls() %in% c("r")])

#### CALCULATE AVERAGE RATE OF MOVEMENT ####
r$Duration<-as.numeric(difftime(r$SpringME,r$SpringMS,units="days"))
r$rate.km.d<-r$distance/r$Duration #km/day

#Remove unnecessary columns and rename columns
names(r)
r<-r[c(1:20,22:24)]
colnames(r)<-c("id_yr","x","y","Lat","Long","date","id","migration","year","SpringMS","SpringME",
               "SpringMSJul","SpringMEJul","jul","timing","c.SpringMSJul","IRGVals","fittedVals",
               "maxIRGdate","km.mark","distance.km","duration.days","rate.km.day")
names(r)
id.list<-r[c(1,15,21)]
str(id.list)
id.list<-id.list[duplicated(id.list$id_yr)==FALSE,]  #gives first row of dataframe for each month
head(id.list)

#Export id_yr list with timing of migration and distance traveled
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Metadata")
write.csv(id.list,"IDList_TimingCat.csv")

#Remove duplicated values to look at summary of duration
dur<-r[duplicated(r$id_yr)==FALSE,]  #gives first row of dataframe for each month
head(dur)

sum = summarySE(dur,measurevar = c("duration.days"),groupvars = c("timing"))
sum$true.ci<-sum$se*1.96
sum

#Clean working environment
rm(list=ls()[! ls() %in% c("r")])

#Ensure that extracted IRG values look correct
#Remove any null values
range(r$maxIRGdate, na.rm=TRUE)
table(is.na(r$maxIRGdate))
names(r)
r<-r[!is.na(r[,17]),]
r<-r[!is.na(r[,18]),]
r<-r[!is.na(r[,19]),]
range(r$maxIRGdate)
range(r$IRGVals) #Should be 0-1
table(is.na(r$maxIRGdate)) #Should all be FALSE

#Calculate DFP (including absolute)
r$DFP<-r$jul-r$maxIRGdate
range(r$DFP)
r$abs.DFP<-abs(r$jul-r$maxIRGdate)
range(r$abs.DFP)

#Ensure that database is still consistent with sample size
u<-unique(r$id_yr)

#Merge data with compensation category
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Metadata")
id<-read.csv("IDList_BehavioralCompensation_20FEB2023.csv",header=TRUE,sep=",")
id<-id[c(2:3)]
colnames(id)<-c("id_yr","compensation")
head(id)

data<-merge(r,id,by=c("id_yr"))
data<-data[order(data$id_yr,data$km.mark),]
head(data)
table(data$id_yr)

#Now, export data as an R.data file and shapefile
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Metadata")
#save(r,file="RD2H_SpringMigrations_2011_thru_2020_WithTimingDistanceIRG_28DEC2021.RData")
save(data,file="RD2H_SpringMigrations_2011_thru_2020_WithTimingDistanceIRG_20FEB2023.RData")

#Make database with GPS points spatially aware
proj <- "+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs+towgs84=0,0,0"   #name a new proj
coordinates(r) <- c("x","y")
proj4string(r) <- CRS(proj)

fileSuffix="SpringMigrations_2011_thru_2020_UsedForAnalyses"
outname<- paste("RD2H",sep="_", fileSuffix)

writeOGR(r, "D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data/SpringMigrations", outname, driver="ESRI Shapefile",overwrite=TRUE)

#Clean working environment
rm(list=ls())

#### IMPORT CLEAN GPS COLLAR DATA WITH TIMING AND DISTANCE INFORMATION ####
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data/SpringMigrations")
load("RD2H_SpringMigrations_2011_thru_2020_WithTimingDistanceIRG_28DEC2021.RData")
head(r)

#Determine sample size of early, mid, and late migrants
early<-subset(r,r$timing=="early")
e1<-unique(early$id_yr) #n=47
e2<-unique(early$id) #n=32

mid<-subset(r,r$timing=="mid")
m1<-unique(mid$id_yr) #n=58
m2<-unique(mid$id) #n=40

late<-subset(r,r$timing=="late")
l1<-unique(late$id_yr) #n=47
l2<-unique(late$id) #n=34

#First, determine mean DFP for each id_yr and km marker 
# d1<-aggregate(r$DFP, by=list(r$id_yr,r$timing,r$km.mark), FUN=mean, na.rm=TRUE)
# colnames(d1)<-c("id_yr","timing","km","DFP")
# d1<-d1[order(d1$id_yr,d1$km),]
# head(d1)
# tail(d1)

#First, determine absolute DFP for each id_yr and km marker (abs DFP for surfing score) 
d1<-aggregate(r$abs.DFP, by=list(r$id_yr,r$timing,r$km.mark), FUN=mean, na.rm=TRUE)
colnames(d1)<-c("id_yr","timing","km","abs.DFP")
d1<-d1[order(d1$id_yr,d1$km),]
head(d1)
tail(d1)

#Second, determine mean surfing score each id_yr
d2<-aggregate(d1$abs.DFP, by=list(d1$id_yr,d1$timing), FUN=mean, na.rm=TRUE)
colnames(d2)<-c("id_yr","timing","abs.DFP")
d2<-d2[order(d2$id_yr),]
head(d2)

#Compare surfing scores among early, mid, and late migrants
library(Rmisc)
sum=summarySE(d2,measurevar = c("abs.DFP"),groupvars = c("timing"))
sum$true.ci<-sum$se*1.96
sum

mod1<-aov(d2$abs.DFP~d2$timing)
summary(mod1)

#Second, determine mean DFP for each id_yr and km marker
d3<-aggregate(r$DFP, by=list(r$id_yr,r$timing,r$km.mark), FUN=mean, na.rm=TRUE)
colnames(d3)<-c("id_yr","timing","km","DFP")
d3<-d3[order(d3$id_yr,d3$km),]
head(d3)
tail(d3)

#Third, determine mean DFP for each id_yr and km marker 
e<-summarySE(d3,measurevar = c("DFP"),groupvars = c("timing","km"))
e$true.ci<-e$se*1.96
e<-e[c("timing","km","DFP","true.ci")]
colnames(e)<-c("timing","km","DFP","CI_95")
e<-e[order(e$timing,e$km),]
head(e)
tail(e)

#Subset km with all three timing categories (early, mid, and late)
table(e$km)
e<-subset(e,e$km<=222) #Early, mid, and late migrants all make it to km 222
table(e$km)

head(e)
tail(e)

#Export data for Source Data File in Nature Communications
all<-e
all$timing<-factor(all$timing, levels=c("early","mid","late"))
all<-all[order(all$timing),]
head(all)
tail(all)

setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Submission/NatureCommunicationsTransfer/SourceData")
write.csv(all,"SourceData_DFP_DuringSpringMigration.csv")

#### EXPORT DATA FOR UO INFOGRAPHICS LAB ####
e.early<-subset(e,e$timing=="early")
e.mid<-subset(e,e$timing=="mid")
e.late<-subset(e,e$timing=="late")

#Predict DFP using loess function for early migrants
m.early <- loess(DFP ~ km,group=timing, data=e.early)

df <- data.frame(expand.grid(timing=unique(e.early$timing),
                             km=seq(from = min(e$km),to = max(e$km),length=500)))

df
range(df$km)

pred <- predict(m.early, newdata = df, se=TRUE)
y = pred$fit
ci <- pred$se.fit * 1.96
ymin = y - ci
ymax = y + ci

loess.early <- data.frame(timing=df$timing,x = df$km, y, ymin, ymax, se = pred$se.fit)
head(loess.early)

rm(df,pred,y,ci,ymin,ymax)

#Predict DFP using loess function for mid-migrants
m.mid <- loess(DFP ~ km,group=timing, data=e.mid)

df <- data.frame(expand.grid(timing=unique(e.mid$timing),
                             km=seq(from = min(e$km),to = max(e$km),length=500)))

df

pred <- predict(m.mid, newdata = df, se=TRUE)
y = pred$fit
ci <- pred$se.fit * 1.96
ymin = y - ci
ymax = y + ci

loess.mid <- data.frame(timing=df$timing,x = df$km, y, ymin, ymax, se = pred$se.fit)
head(loess.mid)

rm(df,pred,y,ci,ymin,ymax)

#Predict DFP using loess function for late migrants
m.late <- loess(DFP ~ km,group=timing, data=e.late)

df <- data.frame(expand.grid(timing=unique(e.late$timing),
                             km=seq(from = min(e$km),to = max(e$km),length=500)))

df

pred <- predict(m.late, newdata = df, se=TRUE)
y = pred$fit
ci <- pred$se.fit * 1.96
ymin = y - ci
ymax = y + ci

loess.late <- data.frame(timing=df$timing,x = df$km, y, ymin, ymax, se = pred$se.fit)
head(loess.late)

rm(pred,y,ci,ymin,ymax)

#Combine all prediction for early, mid, and late migrants into one dataframe
all.loess<-rbind(loess.early,loess.mid,loess.late)
head(all.loess)

p = all.loess %>% group_by(timing) %>% mutate(change = y - lag(y, default = y[1]))
p

#Export predictions
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/DataForInfoGraphics/Figure2")
write.csv(all.loess,"RawDataForFigure2.csv")

#Ensure that loess predictions look correct
a1<-ggplot(data=all.loess, aes(x=x, y=y,group=timing,colour=timing)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax, x=x, group=timing,colour=timing), alpha = 0.3)+
  geom_line(aes(x=x,y=y))+
  scale_colour_manual(values=c("#0000FF","#FFCC00",'#339900'))+
  scale_fill_manual(values=c("#0000FF","#FFCC00",'#339900'))+
  geom_hline(yintercept = 0, linetype="dashed", size=1)

a1

#### PLOT DFP FOR EACH KM ALONG THE MIGRATORY ROUTE ####
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/DataForInfoGraphics/Figure2")
jpeg("Ortega_et_al_Figure2.jpeg", width = 20, height = 8, units = 'in', res = 350)
library(ggplot2)

a1<-ggplot(data=e, aes(x=km, y=DFP,group=timing,colour=timing)) +
  geom_errorbar(aes(ymin=DFP-CI_95,ymax=DFP+CI_95),width=.15, size=0.4,colour="#666666",alpha=0.6)+
  geom_point(size=1.25,alpha=0.6)+
  geom_line(size=0.75,alpha=0.6)+
  geom_smooth(method="loess",data=e,aes(group=timing,colour=timing),size=1.5,se=FALSE)+
  geom_hline(yintercept = 0, linetype="dashed", size=1)

a1

b1<-a1 + scale_colour_manual(name  ="Departure from winter range",
                      breaks=c("early", "mid", "late"),
                      labels=c("Early departure", "Mid-departure", "Late departure"),
                      values=c("#9966CC","#339966","#FF9933"))
  

b1

d1 <- b1+xlab("Distance from winter range (km)") + ylab("Green wave surfing \n(days from peak IRG)")
d1

e1<- d1 + theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.line.x = element_line(colour = "black", size = 0.5),
                             axis.line.y = element_line(colour = "black", size = 0.5))
e1
f1 <-e1+ theme(axis.text.x=element_text(size=30,colour="black",family="serif"),
               axis.text.y=element_text(size=30,colour="black",family="serif",margin=margin(15,4,15,15)),
               axis.ticks.length=unit(.25, "cm"),
               axis.title.x=element_text(size=30,colour="black",family="serif",margin=margin(15,15,15,15)),
               axis.title.y=element_text(size=30,colour="black",family="serif",margin=margin(15,15,15,15)))
f1
g1<- f1 + theme (plot.background=element_blank(),
                 plot.margin = unit(c(.5,.5,0.005,.5),"cm"),
                 legend.position="none",
                 axis.line.x = element_line(colour = "black", size = 0.5),
                 axis.line.y = element_line(colour = "black", size = 0.5))

g1

dev.off()

# #Export figure for UO InfoGraphics
setwd("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/DataForInfoGraphics/Figure2")
ggsave("Ortega_et_al_2023_Figure2.PDF",g1, width = 20, height = 8, units = 'in')
