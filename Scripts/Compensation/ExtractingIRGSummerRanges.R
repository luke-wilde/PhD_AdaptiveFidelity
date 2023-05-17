#Code for extracting date of peak IRG on summer seasonal ranges
#Written by Anna Ortega, Wyoming Cooperative Fish and Wildlife Research Unit
#December 28, 2021

#Load required libraries
library(rgdal)
library(stringr)
library(raster)

#Clean working environment
rm(list=ls())

#Import all summer ranges
d<-readOGR("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_CatchingTheGreenWave/Data/SummerRange/KUDs","RD2H_SummerSeasonalRanges_95KUD_2011_thru_2020")
head(d)
d$id_yr<-d$id

d<-as.data.frame(d)

sum=summarySE(d,measurevar = c("area"))
sum$true.ci<-sum$se*1.96
sum

#Create year
d$year<-str_split_fixed(d$id_yr, "_", 2)[,2]
head(d)

#Subset each unique year
d.11<-subset(d,d$year=="2011")
d.12<-subset(d,d$year=="2012")
d.14<-subset(d,d$year=="2014")
d.16<-subset(d,d$year=="2016")
d.17<-subset(d,d$year=="2017")
d.18<-subset(d,d$year=="2018")
d.19<-subset(d,d$year=="2019")
d.20<-subset(d,d$year=="2020")

#Extract date of peak IRG from Bischof calculations
setwd("D:/R-code/NDVI_Data_2020")
IRG<-stack("maxIRGdate.grd")

projIRG<-IRG[[1]]@crs

#Transform KUDs to projected coordinate system of IRG dataset
g.11<-spTransform(d.11, projIRG)
g.12<-spTransform(d.12, projIRG)
g.14<-spTransform(d.14, projIRG)
g.16<-spTransform(d.16, projIRG)
g.17<-spTransform(d.17, projIRG)
g.18<-spTransform(d.18, projIRG)
g.19<-spTransform(d.19, projIRG)
g.20<-spTransform(d.20, projIRG)

#Subset each year of IRG data from 2014-2020
IRG.11<-IRG[[12]]
IRG.12<-IRG[[13]]
IRG.14<-IRG[[15]]
IRG.16<-IRG[[17]]
IRG.17<-IRG[[18]]
IRG.18<-IRG[[19]]
IRG.19<-IRG[[20]]
IRG.20<-IRG[[21]]
names(IRG)

#Loop through unique animal-years and determine mean of max IRG date on winter range
#2011
id.yr<-unique(g.11$id_yr)

df.11 <- data.frame(id_yr = id.yr,
                    mean.IRG = NA)
df.11

for(i in 1:length(id.yr)){ 
  a <- subset(g.11, id_yr==id.yr[i])
  IRG.r1 <- crop(IRG.11, extent(a))
  IRG.r2 <- mask(IRG.r1, a)
  IRG.m<- cellStats(IRG.r2, stat='mean', na.rm=TRUE)
  df.11[i,"mean.IRG"] <- IRG.m
}

df.11

rm(a,IRG.r1,IRG.r2,IRG.m,id.yr)

#Loop through unique animal-years and determine mean of max IRG date on winter range
#2012
id.yr<-unique(g.12$id_yr)

df.12 <- data.frame(id_yr = id.yr,
                    mean.IRG = NA)
df.12

for(i in 1:length(id.yr)){ 
  a <- subset(g.12, id_yr==id.yr[i])
  IRG.r1 <- crop(IRG.12, extent(a))
  IRG.r2 <- mask(IRG.r1, a)
  IRG.m<- cellStats(IRG.r2, stat='mean', na.rm=TRUE)
  df.12[i,"mean.IRG"] <- IRG.m
}

df.12

rm(a,IRG.r1,IRG.r2,IRG.m,id.yr)

#Loop through unique animal-years and determine mean of max IRG date on winter range
#2014
id.yr<-unique(g.14$id_yr)

df.14 <- data.frame(id_yr = id.yr,
                    mean.IRG = NA)
df.14

for(i in 1:length(id.yr)){ 
  a <- subset(g.14, id_yr==id.yr[i])
  IRG.r1 <- crop(IRG.14, extent(a))
  IRG.r2 <- mask(IRG.r1, a)
  IRG.m<- cellStats(IRG.r2, stat='mean', na.rm=TRUE)
  df.14[i,"mean.IRG"] <- IRG.m
}

df.14

rm(a,IRG.r1,IRG.r2,IRG.m,id.yr)

#Loop through unique animal-years and determine mean of max IRG date on winter range
#2016
id.yr<-unique(g.16$id_yr)

df.16 <- data.frame(id_yr = id.yr,
                    mean.IRG = NA)
df.16

for(i in 1:length(id.yr)){ 
  a <- subset(g.16, id_yr==id.yr[i])
  IRG.r1 <- crop(IRG.16, extent(a))
  IRG.r2 <- mask(IRG.r1, a)
  IRG.m<- cellStats(IRG.r2, stat='mean', na.rm=TRUE)
  df.16[i,"mean.IRG"] <- IRG.m
}

df.16

rm(a,IRG.r1,IRG.r2,IRG.m,id.yr)

#Loop through unique animal-years and determine mean of max IRG date on winter range
#2017
id.yr<-unique(g.17$id_yr)

df.17 <- data.frame(id_yr = id.yr,
                    mean.IRG = NA)
df.17

for(i in 1:length(id.yr)){ 
  a <- subset(g.17, id_yr==id.yr[i])
  IRG.r1 <- crop(IRG.17, extent(a))
  IRG.r2 <- mask(IRG.r1, a)
  IRG.m<- cellStats(IRG.r2, stat='mean', na.rm=TRUE)
  df.17[i,"mean.IRG"] <- IRG.m
}

df.17

rm(a,IRG.r1,IRG.r2,IRG.m,id.yr)

#Loop through unique animal-years and determine mean of max IRG date on winter range
#2018
id.yr<-unique(g.18$id_yr)

df.18 <- data.frame(id_yr = id.yr,
                    mean.IRG = NA)
df.18

for(i in 1:length(id.yr)){ 
  a <- subset(g.18, id_yr==id.yr[i])
  IRG.r1 <- crop(IRG.18, extent(a))
  IRG.r2 <- mask(IRG.r1, a)
  IRG.m<- cellStats(IRG.r2, stat='mean', na.rm=TRUE)
  df.18[i,"mean.IRG"] <- IRG.m
}

df.18

rm(a,IRG.r1,IRG.r2,IRG.m,id.yr)

#Loop through unique animal-years and determine mean of max IRG date on winter range
#2019
id.yr<-unique(g.19$id_yr)

df.19 <- data.frame(id_yr = id.yr,
                    mean.IRG = NA)
df.19

for(i in 1:length(id.yr)){ 
  a <- subset(g.19, id_yr==id.yr[i])
  IRG.r1 <- crop(IRG.19, extent(a))
  IRG.r2 <- mask(IRG.r1, a)
  IRG.m<- cellStats(IRG.r2, stat='mean', na.rm=TRUE)
  df.19[i,"mean.IRG"] <- IRG.m
}

df.19

rm(a,IRG.r1,IRG.r2,IRG.m,id.yr)

#Loop through unique animal-years and determine mean of max IRG date on winter range
#2020
id.yr<-unique(g.20$id_yr)

df.20 <- data.frame(id_yr = id.yr,
                    mean.IRG = NA)
df.20

for(i in 1:length(id.yr)){ 
  a <- subset(g.20, id_yr==id.yr[i])
  IRG.r1 <- crop(IRG.20, extent(a))
  IRG.r2 <- mask(IRG.r1, a)
  IRG.m<- cellStats(IRG.r2, stat='mean', na.rm=TRUE)
  df.20[i,"mean.IRG"] <- IRG.m
}

df.20

all.summer.IRG<-rbind(df.11,df.12,df.14,df.16,df.17,df.18,df.19,df.20)

write.csv(all.summer.IRG,"D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/IRG/IRGSummerSeasonalRanges.csv")

#Clean working environment
rm(list=ls())

#Now extract IRG for population scale summer range for each year
d.11<-readOGR("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs","PopLevel_SummerSeasonalRanges_95KUD_2011")
d.12<-readOGR("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs","PopLevel_SummerSeasonalRanges_95KUD_2012")
d.14<-readOGR("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs","PopLevel_SummerSeasonalRanges_95KUD_2014")
d.16<-readOGR("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs","PopLevel_SummerSeasonalRanges_95KUD_2016")
d.17<-readOGR("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs","PopLevel_SummerSeasonalRanges_95KUD_2017")
d.18<-readOGR("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs","PopLevel_SummerSeasonalRanges_95KUD_2018")
d.19<-readOGR("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs","PopLevel_SummerSeasonalRanges_95KUD_2019")
d.20<-readOGR("D:/WYO_COOP/RedDesertToHobackDataMuleDeerProject/DissertationChapters/Chapter1_GWH_AmongMigratoryStrategies/Data/SummerRange/KUDs","PopLevel_SummerSeasonalRanges_95KUD_2020")

#Extract date of peak IRG from Bischof calculations
setwd("D:/R-code/NDVI_Data_2020")
IRG<-stack("maxIRGdate.grd")

projIRG<-IRG[[1]]@crs

#Transform KUDs to projected coordinate system of IRG dataset
g.11<-spTransform(d.11, projIRG)
g.12<-spTransform(d.12, projIRG)
g.14<-spTransform(d.14, projIRG)
g.16<-spTransform(d.16, projIRG)
g.17<-spTransform(d.17, projIRG)
g.18<-spTransform(d.18, projIRG)
g.19<-spTransform(d.19, projIRG)
g.20<-spTransform(d.20, projIRG)

#Subset each year of IRG data from 2014-2019
IRG.11<-IRG[[12]]
IRG.12<-IRG[[13]]
IRG.14<-IRG[[15]]
IRG.16<-IRG[[17]]
IRG.17<-IRG[[18]]
IRG.18<-IRG[[19]]
IRG.19<-IRG[[20]]
IRG.20<-IRG[[21]]
names(IRG)

IRG.r1 <- crop(IRG.20, extent(g.20))
IRG.r2 <- mask(IRG.r1, g.20)
plot(IRG.r2)
IRG.m<- cellStats(IRG.r2, stat='mean', na.rm=TRUE)
IRG.sd<- cellStats(IRG.r2, stat='sd', na.rm=TRUE)
IRG.n<-ncell(IRG.r2)

IRG.se<-IRG.sd/(sqrt(IRG.n))

#Mean date of peak IRG:
#160 for 2011
#132 for 2012
#151 for 2014
#141 for 2016
#148 for 2017
#137 for 2018
#154 for 2019
#155 for 2020
