#The following is code for the manuscript
#"Global plant trait relationships extend to the climatic extremes of the tundra biome"
#Part 3: World Map
#20th Jan 2020

## COMBINE GEOSPATIAL INFORMATION FOR ALL OBSERVATIONS IN TRY FOR FOUR MAJOR TRAITS ##

#Detach packages####
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()

#Librarys####
library(maps)
library(mapdata)
library(mapproj)
library(raster)
library(dplyr)
library(ggplot2)
library(gridExtra)
`%notin%` <- function(x,y) !(x %in% y)

#read in TTT data
load("...") #Data available from github.com/TundraTraitTeam/TraitHub


ttt.sla <- subset(try.ttt.clean, TraitShort=="SLA" & Source=="TTT")
ttt.ht <- subset(try.ttt.clean, TraitShort=="PlantHeight" & Source=="TTT")
ttt.leafN <- subset(try.ttt.clean, TraitShort=="LeafN" & Source=="TTT")
ttt.ssd <- subset(try.ttt.clean, TraitShort=="StemSpecificDensity" & Source=="TTT")
ttt.sm <- subset(try.ttt.clean, TraitShort=="SeedMass" & Source=="TTT")
ttt.la <- subset(try.ttt.clean, TraitShort=="LeafArea" & Source=="TTT")
ttt.ldmc <- subset(try.ttt.clean, TraitShort=="LDMC" & Source=="TTT")

#read in TRY (cleaned)
ttt.try.sla <- subset(try.ttt.clean, TraitShort=="SLA" & Source=="TRY")
ttt.try.ht <- subset(try.ttt.clean, TraitShort=="PlantHeight" & Source=="TRY")
ttt.try.leafN <- subset(try.ttt.clean, TraitShort=="LeafN" & Source=="TRY")
ttt.try.ssd <- subset(try.ttt.clean, TraitShort=="StemSpecificDensity" & Source=="TRY")
ttt.try.sm <- subset(try.ttt.clean, TraitShort=="SeedMass" & Source=="TRY")
ttt.try.la <- subset(try.ttt.clean, TraitShort=="LeafArea" & Source=="TRY")
ttt.try.ldmc <- subset(try.ttt.clean, TraitShort=="LDMC" & Source=="TRY")

#All TTT and try
ttt.all.sla <- subset(try.ttt.clean, TraitShort=="SLA")
ttt.all.ht <- subset(try.ttt.clean, TraitShort=="PlantHeight")
ttt.all.leafN <- subset(try.ttt.clean, TraitShort=="LeafN")
ttt.all.ssd <- subset(try.ttt.clean, TraitShort=="StemSpecificDensity")
ttt.all.sm <- subset(try.ttt.clean, TraitShort=="SeedMass")
ttt.all.la <- subset(try.ttt.clean, TraitShort=="LeafArea")
ttt.all.ldmc <- subset(try.ttt.clean, TraitShort=="LDMC")

#Read in TRY (everything)

#read in observations of four traits separately
try.sla <- read.delim(file="...",stringsAsFactors = F, strip.white=T,header=T,sep="\t", skip=3)[1:11] #Data available from www.try-db.org 
try.sla$TraitShort <- "SLA"

try.leafN<-read.delim(file="...",stringsAsFactors = F, strip.white=T,header=T,sep="\t", skip=3)[1:11] #Data available from www.try-db.org 
try.leafN$TraitShort <- "LeafN"

try.ht<-read.delim(file="...",stringsAsFactors = F, strip.white=T,header=T,sep="\t", skip=3)[1:11] #Data available from www.try-db.org 
try.ht$TraitShort <- "PlantHeight"

try.ssd<-read.delim(file="...",stringsAsFactors = F, strip.white=T,header=T,sep="\t", skip=3)[1:11] #Data available from www.try-db.org 
try.ssd$TraitShort <- "StemSpecificDensity"

try.sm<-read.delim(file="...",stringsAsFactors = F, strip.white=T,header=T,sep="\t", skip=3)[1:11] #Data available from www.try-db.org 
try.sm$TraitShort <- "SeedMass"

try.la<-read.delim(file="...",stringsAsFactors = F, strip.white=T,header=T,sep="\t", skip=3)[1:11] #Data available from www.try-db.org 
try.la$TraitShort <- "LeafArea"

try.ldmc<-read.delim(file="...",stringsAsFactors = F, strip.white=T,header=T,sep="\t", skip=3)[1:11] #Data available from www.try-db.org 
try.ldmc$TraitShort <- "LDMC"


#combine all TRY traits
#currently written for seven traits but can delete number of traits as required
try.all<-rbind(try.sla, try.ht, try.leafN, try.ssd, try.sm, try.la, try.ldmc)
head(try.all)

#combine all TTT traits
#currently written for seven traits but can delete number of traits as required
ttt.all<-rbind(ttt.sla, ttt.ht, ttt.leafN, ttt.ssd, ttt.sm, ttt.la, ttt.ldmc)
head(ttt.all)

#combine all TTT.TRY traits
#currently written for seven traits but can delete number of traits as required
ttt.try.all<-rbind(ttt.try.sla, ttt.try.ht, ttt.try.leafN, ttt.try.ssd, ttt.try.sm, ttt.try.la, ttt.try.ldmc)
head(ttt.try.all)

#combine all TTT.TRY traits
#currently written for seven traits but can delete number of traits as required
ttt.all.all<-rbind(ttt.all.sla, ttt.all.ht, ttt.all.leafN, ttt.all.ssd, ttt.all.sm, ttt.all.la, ttt.all.ldmc)
head(ttt.try.all)

length(unique(ttt.all.all$AccSpeciesName))


#write.csv(try.all4,file="...") #the resulting file is too large for GitHub

#Denote which of these species are in our dataset
#Select species####

#1) TTT species list

load("...") #Data available from github.com/TundraTraitTeam/TraitHub
cover_species<-as.data.frame(unique(coverc.sub$Name))

names(cover_species)[1]<-"Name"

#2) Species from sites below 10oC July temp

ttt_only<-try.ttt %>%
  filter(Source == "TTT") %>%
  dplyr::select(Lon,Lat,AccSpeciesName) 

ttt_coords<-try.ttt %>%
  filter(Source == "TTT") %>%
  dplyr::select(Lon,Lat) %>%
  distinct()

ttt_only$latlon<-paste(ttt_only$Lat,ttt_only$Lon,sep="")

#Get temp data
July_temp <- raster("...") #Data available from chelsa-climate.org

#Extract cooridnates
ttt_climate <- raster::extract(July_temp, ttt_coords)

#Combine climate
ttt_coords$temp<-ttt_climate/10
ttt_coords$latlon<-paste(ttt_coords$Lat,ttt_coords$Lon,sep="")

#Combine with full object
ttt_only$Jul_temp<-ttt_coords$temp[match(ttt_only$latlon, ttt_coords$latlon)]

cover_species_ttt<-ttt_only %>%
  filter(Jul_temp <10) %>%
  dplyr::select(AccSpeciesName) %>%
  distinct(AccSpeciesName)

cover_species_ttt<-as.data.frame(cover_species_ttt)
names(cover_species_ttt)[1]<-"Name"

#3) ITEX Species list
load("...") #Data available from polardata.ca
cover_species_ITEX<-as.data.frame(unique(itex.vasccover.full$name))
names(cover_species_ITEX)[1]<-"Name"

#4) summit div species list
sdiv<-read.csv("...") #Data from Steinbauer M.J., Grytnes J.-A., Wipf S. (2018) 
#Accelerated increase in plant species richness on mountain summits is linked to warming. Nature v. 556, pages 231–234
cover_species_sdiv<-as.data.frame(unique(sdiv$Accepted_species_name))
names(cover_species_sdiv)[1]<-"Name"

#COMBINE species
species_list<-rbind(cover_species_ITEX, cover_species_ttt, cover_species, cover_species_sdiv)

species_list<-unique(species_list$Name)

#O's and 1's to denote which species are in our dataset (1's are species in our dataset)
#note that there are still some name correction issues that we should deal with before the final version (i.e. names that don't match between TRY and ITEX)
try.all$TundraSpecies<-0
try.all$TundraSpecies[try.all$AccSpeciesName %in% species_list]<-1

ttt.all$TundraSpecies<-0
ttt.all$TundraSpecies[ttt.all$AccSpeciesName %in% species_list]<-1

ttt.all<-ttt.all[ttt.all$AccSpeciesName %in% species_list,]

# Map Figure --------------------------------------------------------------
coords.try <- try.all
coords.ttt <- ttt.all
coords.ttt.try<-ttt.try.all

#Only use unique coordinates
coords.ttt$Lat<-round(coords.ttt$Lat/0.1)*0.1
coords.ttt$Lon<-round(coords.ttt$Lon/0.1)*0.1

coords.ttt$Latlon<-paste(coords.ttt$Lat,coords.ttt$Lon,sep = "_")
coords.ttt<-coords.ttt[!duplicated(coords.ttt$Latlon),]

coords.try$Lat<-round(coords.try$Lat/0.1)*0.1
coords.try$Lon<-round(coords.try$Lon/0.1)*0.1

coords.try$Latlon<-paste(coords.try$Lat,coords.try$Lon,sep = "_")
coords.try<-coords.try[!duplicated(coords.try$Latlon),]

coords.try$Lat<-as.numeric(as.character(coords.try$Lat))
coords.try$Lon<-as.numeric(as.character(coords.try$Lon))

#Remove incorrect coordinates
weird_lat<-unique(subset(coords.try,Lat<5&Lat>0&Lon<(52)&Lon>(45))[,c(5,6)]) #Somalia
weird_lat<-rbind(weird_lat,unique(subset(coords.try,Lat<3&Lat>0&Lon<7&Lon>0)[,c(5,6)])) #Nigeria
weird_lat<-rbind(weird_lat,unique(subset(coords.try,Lat<3&Lat>2&Lon<8&Lon>7)[,c(5,6)])) #Nigeria
weird_lat<-rbind(weird_lat,unique(subset(coords.try,Lat<60&Lat>10&Lon<(-20)&Lon>(-50))[,c(5,6)])) #Atlantic
weird_lat<-rbind(weird_lat,unique(subset(coords.try,Lat<(-12)&Lat>(-30)&Lon<(90)&Lon>(60))[,c(5,6)]))
coords.try<-subset(coords.try,Lat%notin%weird_lat$Lat | Lon%notin%weird_lat$Lon)

coords.try.tundra <- coords.try[coords.try$TundraSpecies ==1, ]


##Script for one map#####

#Different colours####

png(file="...",width=2500,height=1500)
par(mfrow=c(1,1), mar=c(0, 0, 0, 0), oma=c(0, 0, 0, 0), mgp=c(0, 0, 0))

# All traits map
map("world", col="black", proj="gall", param = 0, xlim=c(-180,180) ,ylim=c(-60,88), fill=FALSE, plot=TRUE, wrap=TRUE, boundary=TRUE, myborder = 0.01, lwd = 2)

All_Traits_Map <- mapproject(jitter(coords.try$Lon,amount = 1), jitter(coords.try$Lat,amount = 1), proj="", param = 0)
points(All_Traits_Map, pch=19, col="gray65", bg="#00000000", cex=3.5)

All_Traits_TRY <- subset(coords.ttt.try)
#Only use unique coordinates
All_Traits_TRY$Lat<-round(All_Traits_TRY$Lat,digits=1)
All_Traits_TRY$Lon<-round(All_Traits_TRY$Lon,digits=1)

All_Traits_TRY$Latlon<-paste(All_Traits_TRY$Lat,All_Traits_TRY$Lon,sep = "_")
All_Traits_TRY<-All_Traits_TRY[!duplicated(All_Traits_TRY$Latlon),]
#
All_Traits_TRY_Map <- mapproject(jitter(coords.try.tundra$Lon,amount=1.5), jitter(coords.try.tundra$Lat,amount=1.5), proj="", param = 0) 
points(All_Traits_TRY_Map, pch=21, col="#FF000000", bg="#FF7F0085", cex=5)

All_Traits_TTT_Map <- mapproject(jitter(coords.ttt$Lon,amount=1.5), jitter(coords.ttt$Lat,amount=1.5), proj="", param = 0) 
points(All_Traits_TTT_Map, pch=21, col="#9932CC00", bg="#9932CC97", cex=5)

leg.txt <- c("All TRY data","Tundra TRY data","TTT data")
legend("bottomleft", legend=leg.txt, x.intersp=1, y.intersp=1, col=c("gray55", "#FF7F00", "#9932CC"), pch=c(19,19,19), pt.cex=3, bty="n", cex=3)

dev.off() 

#Only TRY


#Draw distribution plots

#Combine try & ttt

try.all<-dplyr::select(try.all,AccSpeciesName,Lat,Lon,TraitShort,TundraSpecies)
try.all$Source<-"TRY_world"

ttt.all.all<-dplyr::select(ttt.all.all,AccSpeciesName,Lat,Lon,TraitShort,Source)
ttt.all.all$TundraSpecies<-"1"

all.all<-rbind(try.all,ttt.all.all)
tundra.TRY<-subset(all.all,TundraSpecies=="1")
all.all$Source<-"All"

tundra.TRY$Source<-"Tundra"

all.all<-rbind(all.all,try.all,tundra.TRY)


(Observations<-ggplot(all.all,aes(Lat,..count..,fill=Source,alpha=Source,linetype=Source))+
  geom_density(adjust=4)+
  theme_bw()+
  scale_x_reverse(limits = c(90,-90),
                  breaks=c(90,50,0,-50,-90))+
  scale_fill_manual(values = c("red","grey80","blue"),
                    breaks=c("All","TRY", "Tundra"),
                    labels=c("All","TRY","Tundra"))+
  scale_alpha_manual(values = c(0,0.8,0.8),
                     breaks=c("All","TRY", "Tundra"),
                     labels=c("All","TRY","Tundra"))+
  scale_linetype_manual(values = c("dashed","solid","solid"),
                        breaks=c("All","TRY", "Tundra"),
                        labels=c("All","TRY","Tundra"))+
  theme(legend.position="none")+
  labs(x="Latitude (Degrees)",y="Observations")+
  theme(axis.title.x=element_blank())+
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"), axis.text=element_text(size=12),axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))
  

#Extract unique species
test<-all.all %>%
  mutate(Lat = round(Lat,0)) %>%
  group_by(Lat,Source) %>%
  distinct(AccSpeciesName)

(Species<-ggplot(test,aes(x=Lat,..count..,fill=Source,alpha=Source,linetype=Source))+
  geom_density(adjust=4)+
  theme_bw()+
  scale_x_reverse(limits = c(90,-90),
                  breaks=c(90,50,0,-50,-90))+
  scale_fill_manual(values = c("red","grey80","blue"),
                    breaks=c("All","TRY", "Tundra"),
                    labels=c("All","TRY","Tundra"))+
  scale_alpha_manual(values = c(0,0.8,0.8),
                     breaks=c("All","TRY", "Tundra"),
                     labels=c("All","TRY","Tundra"))+
  scale_linetype_manual(values = c("dashed","solid","solid"),
                        breaks=c("All","TRY", "Tundra"),
                        labels=c("All","TRY","Tundra"))+
  theme(legend.position="none")+
  labs(x="Latitude (Degrees)",y="Species")+
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"), axis.text=element_text(size=12),axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))

  
pdf(file="...",width=6,height=2)
Observations
dev.off()

pdf(file="...",width=6,height=2)
Species
dev.off()

grid.arrange(Observations,Species,ncol=2)

#Create Whittacker plot####
#+ laoding_packages2, message = FALSE
#Remove NAs
all.all<-subset(all.all,!is.na(Lat)&!is.na(Lon))
coords<-dplyr::select(all.all,Lat,Lon)
names(coords)<-c("Lat","Lon")

require(slam)
require(rgdal)
require(biomod2)
require(ggplot2)
require(gridExtra)
library(stringr)

coords$Lat<-as.numeric(as.character(coords$Lat))
coords$Lon<-as.numeric(as.character(coords$Lon))
coords<-cbind(coords$Lon,coords$Lat)

#Convert to coordinates
coords<-coordinates(coords)

CHELSA_temp<- raster("...") #Data available from chelsa-climate.org
CHELSA_precip<- raster("...") #Data available from chelsa-climate.org
# This extracts the data for bioclim grid cells where the coordinates where Bnana is found
temp <- extract(CHELSA_temp, coords, df=TRUE) ## Temperature data are in °C × 10 to reduce the file sizes
precip <- extract(CHELSA_precip, coords, df=TRUE) ## Rainfall data are in ??


all.all$temp<-temp[,2]/10
all.all$precip<-precip[,2]

all.all$latlonsp<-paste(all.all$Lat,all.all$Lon,all.all$AccSpeciesName,sep=",")

all.all2<-all.all[!duplicated(all.all$latlonsp),]

#Plot

whit_plot<-ggplot(all.all2[all.all2$TundraSpecies==0,])+
  geom_point(aes(jitter(precip,amount=25),jitter(temp,amount=0.5)),colour="grey80",alpha=0.15,shape=16,size=1)+
  geom_point(data=all.all2[all.all2$TundraSpecies==1,],aes(jitter(precip,amount=25),jitter(temp,amount=0.5)),colour="blue",alpha=0.15,shape=16,size=1)+
  scale_y_reverse()+
  theme_bw()+
  labs(x="Mean annual precipitation (mm)",y = "Mean annual temperature (°C)")+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"), axis.text=element_text(size=12),axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))
  
png(file="...",width=300,height=300)
whit_plot
dev.off()

