#The following is code for the manuscript
#"Global plant trait relationships extend to the climatic extremes of the tundra biome"
#Part 2: Variance Partitioning
#20th Jan 2020

#Detach packages####
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()


#Load Packages####
library(lme4)
require(ape)
library(ggplot2)
library(gtable)
library(ggbiplot)
library(reshape2)
library(MuMIn)
library(stringr)
library(scales)
library(plyr)
library(raster)
library(dplyr)
library(mgcv)
require(vegan)
require(gridExtra)
library(truncnorm)
require(magrittr)
library(splines)
library(segmented)
`%notin%` <- function(x,y) !(x %in% y)
se <- function(x) sd(x)/sqrt(length(x))
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

####Import Data####
load("...") #Data available from github.com/TundraTraitTeam/TraitHub
try.ttt<-try.ttt.clean
try.ttt<-subset(try.ttt,TraitShort=="SeedMass"|TraitShort=="StemSpecificDensity"|TraitShort=="LeafN"|TraitShort=="SLA"|TraitShort=="LDMC"|TraitShort=="LeafArea"|TraitShort=="PlantHeight")
try.ttt[try.ttt$TraitShort=="StemSpecificDensity",2]<-"SSD"
#Traits of Interest:
#Seed dry mass
#Stem dry mass per stem fresh volume (stem specific density, SSD, wood density)
#Leaf nitrogen (N) content per leaf dry mass
#Leaf area per leaf dry mass (specific leaf area, SLA)
#Leaf area
#Leaf dry mass per leaf fresh mass (Leaf dry matter content, LDMC)
#Plant height


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
try.ttt_yes<-filter(try.ttt, AccSpeciesName %in% species_list)
try.ttt_no<-filter(try.ttt, AccSpeciesName %notin% species_list)

this_species_list<-as.data.frame(unique(try.ttt_yes$AccSpeciesName))

try.ttt_all<-try.ttt
try.ttt<-try.ttt_yes

####Clean out extreme values####

#Correct a few final synonyms
try.ttt[try.ttt$AccSpeciesName=="Myosotis alpina",]$AccSpeciesName<-"Myosotis alpestris"
try.ttt[try.ttt$AccSpeciesName=="Alnus alnobetula",]$AccSpeciesName<-"Alnus viridis"

allspecies[allspecies$TraitName=="Leaf nitrogen (N) content per leaf dry mass" & allspecies$AccSpeciesName=="Astragalus armatus",]$StdValue<-allspecies[allspecies$TraitName=="Leaf nitrogen (N) content per leaf dry mass" & allspecies$AccSpeciesName=="Astragalus armatus",]$StdValue/10
try.ttt<-subset(try.ttt,DataContributor!="John Dickie"|StdValue<50)

# #Exclude trees
PH<-subset(try.ttt,TraitShort=="PlantHeight")
PHs<-PH %>%
  group_by(AccSpeciesName) %>%
  summarise(height = mean(StdValue))
PM_species<-unique(subset(PHs,height>=5)$AccSpeciesName)

#Note tree species
PM_species
#Add boreal trees
PM_all<-c(PM_species,c("Alnus viridis", "Picea glauca", "Picea abies", "Picea mariana", "Betula pendula", "Betula glandulosa"))

try.ttt.PM<-subset(try.ttt,AccSpeciesName%in%PM_all)
try.ttt<-subset(try.ttt,AccSpeciesName%notin%PM_all)

#Add metadata
metadata<-unique(as.data.frame(try.ttt.full[,c(1,13)]))

try.ttt$F_Group<-metadata$GrowthForm[match(try.ttt$AccSpeciesName,metadata$Name)]

missing<-subset(try.ttt,is.na(F_Group))
try.ttt<-subset(try.ttt,!is.na(F_Group))

F_Groups<-read.csv("...") #Data from github.com/hjdthomas/Tundra_functional_groups

metadata2<-F_Groups[,c(2,12)]
missing$F_Group<-metadata2$F_Group[match(missing$AccSpeciesName,metadata2$AccSpeciesName)]

try.ttt<-rbind(try.ttt,missing)

try.ttt$unique<-paste(try.ttt$AccSpeciesName,try.ttt$TraitShort,try.ttt$DataContributor,try.ttt$StdValue,sep="_")

try.ttt_del<-try.ttt %>% 
  group_by(TraitShort,Genus) %>% 
  filter(StdValue>mean(StdValue)+sd(StdValue)*4|StdValue<mean(StdValue)-sd(StdValue)*4)

try.ttt_del$unique<-paste(try.ttt_del$AccSpeciesName,try.ttt_del$TraitShort,try.ttt_del$DataContributor,try.ttt_del$StdValue,sep="_")
try.ttt<-subset(try.ttt,unique%notin%try.ttt_del$unique)

try.ttt[try.ttt$TraitShort=="SLA",]$TraitShort<-"LMA"

trait.counts<-try.ttt %>% group_by(TraitShort,F_Group) %>% summarise((length(unique(StdValue))))

#Categorise by scale####

#Add local sites#
try.ttt$Local<-paste(round(try.ttt$Lat/0.1)*0.1,round(try.ttt$Lon/0.1)*0.1,sep=",")
try.ttt$Local<-paste(try.ttt$Lat,try.ttt$Lon,sep=",")

#Add region (1deg x 1deg####
try.ttt$Region<-paste(round(try.ttt$Lat/1)*1,round(try.ttt$Lon/1)*1,sep=",")

#Add continent
try.ttt$Continental<-paste(round(try.ttt$Lat/10)*10,round(try.ttt$Lon/10)*10,sep=",")

#Remove values that will be kicked out at site level
#try.ttt<- try.ttt %>% group_by(TraitShort,Local) %>% filter(length(StdValue)>9)


#Conduct variance partitioning####
Table.mm<-NULL
for (t in unique(try.ttt$TraitShort)){
  trait<-subset(try.ttt,TraitShort==t&!is.na(StdValue)&!is.na(Lat)&!is.na(Lon))
  #Remove species with few observations
  trait.nos <-trait %>% group_by(AccSpeciesName) %>% summarise(Number = length(unique(StdValue)))
  trait.nos<-subset(trait.nos,Number>9)$AccSpeciesName
  trait<-subset(trait,AccSpeciesName%in%trait.nos)
  trait<-subset(trait,StdValue!=0)
  
  #Conduct partition
  mm.tundra<-as.data.frame(varcomp(lme(log10(StdValue) ~1, random= ~1|F_Group/AccSpeciesName, data=trait,na.action=na.omit),1)[1:3])
  mm.tundra<- cbind(Taxonomic_hierarchy = rownames(mm.tundra), mm.tundra)
  names(mm.tundra)[2]<-c("Biome")
  
  ####Continental####
  points<-subset(as.data.frame(table(trait$Continental)),Freq>9)
  trait<-trait[which(trait$Continental %in% points$Var1),]
  trait<-subset(trait,Continental!="NA,NA")
  
  Out<-NULL
  for (i in unique(trait$Continental)){
    subset<-subset(trait,Continental==i)
    mm.continent<-tryCatch(as.data.frame(varcomp(lme(StdValue ~1, random= ~1|F_Group/AccSpeciesName, data=subset),1)[1:3]), error=function(e) data.frame(c(NA,NA,NA)))
    names(mm.continent)<-c("Value")
    row.names(mm.continent) <- c("F_Group","AccSpeciesName","Within")
    mm.continent$continent<-i
    mm.continent<- cbind(Taxonomic_hierarchy = rownames(mm.continent), mm.continent)
    Out<-rbind(Out, mm.continent)
  }
  continental.r2<-data.frame(Out)
  
  mm.continental.meanstable<-ddply(continental.r2,.(Taxonomic_hierarchy),summarise,
                                   Continental = mean(Value, na.rm=TRUE))
  mm.continental.meanstable$Taxonomic_hierarchy<- factor(mm.continental.meanstable$Taxonomic_hierarchy , levels = c("F_Group","AccSpeciesName","Within"))
  mm.continental.meanstable<-mm.continental.meanstable[order(mm.continental.meanstable$Taxonomic_hierarchy), ]
  
  ####Regional####
  points<-subset(as.data.frame(table(trait$Region)),Freq>9)
  trait<-trait[which(trait$Region %in% points$Var1),]
  trait<-subset(trait,Region!="NA,NA")
  
  Out<-NULL
  for (j in unique(trait$Region)){
    subset<-subset(trait,Region==j)
    mm.region<-tryCatch(as.data.frame(varcomp(lme(log10(StdValue) ~1, random= ~1|F_Group/AccSpeciesName, data=subset,na.action=na.omit),1)[1:3]), error=function(e) data.frame(c(NA,NA,NA)))
    names(mm.region)<-c("Value")
    row.names(mm.region) <- c("F_Group","AccSpeciesName","Within")
    mm.region$region<-j
    mm.region<- cbind(Taxonomic_hierarchy = rownames(mm.region), mm.region)
    Out<-rbind(Out, mm.region)
  }
  regional.r2<-data.frame(Out)
  
  mm.regional.meanstable<-ddply(regional.r2,.(Taxonomic_hierarchy),summarise,
                                Regional = mean(Value, na.rm=TRUE))
  mm.regional.meanstable$Taxonomic_hierarchy<- factor(mm.regional.meanstable$Taxonomic_hierarchy , levels = c("F_Group","AccSpeciesName","Within"))
  mm.regional.meanstable<-mm.regional.meanstable[order(mm.regional.meanstable$Taxonomic_hierarchy), ]
  
  ####Local####
  points<-subset(as.data.frame(table(trait$Local)),Freq>9)
  trait<-trait[which(trait$Local %in% points$Var1),]
  trait<-subset(trait,Local!="NA,NA")
  
  Out<-NULL
  for (k in unique(trait$Local)){
    subset<-subset(trait,Local==k)
    mm.local<-tryCatch(as.data.frame(varcomp(lme(log10(StdValue) ~1, random= ~1|F_Group/AccSpeciesName, data=subset,na.action=na.omit),1)[1:3]), error=function(e) data.frame(c(NA,NA,NA)))
    names(mm.local)<-c("Value")
    row.names(mm.local) <- c("F_Group","AccSpeciesName","Within")
    mm.local$local<-k
    mm.local<- cbind(Taxonomic_hierarchy = rownames(mm.local), mm.local)
    Out<-rbind(Out, mm.local)
  }
  local.r2<-data.frame(Out)
  
  mm.local.meanstable<-ddply(local.r2,.(Taxonomic_hierarchy),summarise,
                             Local = mean(Value, na.rm=TRUE))
  mm.local.meanstable$Taxonomic_hierarchy<- factor(mm.local.meanstable$Taxonomic_hierarchy , levels = c("F_Group","AccSpeciesName","Within"))
  mm.local.meanstable<-mm.local.meanstable[order(mm.local.meanstable$Taxonomic_hierarchy), ]
  
  
  trait.mm<-cbind(mm.tundra,mm.continental.meanstable,mm.regional.meanstable,mm.local.meanstable)
  trait.mm$Trait=t
  trait.mm<-trait.mm[,cbind(2,4,6,8,9)]
  Table.mm<-rbind(Table.mm,trait.mm)
}


####Coefficient of variation####
#cv=sd/mean


Table.cv<-NULL
for (t in unique(try.ttt$TraitShort)){
  trait<-subset(try.ttt,TraitShort==t&!is.na(StdValue)&!is.na(Lat)&!is.na(Lon))
  #Remove species with few observations
  trait.nos <-trait %>% group_by(AccSpeciesName) %>% summarise((length(unique(StdValue))))
  names(trait.nos)[2]<-"Number"
  trait.nos<-subset(trait.nos,Number>9)$AccSpeciesName
  trait<-subset(trait,AccSpeciesName%in%trait.nos)
  
  #Conduct partition
  cv.tundra<-as.data.frame(sd(trait$StdValue)/mean(trait$StdValue))
  
  ####Continental####
  points<-subset(as.data.frame(table(trait$Continental)),Freq>9)
  trait<-trait[which(trait$Continental %in% points$Var1),]
  trait<-subset(trait,Continental!="NA,NA")
  
  Out<-NULL
  for (j in unique(trait$Continental)){
    subset<-subset(trait,Continental==j)
    cv<-as.data.frame(sd(subset$StdValue)/mean(subset$StdValue))
    Out<-rbind(Out, cv)
  }
  cv.continent<-data.frame(Out)
  names(cv.continent)<-"Value"
  
  cv.continental.meanstable<-as.data.frame(mean(cv.continent$Value))
  
  ####Regional####
  points<-subset(as.data.frame(table(trait$Region)),Freq>9)
  trait<-trait[which(trait$Region %in% points$Var1),]
  trait<-subset(trait,Region!="NA,NA")
  
  Out<-NULL
  for (j in unique(trait$Region)){
    subset<-subset(trait,Region==j)
    cv<-as.data.frame(sd(subset$StdValue)/mean(subset$StdValue))
    Out<-rbind(Out, cv)
  }
  cv.region<-data.frame(Out)
  names(cv.region)<-"Value"
  
  cv.regional.meanstable<-as.data.frame(mean(cv.region$Value))
  
  ####Local####
  points<-subset(as.data.frame(table(trait$Local)),Freq>9)
  trait<-trait[which(trait$Local %in% points$Var1),]
  trait<-subset(trait,Local!="NA,NA")
  
  Out<-NULL
  for (k in unique(trait$Local)){
    subset<-subset(trait,Local==k)
    cv<-as.data.frame(sd(subset$StdValue)/mean(subset$StdValue))
    Out<-rbind(Out, cv)
  }
  cv.local<-data.frame(Out)
  names(cv.local)<-"Value"
  
  cv.local.meanstable<-as.data.frame(mean(cv.local$Value))
  
  trait.cv<-NULL
  trait.cv<-cbind(cv.tundra,cv.continental.meanstable,cv.regional.meanstable,cv.local.meanstable)
  trait.cv$Trait<-t
  names(trait.cv)<-c("4","3","2","1","Trait")
  Table.cv<-rbind(Table.cv,trait.cv)
}

Table.cv<-melt(Table.cv,id=c("Trait"))




#Reshape data####

Table.mm$Taxonomic_hierarchy<-(rep(c("F_Group","AccSpeciesName","Within"),7))
Table.mm<-melt(Table.mm, id=c("Trait","Taxonomic_hierarchy"))
names(Table.mm)[3]<-"Scale"
Table.mm$value<-as.numeric(as.character(Table.mm$value))
Table.mm<-subset(Table.mm, !is.na(value))


####Individual Trait Figures####
#Add up to 100%
plots_100<-NULL
for(i in unique(Table.mm$Trait)){
  mm.trait<-subset(Table.mm, Trait==i)
  mm.trait$Taxonomic_hierarchy<-as.character(mm.trait$Taxonomic_hierarchy)
  mm.trait[mm.trait$Taxonomic_hierarchy=="AccSpeciesName",]$Taxonomic_hierarchy<-"Species"
  mm.trait[mm.trait$Taxonomic_hierarchy=="F_Group",]$Taxonomic_hierarchy<-"Functional Group"
  mm.trait[mm.trait$Taxonomic_hierarchy=="Within",]$Taxonomic_hierarchy<-"Within"
  mm.trait$Scale<-as.character(mm.trait$Scale)
  mm.trait[mm.trait$Scale=="Local",]$Scale<-1
  mm.trait[mm.trait$Scale=="Regional",]$Scale<-2
  mm.trait[mm.trait$Scale=="Continental",]$Scale<-3
  mm.trait[mm.trait$Scale=="Biome",]$Scale<-4
  mm.trait$Scale<-as.numeric(mm.trait$Scale)
  mm.trait$order<-rep((c(1:3)),4)
  mm.trait<-mm.trait[order(mm.trait$order),]
  Table.cv.trait<-subset(Table.cv,Trait==i)
  mm.trait$cv<-Table.cv.trait$value[match(mm.trait$Scale,Table.cv.trait$variable)]
  mm.trait$value<-mm.trait$value*100
  
  plot<-ggplot(mm.trait, aes(x=Scale, y= value, fill=Taxonomic_hierarchy,linetype=Taxonomic_hierarchy)) + 
    geom_bar(stat="identity",colour="black")+
    scale_x_continuous(breaks=c(1,2,3,4), labels=c("Local", "Region", "Continent", "Biome"))+
    labs(x=NULL,
         y=NULL)+
    theme_bw()+
    #ylim(0,ifelse(max(mm.trait$value)>3,5.2,ifelse(max(mm.trait$value)>0.5,2.7,0.6)))+
    ggtitle(i)+
    theme(plot.title = element_text(face="bold"))+
    theme(plot.title = element_text(size=14, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.text=element_text(size=11), axis.title=element_text(size=14),axis.text.x=element_text(angle=-45,hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    scale_linetype_manual(values=c("solid", "solid","solid"),
                          guide=FALSE)+
    scale_fill_manual(values=alpha(c("darkgreen", "#F8766D","#619CFF" ),c(0.9,0.9,0.9)),
                      guide=FALSE,
                      name="Taxonomic Hierarchy",
                      breaks=c("Within", "Species","Functional Group"),
                      labels=c("Within species / individuals","Species","Functional Group"))
  plots_100[[i]]<-plot
}

png(file="Traits_1/var_part_pc.png", width = 500, height = 500)
grid.arrange(plots_100[["LeafArea"]],plots_100[["SeedMass"]],plots_100[["PlantHeight"]],plots_100[["LMA"]],plots_100[["LeafN"]],plots_100[["LDMC"]],ncol=3)
dev.off()

####Individual Trait Figures####
#Add up to CV
plots_cv<-NULL
for(i in unique(Table.mm$Trait)){
  mm.trait<-subset(Table.mm, Trait==i)
  mm.trait$Taxonomic_hierarchy<-as.character(mm.trait$Taxonomic_hierarchy)
  mm.trait[mm.trait$Taxonomic_hierarchy=="AccSpeciesName",]$Taxonomic_hierarchy<-"Species"
  mm.trait[mm.trait$Taxonomic_hierarchy=="F_Group",]$Taxonomic_hierarchy<-"Functional Group"
  mm.trait[mm.trait$Taxonomic_hierarchy=="Within",]$Taxonomic_hierarchy<-"Within"
  mm.trait$Scale<-as.character(mm.trait$Scale)
  mm.trait[mm.trait$Scale=="Local",]$Scale<-1
  mm.trait[mm.trait$Scale=="Regional",]$Scale<-2
  mm.trait[mm.trait$Scale=="Continental",]$Scale<-3
  mm.trait[mm.trait$Scale=="Biome",]$Scale<-4
  mm.trait$Scale<-as.numeric(mm.trait$Scale)
  mm.trait$order<-rep((c(1:3)),4)
  mm.trait<-mm.trait[order(mm.trait$order),]
  Table.cv.trait<-subset(Table.cv,Trait==i)
  mm.trait$cv<-Table.cv.trait$value[match(mm.trait$Scale,Table.cv.trait$variable)]
  mm.trait$value<-mm.trait$value*mm.trait$cv
  
  plot<-ggplot(mm.trait, aes(x=Scale, y= value, fill=Taxonomic_hierarchy,linetype=Taxonomic_hierarchy)) + 
    geom_bar(stat="identity",colour="black")+
    scale_x_continuous(breaks=c(1,2,3,4), labels=c("Local", "Region", "Continent", "Biome"))+
    labs(x=NULL,
         y=NULL)+
    theme_bw()+
    #ylim(0,ifelse(max(mm.trait$value)>3,5.2,ifelse(max(mm.trait$value)>0.5,2.7,0.6)))+
    ggtitle(i)+
    theme(plot.title = element_text(face="bold"))+
    theme(plot.title = element_text(size=14, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.text=element_text(size=11), axis.title=element_text(size=14),axis.text.x=element_text(angle=-45,hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    scale_linetype_manual(values=c("solid", "solid","solid"),
                          guide=FALSE)+
    scale_fill_manual(values=alpha(c("darkgreen", "#F8766D","#619CFF" ),c(0.9,0.9,0.9)),
                      guide=FALSE,
                      name="Taxonomic Hierarchy",
                      breaks=c("Within", "Species","Functional Group"),
                      labels=c("Within species / individuals","Species","Functional Group"))
  plots_cv[[i]]<-plot
}

png(file="Traits_1/var_part_cv.png", width = 500, height = 500)
grid.arrange(plots_cv[["LeafArea"]],plots_cv[["SeedMass"]],plots_cv[["PlantHeight"]],plots_cv[["LMA"]],plots_cv[["LeafN"]],plots_cv[["LDMC"]],ncol=3)
dev.off()

#Combine figures just for tundra

Table.mm$Scale<-as.character(Table.mm$Scale)
Table.mm[Table.mm$Scale=="Local",]$Scale<-1
Table.mm[Table.mm$Scale=="Regional",]$Scale<-2
Table.mm[Table.mm$Scale=="Continental",]$Scale<-3
Table.mm[Table.mm$Scale=="Biome",]$Scale<-4
names(Table.cv)[c(2,3)]<-c("Scale","cv")

Table.mm<-merge(Table.mm,Table.cv, by = c('Trait', 'Scale'))
tundra.only<-subset(Table.mm,Scale=="4")
tundra.only$value.cv<-tundra.only$value*tundra.only$cv
tundra.only$value.100<-tundra.only$value*100
tundra.only<-subset(tundra.only,Trait!="SSD")

tundra.only$Trait <- factor(tundra.only$Trait, levels = c("LeafArea","SeedMass","PlantHeight","LMA","LeafN","LDMC"))
tundra.only$Taxonomic_hierarchy <- factor(tundra.only$Taxonomic_hierarchy, levels = c("Within","AccSpeciesName","F_Group"))

#Mean of variance explained
mean(tundra.only[tundra.only$Taxonomic_hierarchy=="Within",]$value)
mean(tundra.only[tundra.only$Taxonomic_hierarchy=="AccSpeciesName",]$value)
mean(tundra.only[tundra.only$Taxonomic_hierarchy=="F_Group",]$value)

(a<-ggplot(tundra.only, aes(x=Trait, y= value.100, fill=factor(Taxonomic_hierarchy))) + 
    geom_bar(stat="identity",colour="black")+
    theme_bw()+
    xlab("\nTrait")+
    ylab("Variation explained (%)\n")+
    scale_fill_manual(values=alpha(c("#619CFF" ,"#F8766D","gold1"),c(0.9,0.9,0.9)),
                      name="Variance Hierarchy",
                      guide=FALSE,
                      breaks=c("Within", "AccSpeciesName","F_Group"),
                      labels=c("Within species","Between Species","Between Functional Groups"))+
    scale_x_discrete(labels = c("Leaf\nArea", "Seed\nMass", "Height", "LMA", "Leaf N", "LDMC"))+
    theme(plot.title = element_text(size=14, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.text=element_text(size=11), axis.title=element_text(size=14),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")))

(b<-ggplot(tundra.only, aes(x=Trait, y= value.cv, fill=Taxonomic_hierarchy)) + 
    geom_bar(stat="identity",colour="black")+
    theme_bw()+
    xlab("\nTrait")+
    ylab("Coefficient of variation\n")+
    scale_fill_manual(values=alpha(c("#619CFF" ,"#F8766D","gold1"),c(0.9,0.9,0.9)),
                      name="Variance Hierarchy\n",
                      breaks=c("Within", "AccSpeciesName","F_Group"),
                      labels=c("Within species","Between Species","Between Functional Groups"))+
    theme(plot.margin=unit(c(0.5,0.5,0,0.5), "cm"))+
    theme(legend.title = element_text(face="bold", size = 13),
          legend.position=c(0.70, 0.80),
          legend.text=element_text(size=12))+
    scale_x_discrete(labels = c("Leaf\nArea", "Seed\nMass", "Height", "LMA", "Leaf N", "LDMC"))+
    theme(plot.title = element_text(size=14, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.text=element_text(size=11), axis.title=element_text(size=14),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")))

a <- ggplot_gtable(ggplot_build(a))
b <- ggplot_gtable(ggplot_build(b))
maxWidth = unit.pmax(a$widths[2:3], b$widths[2:3])
a$widths[2:3] <- maxWidth
b$widths[2:3] <- maxWidth


pdf(file="scripts/users/hthomas/Output_Images/Traits_1/tundra_bars.pdf", width = 10, height = 4.5)
grid.arrange(a,b,ncol=2)
dev.off()

ddply(tundra.only,.(Taxonomic_hierarchy),summarise,
      mean = mean(value.100))

alt_sums<-subset(tundra.only,Taxonomic_hierarchy!="Within")
alt_values<-ddply(alt_sums,.(Trait),summarise,
                  species = sum(value.100))
mean(alt_values$species)

tundra.only$Trait <- factor(tundra.only$Trait, levels = c("LeafN","PlantHeight","LDMC","LMA","LeafArea","SeedMass"))
tundra.only$Taxonomic_hierarchy <- factor(tundra.only$Taxonomic_hierarchy, levels = c("Within","AccSpeciesName","F_Group"))

#Figure for ArcticNet#
ggplot(tundra.only[tundra.only$Taxonomic_hierarchy=="Within",], aes(x=Trait, y= value)) + 
  geom_bar(stat="identity", fill="white", colour="black")+
  theme_bw()+
  scale_y_continuous(limits = c(0, 0.72), breaks = seq(0,0.8, by = 0.1))+
  xlab("\nTrait")+
  ylab("ITV (%)\n")+
  scale_fill_manual(values=alpha(c("#619CFF" ,"#F8766D","darkgreen"),c(0.9,0.9,0.9)),
                    name="Taxonomic Hierarchy",
                    guide=FALSE,
                    breaks=c("Within", "AccSpeciesName","F_Group"),
                    labels=c("Within species","Between Species","Between Functional Groups"))+
  theme(plot.title = element_text(size=14, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.text=element_text(size=11), axis.title=element_text(size=14),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))

#Continuous Scale####
#Calculate distance between coordinates
library("geosphere")

#10 species / sites####

#Remove NAs
try.ttt2<-try.ttt
try.ttt<-subset(try.ttt,!is.na(Lat)&!is.na(Lon))

#Remove sites that only <10 species 
one.species<- try.ttt %>% 
  group_by(TraitShort,Local) %>% 
  summarise(n.sp = length(unique(AccSpeciesName))) %>%
  filter(n.sp<10)

try.ttt$Site.Trait<-paste(try.ttt$TraitShort,try.ttt$Local)
one.species$Site.Trait<-paste(one.species$TraitShort,one.species$Local)
one.species.list<-as.list(one.species$Site.Trait)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait%notin%one.species.list)

#& species <10obs/site
one.obs<- try.ttt %>% 
  group_by(TraitShort,Local,AccSpeciesName) %>% 
  summarise(n.obs = length(AccSpeciesName)) %>%
  filter(n.obs<10)

try.ttt$Site.Trait.Sp<-paste(try.ttt$TraitShort,try.ttt$Local,try.ttt$AccSpeciesName)
one.obs$Site.Trait.Sp<-paste(one.obs$TraitShort,one.obs$Local,one.obs$AccSpeciesName)
one.obs.list<-as.list(one.obs$Site.Trait.Sp)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait.Sp%notin%one.obs.list)


#Subset trait
all_traits<-NULL
for (t in unique(try.ttt$TraitShort)){
  trait<-subset(try.ttt,TraitShort==t)
  
  #pick unique site
  
  distance_frame<-trait
  all_sites<-NULL
  
  for(i in unique(trait$Local)){
    distance_frame<-trait
    #Identify coordinates
    pt<-as.matrix(subset(distance_frame,Local==i)[1,c(4,5)])
    #Matrix of all coordinates
    pts<-as.matrix(distance_frame[,c(4,5)])
    #Calculate distance
    distance_frame$distance<-round(spDistsN1(pts, pt, longlat = TRUE),5)
    
    #create object of lenght of unique distances
    length<-length(unique(distance_frame$distance))
    
    Scale_table<-NULL
    Out<-NULL
    repeat {
      #Add next site
      subset<-subset(distance_frame,distance==min(distance))
      distance_frame<-subset(distance_frame,distance!=min(distance))
      Scale_table<-rbind(Scale_table,subset)
      
      #Sample no. values in site to remove additional data problem
      no_samples<-nrow(pt)
      #Number of species in analysis
      no_species<-length(unique(Scale_table$AccSpeciesName))
      #Random sample from Scale_table no. values of initial site
      #Scale_table[sample(nrow(Scale_table), no_samples), ]
      
      
      #Conduct partition
      mm.scale<-tryCatch(as.data.frame(varcomp(lme(log10(StdValue) ~1, random= ~1|AccSpeciesName, data=Scale_table,na.action=na.omit),1)[1:2]), error=function(e) data.frame(c(NA,NA)))
      names(mm.scale)<-c("Value")
      cv<-as.data.frame(sd(Scale_table$StdValue)/mean(Scale_table$StdValue))
      names(cv)<-c("CV")
      mm.scale$scale<-max(Scale_table$distance)
      mm.scale$ndata<-nrow(Scale_table)
      mm.scale<- cbind(Taxonomic_hierarchy = rownames(mm.scale), cv,mm.scale,i,no_species)
      Out<-rbind(Out,mm.scale)
      if (length(unique(Scale_table$distance)) == length){
        break
      }
    }
    site_scale<-Out
    all_sites<-rbind(all_sites,site_scale)
  }
  all_sites$trait<-t
  all_traits<-rbind(all_traits,all_sites)
}

save(all_traits,file="...")

#Continuous Scale - larger dataset
#Calculate distance between coordinates
#5 species / sites####

#Remove NAs
try.ttt<-try.ttt2
try.ttt<-subset(try.ttt,!is.na(Lat)&!is.na(Lon))

#Remove sites that only <5 species 
one.species<- try.ttt %>% 
  group_by(TraitShort,Local) %>% 
  summarise(n.sp = length(unique(AccSpeciesName))) %>%
  filter(n.sp<5)

try.ttt$Site.Trait<-paste(try.ttt$TraitShort,try.ttt$Local)
one.species$Site.Trait<-paste(one.species$TraitShort,one.species$Local)
one.species.list<-as.list(one.species$Site.Trait)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait%notin%one.species.list)

#& species <5obs/site
one.obs<- try.ttt %>% 
  group_by(TraitShort,Local,AccSpeciesName) %>% 
  summarise(n.obs = length(AccSpeciesName)) %>%
  filter(n.obs<5)

try.ttt$Site.Trait.Sp<-paste(try.ttt$TraitShort,try.ttt$Local,try.ttt$AccSpeciesName)
one.obs$Site.Trait.Sp<-paste(one.obs$TraitShort,one.obs$Local,one.obs$AccSpeciesName)
one.obs.list<-as.list(one.obs$Site.Trait.Sp)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait.Sp%notin%one.obs.list)

#Subset trait
all_traits<-NULL
for (t in unique(try.ttt$TraitShort)){
  trait<-subset(try.ttt,TraitShort==t)
  
  #pick unique site
  
  distance_frame<-trait
  all_sites<-NULL
  
  for(i in unique(trait$Local)){
    distance_frame<-trait
    #Identify coordinates
    pt<-as.matrix(subset(distance_frame,Local==i)[1,c(4,5)])
    #Matrix of all coordinates
    pts<-as.matrix(distance_frame[,c(4,5)])
    #Calculate distance
    distance_frame$distance<-round(spDistsN1(pts, pt, longlat = TRUE),5)
    
    #create object of lenght of unique distances
    length<-length(unique(distance_frame$distance))
    
    Scale_table<-NULL
    Out<-NULL
    repeat {
      #Add next site
      subset<-subset(distance_frame,distance==min(distance))
      distance_frame<-subset(distance_frame,distance!=min(distance))
      Scale_table<-rbind(Scale_table,subset)
      
      #Sample no. values in site to remove additional data problem
      no_samples<-nrow(pt)
      #Number of species in analysis
      no_species<-length(unique(Scale_table$AccSpeciesName))
      #Random sample from Scale_table no. values of initial site
      #Scale_table[sample(nrow(Scale_table), no_samples), ]
      
      
      #Conduct partition
      mm.scale<-tryCatch(as.data.frame(varcomp(lme(log10(StdValue) ~1, random= ~1|AccSpeciesName, data=Scale_table,na.action=na.omit),1)[1:2]), error=function(e) data.frame(c(NA,NA)))
      names(mm.scale)<-c("Value")
      cv<-as.data.frame(sd(Scale_table$StdValue)/mean(Scale_table$StdValue))
      names(cv)<-c("CV")
      mm.scale$scale<-max(Scale_table$distance)
      mm.scale$ndata<-nrow(Scale_table)
      mm.scale<- cbind(Taxonomic_hierarchy = rownames(mm.scale), cv, mm.scale,i,no_species)
      Out<-rbind(Out,mm.scale)
      if (length(unique(Scale_table$distance)) == length){
        break
      }
    }
    site_scale<-Out
    all_sites<-rbind(all_sites,site_scale)
  }
  all_sites$trait<-t
  all_traits<-rbind(all_traits,all_sites)
}

save(all_traits,file="...")

#3 species / sites####
try.ttt<-try.ttt2
try.ttt<-subset(try.ttt,!is.na(Lat)&!is.na(Lon))

one.species<- try.ttt %>% 
  group_by(TraitShort,Local) %>% 
  summarise(n.sp = length(unique(AccSpeciesName))) %>%
  filter(n.sp<3)

try.ttt$Site.Trait<-paste(try.ttt$TraitShort,try.ttt$Local)
one.species$Site.Trait<-paste(one.species$TraitShort,one.species$Local)
one.species.list<-as.list(one.species$Site.Trait)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait%notin%one.species.list)

#& species <5obs/site
one.obs<- try.ttt %>% 
  group_by(TraitShort,Local,AccSpeciesName) %>% 
  summarise(n.obs = length(AccSpeciesName)) %>%
  filter(n.obs<3)

try.ttt$Site.Trait.Sp<-paste(try.ttt$TraitShort,try.ttt$Local,try.ttt$AccSpeciesName)
one.obs$Site.Trait.Sp<-paste(one.obs$TraitShort,one.obs$Local,one.obs$AccSpeciesName)
one.obs.list<-as.list(one.obs$Site.Trait.Sp)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait.Sp%notin%one.obs.list)

#Subset trait
all_traits<-NULL
for (t in unique(try.ttt$TraitShort)){
  trait<-subset(try.ttt,TraitShort==t)
  
  #pick unique site
  
  distance_frame<-trait
  all_sites<-NULL
  
  for(i in unique(trait$Local)){
    distance_frame<-trait
    #Identify coordinates
    pt<-as.matrix(subset(distance_frame,Local==i)[1,c(4,5)])
    #Matrix of all coordinates
    pts<-as.matrix(distance_frame[,c(4,5)])
    #Calculate distance
    distance_frame$distance<-round(spDistsN1(pts, pt, longlat = TRUE),5)
    
    #create object of lenght of unique distances
    length<-length(unique(distance_frame$distance))
    
    Scale_table<-NULL
    Out<-NULL
    repeat {
      #Add next site
      subset<-subset(distance_frame,distance==min(distance))
      distance_frame<-subset(distance_frame,distance!=min(distance))
      Scale_table<-rbind(Scale_table,subset)
      
      #Sample no. values in site to remove additional data problem
      no_samples<-nrow(pt)
      #Number of species in analysis
      no_species<-length(unique(Scale_table$AccSpeciesName))
      #Random sample from Scale_table no. values of initial site
      #Scale_table[sample(nrow(Scale_table), no_samples), ]
      
      
      #Conduct partition
      mm.scale<-tryCatch(as.data.frame(varcomp(lme(log10(StdValue) ~1, random= ~1|AccSpeciesName, data=Scale_table,na.action=na.omit),1)[1:2]), error=function(e) data.frame(c(NA,NA)))
      names(mm.scale)<-c("Value")
      cv<-as.data.frame(sd(Scale_table$StdValue)/mean(Scale_table$StdValue))
      names(cv)<-c("CV")
      mm.scale$scale<-max(Scale_table$distance)
      mm.scale$ndata<-nrow(Scale_table)
      mm.scale<- cbind(Taxonomic_hierarchy = rownames(mm.scale), cv, mm.scale,i,no_species)
      Out<-rbind(Out,mm.scale)
      if (length(unique(Scale_table$distance)) == length){
        break
      }
    }
    site_scale<-Out
    all_sites<-rbind(all_sites,site_scale)
  }
  all_sites$trait<-t
  all_traits<-rbind(all_traits,all_sites)
}

save(all_traits,file="...")

####Arctic Only####
#REMOVE NON-ARCTIC DATA
try.ttt<-try.ttt2
try.ttt<-subset(try.ttt,!is.na(Lat)&!is.na(Lon))

try.ttt<-subset(try.ttt,Lat>=60)

#Remove sites that only <3 species 
one.species<- try.ttt %>% 
  group_by(TraitShort,Local) %>% 
  summarise(n.sp = length(unique(AccSpeciesName))) %>%
  filter(n.sp<5)

try.ttt$Site.Trait<-paste(try.ttt$TraitShort,try.ttt$Local)
one.species$Site.Trait<-paste(one.species$TraitShort,one.species$Local)
one.species.list<-as.list(one.species$Site.Trait)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait%notin%one.species.list)

#& species <5obs/site
one.obs<- try.ttt %>% 
  group_by(TraitShort,Local,AccSpeciesName) %>% 
  summarise(n.obs = length(AccSpeciesName)) %>%
  filter(n.obs<5)

try.ttt$Site.Trait.Sp<-paste(try.ttt$TraitShort,try.ttt$Local,try.ttt$AccSpeciesName)
one.obs$Site.Trait.Sp<-paste(one.obs$TraitShort,one.obs$Local,one.obs$AccSpeciesName)
one.obs.list<-as.list(one.obs$Site.Trait.Sp)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait.Sp%notin%one.obs.list)

#Subset trait
all_traits<-NULL
for (t in unique(try.ttt$TraitShort)){
  trait<-subset(try.ttt,TraitShort==t)
  
  #pick unique site
  
  distance_frame<-trait
  all_sites<-NULL
  
  for(i in unique(trait$Local)){
    distance_frame<-trait
    #Identify coordinates
    pt<-as.matrix(subset(distance_frame,Local==i)[1,c(4,5)])
    #Matrix of all coordinates
    pts<-as.matrix(distance_frame[,c(4,5)])
    #Calculate distance
    distance_frame$distance<-round(spDistsN1(pts, pt, longlat = TRUE),5)
    
    #create object of lenght of unique distances
    length<-length(unique(distance_frame$distance))
    
    Scale_table<-NULL
    Out<-NULL
    repeat {
      #Add next site
      subset<-subset(distance_frame,distance==min(distance))
      distance_frame<-subset(distance_frame,distance!=min(distance))
      Scale_table<-rbind(Scale_table,subset)
      
      #Sample no. values in site to remove additional data problem
      no_samples<-nrow(pt)
      #Number of species in analysis
      no_species<-length(unique(Scale_table$AccSpeciesName))
      #Random sample from Scale_table no. values of initial site
      #Scale_table[sample(nrow(Scale_table), no_samples), ]
      
      
      #Conduct partition
      mm.scale<-tryCatch(as.data.frame(varcomp(lme(log10(StdValue) ~1, random= ~1|AccSpeciesName, data=Scale_table,na.action=na.omit),1)[1:2]), error=function(e) data.frame(c(NA,NA)))
      names(mm.scale)<-c("Value")
      cv<-as.data.frame(sd(Scale_table$StdValue)/mean(Scale_table$StdValue))
      names(cv)<-c("CV")
      mm.scale$scale<-max(Scale_table$distance)
      mm.scale$ndata<-nrow(Scale_table)
      mm.scale<- cbind(Taxonomic_hierarchy = rownames(mm.scale), cv, mm.scale,i,no_species)
      Out<-rbind(Out,mm.scale)
      if (length(unique(Scale_table$distance)) == length){
        break
      }
    }
    site_scale<-Out
    all_sites<-rbind(all_sites,site_scale)
  }
  all_sites$trait<-t
  all_traits<-rbind(all_traits,all_sites)
}

save(all_traits,file="...")

####Arctic Only - Temp####
#REMOVE NON-ARCTIC DATA - Low temp
try.ttt<-try.ttt2
try.ttt<-subset(try.ttt,!is.na(Lat)&!is.na(Lon))

try.ttt$Lat<-as.numeric(as.character(try.ttt$Lat))
try.ttt$Lon<-as.numeric(as.character(try.ttt$Lon))

#Remove NAs
try.ttt<-subset(try.ttt,!is.na(Lat)&!is.na(Lon))

CHELSA_temp<- raster("Z:/Climate_Data/Chelsa/CHELSA_bio10_1.tif")

try.ttt$temp<-raster::extract(CHELSA_temp, try.ttt[,c(5,4)],df = TRUE)[,2] ## Temperature data are in °C × 10 to reduce the file sizes
try.ttt<-subset(try.ttt,temp<=0)

#Remove sites that only <3 species 
one.species<- try.ttt %>% 
  group_by(TraitShort,Local) %>% 
  summarise(n.sp = length(unique(AccSpeciesName))) %>%
  filter(n.sp<3)

try.ttt$Site.Trait<-paste(try.ttt$TraitShort,try.ttt$Local)
one.species$Site.Trait<-paste(one.species$TraitShort,one.species$Local)
one.species.list<-as.list(one.species$Site.Trait)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait%notin%one.species.list)

#& species <5obs/site
one.obs<- try.ttt %>% 
  group_by(TraitShort,Local,AccSpeciesName) %>% 
  summarise(n.obs = length(AccSpeciesName)) %>%
  filter(n.obs<3)

try.ttt$Site.Trait.Sp<-paste(try.ttt$TraitShort,try.ttt$Local,try.ttt$AccSpeciesName)
one.obs$Site.Trait.Sp<-paste(one.obs$TraitShort,one.obs$Local,one.obs$AccSpeciesName)
one.obs.list<-as.list(one.obs$Site.Trait.Sp)
try.ttt<-subset(try.ttt,try.ttt$Site.Trait.Sp%notin%one.obs.list)

#Subset trait
all_traits<-NULL
for (t in unique(try.ttt$TraitShort)){
  trait<-subset(try.ttt,TraitShort==t)
  
  #pick unique site
  
  distance_frame<-trait
  all_sites<-NULL
  
  for(i in unique(trait$Local)){
    distance_frame<-trait
    #Identify coordinates
    pt<-as.matrix(subset(distance_frame,Local==i)[1,c(4,5)])
    #Matrix of all coordinates
    pts<-as.matrix(distance_frame[,c(4,5)])
    #Calculate distance
    distance_frame$distance<-round(spDistsN1(pts, pt, longlat = TRUE),5)
    
    #create object of lenght of unique distances
    length<-length(unique(distance_frame$distance))
    
    Scale_table<-NULL
    Out<-NULL
    repeat {
      #Add next site
      subset<-subset(distance_frame,distance==min(distance))
      distance_frame<-subset(distance_frame,distance!=min(distance))
      Scale_table<-rbind(Scale_table,subset)
      
      #Sample no. values in site to remove additional data problem
      no_samples<-nrow(pt)
      #Number of species in analysis
      no_species<-length(unique(Scale_table$AccSpeciesName))
      #Random sample from Scale_table no. values of initial site
      #Scale_table[sample(nrow(Scale_table), no_samples), ]
      
      
      #Conduct partition
      mm.scale<-tryCatch(as.data.frame(varcomp(lme(log10(StdValue) ~1, random= ~1|AccSpeciesName, data=Scale_table,na.action=na.omit),1)[1:2]), error=function(e) data.frame(c(NA,NA)))
      names(mm.scale)<-c("Value")
      cv<-as.data.frame(sd(Scale_table$StdValue)/mean(Scale_table$StdValue))
      names(cv)<-c("CV")
      mm.scale$scale<-max(Scale_table$distance)
      mm.scale$ndata<-nrow(Scale_table)
      mm.scale<- cbind(Taxonomic_hierarchy = rownames(mm.scale), cv, mm.scale,i,no_species)
      Out<-rbind(Out,mm.scale)
      if (length(unique(Scale_table$distance)) == length){
        break
      }
    }
    site_scale<-Out
    all_sites<-rbind(all_sites,site_scale)
  }
  all_sites$trait<-t
  all_traits<-rbind(all_traits,all_sites)
}

save(all_traits,file="...")


#Load in processed data####


####With species only#####

#####Breakpoint analysis####

# -------------------
# analyse breakpoints
# -------------------

#Load data
load("...")
max(all_traits$scale)

#Remove runs with fewer than 3 species (should be removed but for some reasons a few aren't)
all_traits<-subset(all_traits,no_species>2)
all_traits<-subset(all_traits,ndata>2)

#Remove odd extra tax hierarchies
all_traits<-subset(all_traits,Taxonomic_hierarchy=="Within"|Taxonomic_hierarchy=="AccSpeciesName")


#Binning - scale
all_traits$bin<- cut(all_traits$scale, breaks=c(seq(0,20000,10)))

bin0<-subset(all_traits,is.na(bin))
bin0$bin<-"(1,1)"
bin1<-subset(all_traits,!is.na(bin))
all_traits<-rbind(bin0,bin1)

##Add median bins - scale##
medians<-ddply(all_traits,.(trait,Taxonomic_hierarchy,bin),summarise,
               median = median(Value,na.rm=TRUE),
               UP = quantile(Value,0.975,na.rm=TRUE),
               LP = quantile(Value,0.025,na.rm=TRUE))

medians$bin.nos<-gsub('^.|.$', '', medians$bin)
medians$bin.nos.a<-as.numeric(do.call( rbind, strsplit(medians$bin.nos, ","))[,1])
medians$bin.nos.b<-as.numeric(do.call( rbind, strsplit(medians$bin.nos, ","))[,2])
medians$bin<-medians$bin.nos.a+((medians$bin.nos.b-medians$bin.nos.a)/2)

#Binning - species
max_sp<-max(all_traits$no_species)
all_traits$sp_bin<- cut(all_traits$no_species, breaks=c(seq(0,max_sp,1)))

##Add median bins - species##
medians_sp<-ddply(all_traits,.(trait,Taxonomic_hierarchy,sp_bin),summarise,
                  median = median(Value,na.rm=TRUE),
                  UP = quantile(Value,0.975,na.rm=TRUE),
                  LP = quantile(Value,0.025,na.rm=TRUE))


medians_sp$bin.nos<-gsub('^.|.$', '', medians_sp$sp_bin)
medians_sp$bin.nos.a<-as.numeric(do.call( rbind, strsplit(medians_sp$bin.nos, ","))[,1])
medians_sp$bin.nos.b<-as.numeric(do.call( rbind, strsplit(medians_sp$bin.nos, ","))[,2])
medians_sp$bin<-medians_sp$bin.nos.a+((medians_sp$bin.nos.b-medians_sp$bin.nos.a)/2)

#Binning - CV
all_traits$bin<- cut(all_traits$scale, breaks=c(seq(0,20000,50)))

bin0<-subset(all_traits,is.na(bin))
bin0$bin<-"(1,1)"
bin1<-subset(all_traits,!is.na(bin))
all_traits<-rbind(bin0,bin1)

##Add median bins - cv##
medians_cv<-ddply(all_traits,.(trait,Taxonomic_hierarchy,bin),summarise,
               median = median(CV,na.rm=TRUE),
               UP = quantile(CV,0.975,na.rm=TRUE),
               LP = quantile(CV,0.025,na.rm=TRUE))

medians_cv$bin.nos<-gsub('^.|.$', '', medians_cv$bin)
medians_cv$bin.nos.a<-as.numeric(do.call( rbind, strsplit(medians_cv$bin.nos, ","))[,1])
medians_cv$bin.nos.b<-as.numeric(do.call( rbind, strsplit(medians_cv$bin.nos, ","))[,2])
medians_cv$bin<-medians_cv$bin.nos.a+((medians_cv$bin.nos.b-medians_cv$bin.nos.a)/2)

#Remove zeros so logs to zero 
all_traits$scale<-all_traits$scale+1

#Log scale
all_traits$scale<-log10(all_traits$scale)
nrow(subset(all_traits,is.na(scale)))==0#check for nans

#Log species
all_traits$no_species<-log10(all_traits$no_species)

#Add significances
all_traits$log_bin<- cut(all_traits$scale, breaks=c(seq(0,max(all_traits$scale),
                                                        length.out = 10)))

bin0<-subset(all_traits,is.na(log_bin))
bin0$log_bin<-"(0,0)"
bin1<-subset(all_traits,!is.na(log_bin))
all_traits<-rbind(bin0,bin1)

a<-as.data.frame(as.character(do.call( rbind, strsplit(all_traits$log_bin, ","))[,1]))
names(a)<-"bin"
lower<-as.numeric(as.character(str_sub(a$bin, 2)))

a<-as.data.frame(as.character(do.call( rbind, strsplit(all_traits$log_bin, ","))[,2]))
names(a)<-"bin"
upper<-as.numeric(as.character(substr(a$bin,1,nchar(as.vector(a$bin))-1)))

all_traits$log_bin_centroid<-lower+(upper-lower)/2
all_traits$upper<-upper
all_traits$lower<-lower
  

medians_sp$bin.nos.b<-as.numeric(do.call( rbind, strsplit(medians_sp$bin.nos, ","))[,2])

all_traits$log_sp_bin<- cut(all_traits$no_species, breaks=c(seq(0,max(all_traits$no_species),
                                                        length.out = 10)))

a<-as.data.frame(as.character(do.call( rbind, strsplit(as.character(all_traits$log_sp_bin), ","))[,1]))
names(a)<-"bin"
lower<-as.numeric(as.character(str_sub(a$bin, 2)))

a<-as.data.frame(as.character(do.call( rbind, strsplit(as.character(all_traits$log_sp_bin), ","))[,2]))
names(a)<-"bin"
upper<-as.numeric(as.character(substr(a$bin,1,nchar(as.vector(a$bin))-1)))

all_traits$log_sp_bin_centroid<-(upper-lower)/2
all_traits$upper_sp<-upper
all_traits$lower_sp<-lower


library(broom)


t_tests<-all_traits %>%
  group_by(trait,log_bin,upper,lower,log_bin_centroid) %>%
  do(model = tidy(lm(.$Value~.$Taxonomic_hierarchy))$p.value[2])

t_tests$significance<-ifelse(t_tests$model >0.05,"*","")
t_tests$significance<-ifelse(is.na(t_tests$model),"",t_tests$significance)

t_tests_sp<-all_traits %>%
  group_by(trait,log_sp_bin,upper_sp,lower_sp,log_sp_bin_centroid) %>%
  do(model = tidy(lm(.$Value~.$Taxonomic_hierarchy))$p.value[2])

t_tests_sp$significance<-ifelse(t_tests_sp$model >0.05,"*","")
t_tests_sp$significance<-ifelse(is.na(t_tests_sp$model),"",t_tests_sp$significance)

##Using one break point######

scales_one<-NULL#set blank object
all_breaks_one<-NULL
all_mod_breaks_one<-NULL

#Subset to trait####
for(t in unique(all_traits$trait)){
  all_traits_trait<-subset(all_traits,trait==t)
  all_traits_trait<-all_traits_trait[order(all_traits_trait$scale),]
  #Remove odd extra tax hierarchies
  all_traits_trait<-subset(all_traits_trait,Taxonomic_hierarchy=="Within"|Taxonomic_hierarchy=="AccSpeciesName")
  
  ##RUN FOR SCALE####
  
  #Sample random subset based on bins to have equal sample sizes
  my.lm <- lm(Value ~ scale, data = all_traits_trait[all_traits_trait$Taxonomic_hierarchy=="Within",])
  my.coef<-coef(my.lm)
  
  x=0
  break_points<-NULL
  repeat{
    x<-x+1
    repeat {
      while(TRUE){
        my.seg <- try(segmented(my.lm,
                                psi = 1,
                                seg.Z = ~ scale), 
                      silent=TRUE)
        if(!is(my.seg, 'try-error')) break
      }
      breaks<-my.seg$psi[,2]
      # exit if the condition is met
      if (ifelse(t == "LDMC",breaks[1]<3,breaks[1] <3)) break
    }
    breaks_first<-as.data.frame(cbind(run = x,breaks_1 = breaks[1]))
    breaks_run<-breaks_first
    break_points<-as.data.frame(rbind(break_points,breaks_run))
    if (x == 5) break
  }
  
  LQ_1<-as.numeric(quantile(break_points$breaks_1,0.025))
  UQ_1<-as.numeric(quantile(break_points$breaks_1,0.975))
  mean_1<-as.numeric(mean(break_points$breaks_1))
  median_1<-as.numeric(median(break_points$breaks_1))
  
  trait_breaks<-as.data.frame(cbind(t,LQ_1,UQ_1,mean_1,median_1))
  all_breaks_one<-rbind(all_breaks_one,trait_breaks)
  
  
  #Rerun using new bounds
  repeat {
    while(TRUE){
      my.seg <- try(segmented(my.lm,
                              psi = 1,
                              seg.Z = ~ scale),
                    silent=TRUE)
      if(!is(my.seg, 'try-error')) break
    }
    breaks<-my.seg$psi[,2]
    # exit if the condition is met
    if (round(breaks[1],1) == round(median_1,1)) break
  }
  
  # display the summary
  summary(my.seg)
  
  # get the fitted data
  my.fitted <- fitted(my.seg)
  my.model_within <- data.frame(scale = all_traits_trait[all_traits_trait$Taxonomic_hierarchy=="Within",]$scale, Value = my.fitted)
  
  #Make sure break points are plotted
  break_preds<-as.data.frame(cbind(breaks[1],coef(my.seg)[1] + coef(my.seg)[2]*breaks[1]))
  
  names(break_preds)<-c("scale","Value")
  
  my.model_within<-rbind(my.model_within,break_preds)

  
  # Repeat for among (1-within)
  my.model_among <- data.frame(scale = all_traits_trait[all_traits_trait$Taxonomic_hierarchy=="Within",]$scale, Value = 1-my.fitted)
  
  #Make sure break points are plotted
  break_preds2<-as.data.frame(cbind(breaks[1],1-break_preds[1,2]))
  names(break_preds2)<-c("scale","Value")
  
  my.model_among<-rbind(my.model_among,break_preds2)
  
  mod_breaks<-as.data.frame(cbind(t,breaks[1],breaks[2]))
  all_mod_breaks_one<-rbind(all_mod_breaks_one,mod_breaks)

  ggplot(my.model_within, aes(x = scale, y = Value)) + geom_line()
  ggplot(my.model_among, aes(x = scale, y = Value)) + geom_line()
  
  
  #Get median trait
  median_trait<-subset(medians,trait==t)
  
  g1_species <- ggplot(median_trait[median_trait$Taxonomic_hierarchy=="AccSpeciesName",]) + 
    stat_smooth(aes(x = log10(bin), y = LP), method = "lm", formula = y ~ poly(x, degree = 2), se = FALSE) +
    stat_smooth(aes(x = log10(bin), y = UP), method = "lm", formula = y ~ poly(x, degree = 2), se = FALSE)
  g1_species
  
  g1_within <- ggplot(median_trait[median_trait$Taxonomic_hierarchy=="Within",]) + 
    stat_smooth(aes(x = log10(bin), y = LP), method = "lm", formula = y ~ poly(x, degree = 2), se = FALSE) +
    stat_smooth(aes(x = log10(bin), y = UP), method = "lm", formula = y ~ poly(x, degree = 2), se = FALSE)
  g1_within
  
  # build plot object for rendering 
  gg1_species <- ggplot_build(g1_species)
  gg1_within <- ggplot_build(g1_within)
  
  # extract data for the loess lines from the 'data' slot
  df2_species <- data.frame(x = gg1_species$data[[1]]$x,
                            ymin = gg1_species$data[[1]]$y,
                            ymax = gg1_species$data[[2]]$y) 
  
  df2_within <- data.frame(x = gg1_within$data[[1]]$x,
                           ymin = gg1_within$data[[1]]$y,
                           ymax = gg1_within$data[[2]]$y) 
  
  #Add significances
  t_test_trait<-subset(t_tests,trait == t)
  t_test_trait[t_test_trait$significance != "*",]$upper<-NA
  sig_upper<-tryCatch(max(t_test_trait$upper,na.rm = T), error=function(e) as.numeric(NA))
  sig_upper<-ifelse(abs(sig_upper)==Inf,-10,sig_upper)
  sig_lower<-tryCatch(min(t_test_trait$upper,na.rm = T), error=function(e) as.numeric(NA))
  sig_lower<-ifelse(abs(sig_lower)==Inf,-11,sig_lower)
  sig_lower<-ifelse(sig_lower==sig_upper,min(df2_species$x),sig_lower)
  
  (scale<-ggplot(all_traits_trait)+
      annotate("rect", xmin = sig_lower,
               xmax = sig_upper,
               ymin = -0.05, ymax = 1.05,fill = "darkorchid", alpha = 0.05, colour = "darkorchid",linetype = "dotted")+
      annotate("rect", xmin = as.numeric(as.character(all_breaks_one[all_breaks_one$t==t,2])),
               xmax = as.numeric(as.character(all_breaks_one[all_breaks_one$t==t,4])),
               ymin = -0.05, ymax = 1.05,alpha=0.1)+
      #annotate("text",label = "NS", x = ((t_test_trait$upper - 0)/2), y = 1.1, size = 3.5)+
      geom_point(aes(scale,Value,colour=Taxonomic_hierarchy),alpha=0.2,shape=16,size=1.6)+
      geom_ribbon(data = df2_species, aes(x = x, ymin = ymin, ymax = ymax),fill="#F8766D",alpha = 0.4)+
      geom_ribbon(data = df2_within, aes(x = x, ymin = ymin, ymax = ymax),fill="#619CFF",alpha = 0.4)+
      #geom_smooth(aes(scale,Value,colour=Taxonomic_hierarchy),
      #method = "lm",
      #formula = y ~ poly(x, degree = 3), 
      #se = FALSE)+
      scale_x_continuous(breaks = 1:4,
                         labels = c(expression(10^{1}), expression(10^{2}), expression(10^{3}), expression(10^{4})))+
      scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))+
      coord_cartesian(ylim=c(-0.05, 1.15),xlim = c(min(df2_species$x),max(df2_species$x)))+
      geom_line(data = my.model_within, aes(x = scale, y = Value),colour = "blue",size=1)+
      geom_line(data = my.model_among, aes(x = scale, y = Value),colour = "red",size=1)+
      labs(x=expression(paste("Spatial Extent (", km^2,")")),y="Proportion of variance explained")+
      ggtitle(t)+
      scale_colour_manual(values = c("#F8766D","#619CFF","Red","Blue"),
                          breaks=c("AccSpeciesName", "Within"),
                          labels=c("Among species","Within species"),
                          name = "Source of Variation")+
      theme_bw()+
      theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))
  
  scales_one[[t]]<-scale
}

pdf(file = "scripts/users/hthomas/Output_Images/Traits_1/scale_figs_onebreak.pdf", height = 6,width = 10)
grid.arrange(arrangeGrob(scales_one[["LeafArea"]]+
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="LeafArea",2])),xend = as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="LeafArea",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         scales_one[["PlantHeight"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="PlantHeight",2])),xend = as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="PlantHeight",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         scales_one[["LMA"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="LMA",2])),xend = as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="LMA",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         scales_one[["LeafN"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="LeafN",2])),xend = as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="LeafN",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         scales_one[["LDMC"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="LDMC",2])),xend = as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="LDMC",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         scales_one[["SeedMass"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="SeedMass",2])),xend = as.numeric(as.character(all_mod_breaks_one[all_mod_breaks_one$t=="SeedMass",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         nrow=2))
dev.off()


##RUN FOR SPECIES####
all_mod_breaks_sp<-NULL
all_breaks_sp<-NULL
no_species_full<-NULL

#Subset to trait####
for(t in unique(all_traits$trait)){
  t = "LDMC"
  all_traits_trait<-subset(all_traits,trait==t)
  #Remove odd extra tax hierarchies
  all_traits_trait<-subset(all_traits_trait,Taxonomic_hierarchy=="Within"|Taxonomic_hierarchy=="AccSpeciesName")
  
    ##RUN FOR no_species####
  my.lm <- lm(Value ~ no_species, data = all_traits_trait[all_traits_trait$Taxonomic_hierarchy=="Within",])
  my.coef<-coef(my.lm)
  
  x=0
  break_points<-NULL
  repeat{
    x<-x+1
    repeat {
      while(TRUE){
        my.seg <- try(segmented(my.lm,
                                seg.Z = ~ no_species,
                                psi = c(1,1.66)), silent=TRUE)
        if(!is(my.seg, 'try-error')) break
      }
      breaks<-my.seg$psi[,2]
      # exit if the condition is met
      if (breaks[2] > 1 & breaks[1] < 2 & (breaks[2] - breaks[1]>0.5)) break
    }
    breaks_first<-as.data.frame(cbind(run = x,breaks_1 = breaks[1]))
    breaks_second<-as.data.frame(cbind(run = x, breaks_2 = breaks[2]))
    breaks_run<-merge(breaks_first,breaks_second)
    break_points<-as.data.frame(rbind(break_points,breaks_run))
    if (x == 5) break
  }
  
  LQ_1<-as.numeric(quantile(break_points$breaks_1,0.1))
  UQ_1<-as.numeric(quantile(break_points$breaks_1,0.9))
  median_1<-as.numeric(median(break_points$breaks_1))
  mean_1<-as.numeric(mean(break_points$breaks_1))
  
  LQ_2<-as.numeric(quantile(break_points$breaks_2,0.1))
  UQ_2<-as.numeric(quantile(break_points$breaks_2,0.9))
  median_2<-as.numeric(median(break_points$breaks_2))
  mean_2<-as.numeric(mean(break_points$breaks_2))
  
  trait_breaks<-as.data.frame(cbind(t,LQ_1,LQ_2,UQ_1,UQ_2,median_1,median_2,mean_1,mean_2))
  all_breaks_sp<-rbind(all_breaks_sp,trait_breaks)
  
  #Rerun using new bounds
  repeat{
    repeat {
      while(TRUE){
        my.seg <- try(segmented(my.lm, 
                                seg.Z = ~ no_species, 
                                psi = c(1,1.66)), silent=TRUE)
        if(!is(my.seg, 'try-error')) break
      }
      breaks<-my.seg$psi[,2]
      # exit if the condition is met
      if (round(breaks[2],1) == round(median_2,1)) break
    }
    if (round(breaks[1],1) == round(median_1,1)) break
  }
  
  # display the summary
  summary(my.seg)
  
  # get the fitted data
  my.fitted <- fitted(my.seg)
  my.model_within <- data.frame(no_species = all_traits_trait[all_traits_trait$Taxonomic_hierarchy=="Within",]$no_species, Value = my.fitted)
  
  #Make sure break points are plotted
  break_preds<-as.data.frame(rbind(cbind(breaks[1],coef(my.seg)[1] + coef(my.seg)[2]*breaks[1]),
                                   (cbind(breaks[2],coef(my.seg)[1] + coef(my.seg)[2]*breaks[1] + (coef(my.seg)[3]+coef(my.seg)[2])*(breaks[2]-breaks[1])))))
  names(break_preds)<-c("no_species","Value")
  
  my.model_within<-rbind(my.model_within,break_preds)
  
  
  # Repeat for among (1-within)
  my.model_among <- data.frame(no_species = all_traits_trait[all_traits_trait$Taxonomic_hierarchy=="Within",]$no_species, Value = 1-my.fitted)
  
  #Make sure break points are plotted
  break_preds2<-as.data.frame(rbind(cbind(breaks[1],1-break_preds[1,2]),
                                    (cbind(breaks[2],1-break_preds[2,2]))))
  names(break_preds2)<-c("no_species","Value")
  
  my.model_among<-rbind(my.model_among,break_preds2)
  
  mod_breaks<-as.data.frame(cbind(t,breaks[1],breaks[2]))
  all_mod_breaks_sp<-rbind(all_mod_breaks_sp,mod_breaks)

  # plot the fitted models
  ggplot(my.model_within, aes(x = no_species, y = Value)) + geom_line()
  ggplot(my.model_among, aes(x = no_species, y = Value)) + geom_line()
  
  #Get median trait
  median_trait<-subset(medians_sp,trait==t)
  
  g1_species <- ggplot(median_trait[median_trait$Taxonomic_hierarchy=="AccSpeciesName",]) +
    stat_smooth(aes(x = log10(bin), y = LP), method = "loess", se = FALSE) +
    stat_smooth(aes(x = log10(bin), y = UP), method = "loess", se = FALSE)
  g1_species
  
  g1_within <- ggplot(median_trait[median_trait$Taxonomic_hierarchy=="Within",]) +
    stat_smooth(aes(x = log10(bin), y = LP), method = "loess", se = FALSE) +
    stat_smooth(aes(x = log10(bin), y = UP), method = "loess", se = FALSE)
  g1_within
  
  #build plot object for rendering
  gg1_species <- ggplot_build(g1_species)
  gg1_within <- ggplot_build(g1_within)
  
  #extract data for the loess lines from the 'data' slot
  df2_species <- data.frame(x = gg1_species$data[[1]]$x,
                            ymin = gg1_species$data[[1]]$y,
                            ymax = gg1_species$data[[2]]$y)
  
  df2_within <- data.frame(x = gg1_within$data[[1]]$x,
                           ymin = gg1_within$data[[1]]$y,
                           ymax = gg1_within$data[[2]]$y)
  
  
  #Add significances
  t_test_sp_trait<-subset(t_tests_sp,trait == t)
  t_test_sp_trait[t_test_sp_trait$significance != "*",]$upper_sp<-NA
  sig_upper_sp<-tryCatch(max(t_test_sp_trait$upper_sp,na.rm = T), error=function(e) as.numeric(NA))
  sig_upper_sp<-ifelse(abs(sig_upper_sp)==Inf,-10,sig_upper_sp)
  sig_lower_sp<-tryCatch(min(t_test_sp_trait$upper_sp,na.rm = T), error=function(e) as.numeric(NA))
  sig_lower_sp<-ifelse(abs(sig_lower_sp)==Inf,-11,sig_lower_sp)
  sig_lower_sp<-ifelse(sig_lower_sp==sig_upper_sp,min(df2_species$x),sig_lower_sp)
  
  
  (no_species_fig<-ggplot(all_traits_trait)+
      annotate("rect", xmin = sig_lower_sp,
               xmax = sig_upper_sp,
               ymin = -0.05, ymax = 1.05,fill = "darkorchid", alpha = 0.05, colour = "darkorchid",linetype = "dotted")+
      annotate("rect", xmin = as.numeric(as.character(all_breaks_sp[all_breaks_sp$t==t,2])),
               xmax = as.numeric(as.character(all_breaks_sp[all_breaks_sp$t==t,4])),
               ymin = -0.05, ymax = 1.05,alpha=0.1)+
      annotate("rect", xmin = as.numeric(as.character(all_breaks_sp[all_breaks_sp$t==t,3])),
               xmax = as.numeric(as.character(all_breaks_sp[all_breaks_sp$t==t,5])),
               ymin = -0.05, ymax = 1.05,alpha=0.1)+
      #annotate("text",label = "NS", x = ((t_test_sp_trait$upper_sp - (t_test_sp_trait$upper_sp - min(df2_species$x))/2)), y = 1.1, size = 3.5)+
    geom_point(aes(no_species,Value,colour=Taxonomic_hierarchy),alpha=0.2,shape=16,size=1.6)+
    geom_ribbon(data = df2_species, aes(x = x, ymin = ymin, ymax = ymax),fill="#F8766D",alpha = 0.4)+
    geom_ribbon(data = df2_within, aes(x = x, ymin = ymin, ymax = ymax),fill="#619CFF",alpha = 0.4)+
    #geom_smooth(aes(no_species,Value,colour=Taxonomic_hierarchy), method = "lm", formula = y ~ poly(x, degree = 3), se = TRUE)+
    #method = "lm",
    #formula = y ~ poly(x, degree = 3), 
    #se = FALSE)+
      coord_cartesian(ylim=c(-0.05, 1.15))+
    geom_line(data = my.model_within, aes(x = no_species, y = Value),colour = "blue",size=1)+
    geom_line(data = my.model_among, aes(x = no_species, y = Value),colour = "red",size=1)+
    labs(x="Number of species",y="Variance Explained")+
    ggtitle(t)+
    scale_colour_manual(values = c("#F8766D","#619CFF","Red","Blue"),
                        breaks=c("AccSpeciesName", "Within"),
                        labels=c("Among species","Within species"),
                        name = "Source of Variation")+
    scale_x_continuous(breaks = c(1,2),
                       labels = c(10,100))+
      scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))+
      coord_cartesian(ylim=c(-0.05, 1.15),xlim = c(min(df2_species$x),max(df2_species$x)))+
    theme_bw()+
    theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))
  
  no_species_full[[t]]<-no_species_fig
}


#Extract legend
mylegend<-g_legend(no_species_fig)

pdf(file = "scripts/users/hthomas/Output_Images/Traits_1/species_figs_twobreak.pdf", height = 6,width = 10)
grid.arrange(arrangeGrob(no_species_full[["LeafArea"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafArea",2])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafArea",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafArea",3])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafArea",3])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full[["PlantHeight"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="PlantHeight",2])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="PlantHeight",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="PlantHeight",3])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="PlantHeight",3])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full[["LMA"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LMA",2])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LMA",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LMA",3])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LMA",3])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full[["LeafN"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafN",2])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafN",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafN",3])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafN",3])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full[["LDMC"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LDMC",2])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LDMC",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LDMC",3])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LDMC",3])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full[["SeedMass"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="SeedMass",2])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="SeedMass",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="SeedMass",3])),xend = as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="SeedMass",3])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         nrow=2))
dev.off()

#One Break species#########################

all_mod_breaks_sp_one<-NULL
all_breaks_sp_one<-NULL
no_species_full_one<-NULL

#Subset to trait####
for(t in unique(all_traits$trait)){
  all_traits_trait<-subset(all_traits,trait==t)
  #Remove odd extra tax hierarchies
  all_traits_trait<-subset(all_traits_trait,Taxonomic_hierarchy=="Within"|Taxonomic_hierarchy=="AccSpeciesName")
  
  ##RUN FOR no_species####
  my.lm <- lm(Value ~ no_species, data = all_traits_trait[all_traits_trait$Taxonomic_hierarchy=="Within",])
  my.coef<-coef(my.lm)
  
  x=0
  break_points<-NULL
  repeat{
    x<-x+1
    repeat {
      while(TRUE){
        my.seg <- try(segmented(my.lm,
                                seg.Z = ~ no_species),
                                silent=TRUE)
        if(!is(my.seg, 'try-error')) break
      }
      breaks<-my.seg$psi[,2]
      # exit if the condition is met
      if (breaks[1] < 2.5) break
    }
    breaks_first<-as.data.frame(cbind(run = x,breaks_1 = breaks[1]))
    breaks_run<-breaks_first
    break_points<-as.data.frame(rbind(break_points,breaks_run))
    if (x == 5) break
  }
  
  LQ_1<-as.numeric(quantile(break_points$breaks_1,0.1))
  UQ_1<-as.numeric(quantile(break_points$breaks_1,0.9))
  median_1<-as.numeric(median(break_points$breaks_1))
  mean_1<-as.numeric(mean(break_points$breaks_1))
  
  trait_breaks<-as.data.frame(cbind(t,LQ_1,UQ_1,median_1,mean_1))
  all_breaks_sp_one<-rbind(all_breaks_sp_one,trait_breaks)
  
  #Rerun using new bounds

    repeat {
      while(TRUE){
        my.seg <- try(segmented(my.lm, 
                                seg.Z = ~ no_species), 
                                silent=TRUE)
        if(!is(my.seg, 'try-error')) break
      }
      breaks<-my.seg$psi[,2]
      # exit if the condition is met
    if (round(breaks[1],1) == round(median_1,1)) break
  }
  
  # display the summary
  summary(my.seg)
  
  # get the fitted data
  my.fitted <- fitted(my.seg)
  my.model_within <- data.frame(no_species = all_traits_trait[all_traits_trait$Taxonomic_hierarchy=="Within",]$no_species, Value = my.fitted)
  
  #Make sure break points are plotted
  break_preds<-as.data.frame(cbind(breaks[1],coef(my.seg)[1] + coef(my.seg)[2]*breaks[1]))
  
  names(break_preds)<-c("no_species","Value")
  
  my.model_within<-rbind(my.model_within,break_preds)
  
  
  # Repeat for among (1-within)
  my.model_among <- data.frame(no_species = all_traits_trait[all_traits_trait$Taxonomic_hierarchy=="Within",]$no_species, Value = 1-my.fitted)
  
  #Make sure break points are plotted
  break_preds2<-as.data.frame(cbind(breaks[1],1-break_preds[1,2]))
  names(break_preds2)<-c("no_species","Value")
  
  my.model_among<-rbind(my.model_among,break_preds2)
  
  mod_breaks<-as.data.frame(cbind(t,breaks[1],breaks[2]))
  all_mod_breaks_sp_one<-rbind(all_mod_breaks_sp_one,mod_breaks)
  
  # plot the fitted models
  ggplot(my.model_within, aes(x = no_species, y = Value)) + geom_line()
  ggplot(my.model_among, aes(x = no_species, y = Value)) + geom_line()
  
  #Get median trait
  median_trait<-subset(medians_sp,trait==t)
  
  g1_species <- ggplot(median_trait[median_trait$Taxonomic_hierarchy=="AccSpeciesName",]) +
    stat_smooth(aes(x = log10(bin), y = LP), method = "loess", se = FALSE) +
    stat_smooth(aes(x = log10(bin), y = UP), method = "loess", se = FALSE)
  g1_species
  
  g1_within <- ggplot(median_trait[median_trait$Taxonomic_hierarchy=="Within",]) +
    stat_smooth(aes(x = log10(bin), y = LP), method = "loess", se = FALSE) +
    stat_smooth(aes(x = log10(bin), y = UP), method = "loess", se = FALSE)
  g1_within
  
  #build plot object for rendering
  gg1_species <- ggplot_build(g1_species)
  gg1_within <- ggplot_build(g1_within)
  
  #extract data for the loess lines from the 'data' slot
  df2_species <- data.frame(x = gg1_species$data[[1]]$x,
                            ymin = gg1_species$data[[1]]$y,
                            ymax = gg1_species$data[[2]]$y)
  
  df2_within <- data.frame(x = gg1_within$data[[1]]$x,
                           ymin = gg1_within$data[[1]]$y,
                           ymax = gg1_within$data[[2]]$y)
  
  #Add significances
  t_test_sp_trait<-subset(t_tests_sp,trait == t)
  t_test_sp_trait[t_test_sp_trait$significance != "*",]$upper_sp<-NA
  sig_upper_sp<-tryCatch(max(t_test_sp_trait$upper_sp,na.rm = T), error=function(e) as.numeric(NA))
  sig_upper_sp<-ifelse(abs(sig_upper_sp)==Inf,-10,sig_upper_sp)
  sig_lower_sp<-tryCatch(min(t_test_sp_trait$upper_sp,na.rm = T), error=function(e) as.numeric(NA))
  sig_lower_sp<-ifelse(abs(sig_lower_sp)==Inf,-11,sig_lower_sp)
  sig_lower_sp<-ifelse(sig_lower_sp==sig_upper_sp,min(df2_species$x),sig_lower_sp)
  
  (no_species_fig<-ggplot(all_traits_trait)+
      annotate("rect", xmin = sig_lower_sp,
               xmax = sig_upper_sp,
               ymin = -0.05, ymax = 1.05,fill = "darkorchid", alpha = 0.05, colour = "darkorchid",linetype = "dotted")+
      annotate("rect", xmin = as.numeric(as.character(all_breaks_sp_one[all_breaks_sp_one$t==t,2])),
               xmax = as.numeric(as.character(all_breaks_sp_one[all_breaks_sp_one$t==t,4])),
               ymin = -0.05, ymax = 1.05,alpha=0.1)+
      #annotate("text",label = "NS", x = ((t_test_sp_trait$upper_sp - (t_test_sp_trait$upper_sp - min(df2_species$x))/2)), y = 1.1, size = 3.5)+
      geom_point(aes(no_species,Value,colour=Taxonomic_hierarchy),alpha=0.2,shape=16,size=1.6)+
      geom_ribbon(data = df2_species, aes(x = x, ymin = ymin, ymax = ymax),fill="#F8766D",alpha = 0.4)+
      geom_ribbon(data = df2_within, aes(x = x, ymin = ymin, ymax = ymax),fill="#619CFF",alpha = 0.4)+
      #geom_smooth(aes(no_species,Value,colour=Taxonomic_hierarchy), method = "lm", formula = y ~ poly(x, degree = 3), se = TRUE)+
      #method = "lm",
      #formula = y ~ poly(x, degree = 3), 
      #se = FALSE)+
      coord_cartesian(ylim=c(-0.05, 1.15))+
      geom_line(data = my.model_within, aes(x = no_species, y = Value),colour = "blue",size=1)+
      geom_line(data = my.model_among, aes(x = no_species, y = Value),colour = "red",size=1)+
      labs(x="Number of species",y="Variance Explained")+
      ggtitle(t)+
      scale_colour_manual(values = c("#F8766D","#619CFF","Red","Blue"),
                          breaks=c("AccSpeciesName", "Within"),
                          labels=c("Among species","Within species"),
                          name = "Source of Variation")+
      scale_x_continuous(breaks = c(1,2),
                         labels = c(10,100))+
      scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))+
      coord_cartesian(ylim=c(-0.05, 1.15),xlim = c(min(df2_species$x),max(df2_species$x)))+
      theme_bw()+
      theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))
  
  no_species_full_one[[t]]<-no_species_fig
}


#Extract legend
mylegend<-g_legend(no_species_fig)

pdf(file = "scripts/users/hthomas/Output_Images/Traits_1/species_figs_onebreak.pdf", height = 6,width = 10)
grid.arrange(arrangeGrob(no_species_full_one[["LeafArea"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="LeafArea",2])),xend = as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="LeafArea",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full_one[["PlantHeight"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="PlantHeight",2])),xend = as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="PlantHeight",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full_one[["LMA"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="LMA",2])),xend = as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="LMA",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full_one[["LeafN"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="LeafN",2])),xend = as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="LeafN",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full_one[["LDMC"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="LDMC",2])),xend = as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="LDMC",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         no_species_full_one[["SeedMass"]] + 
                           theme(legend.position="none")+
                           geom_segment(aes(x=as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="SeedMass",2])),xend = as.numeric(as.character(all_mod_breaks_sp_one[all_mod_breaks_sp_one$t=="SeedMass",2])),y = -0.05, yend = 1.05),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         nrow=2))
dev.off()


#Coefficient of variation####################################

all_mod_breaks_cv_one<-NULL
all_breaks_cv_one<-NULL
cv_full_one<-NULL

for(t in unique(all_traits$trait)){
  all_traits_trait<-subset(all_traits,trait==t)
  #Remove odd extra tax hierarchies
  all_traits_trait<-subset(all_traits_trait,Taxonomic_hierarchy=="Within"|Taxonomic_hierarchy=="AccSpeciesName")
  
  ##RUN FOR CV####
  
  #Get median trait
  median_trait<-subset(medians_cv,trait==t)
  
  g1_within <- ggplot(median_trait[median_trait$Taxonomic_hierarchy=="Within",]) +
    stat_smooth(aes(x = log10(bin), y = LP), method = "loess", se = FALSE) +
    stat_smooth(aes(x = log10(bin), y = UP), method = "loess", se = FALSE)
  g1_within
  
  #build plot object for rendering
  gg1_within <- ggplot_build(g1_within)
  
  #extract data for the loess lines from the 'data' slot
  df2_within <- data.frame(x = gg1_within$data[[1]]$x,
                           ymin = gg1_within$data[[1]]$y,
                           ymax = gg1_within$data[[2]]$y)
  
  (cv_fig<-ggplot(all_traits_trait)+
      #annotate("text",label = "NS", x = ((t_test_sp_trait$upper_sp - (t_test_sp_trait$upper_sp - min(df2_species$x))/2)), y = 1.1, size = 3.5)+
      geom_point(aes(scale,CV,colour=Taxonomic_hierarchy),alpha=0.2,shape=16,size=1.6)+
      stat_smooth(aes(scale,CV,colour=Taxonomic_hierarchy),method = "lm", se = F)+
      geom_ribbon(data = df2_within, aes(x = x, ymin = ymin, ymax = ymax),fill="#619CFF",alpha = 0.4)+
      #geom_smooth(aes(no_species,Value,colour=Taxonomic_hierarchy), method = "lm", formula = y ~ poly(x, degree = 3), se = TRUE)+
      #method = "lm",
      #formula = y ~ poly(x, degree = 3), 
      #se = FALSE)+
      labs(x=expression(paste("Spatial Extent (", km^2,")")),y="Coefficient of Variation (CV)")+
      ggtitle(t)+
      scale_x_continuous(breaks = 1:4,
                         labels = c(expression(10^{1}), expression(10^{2}), expression(10^{3}), expression(10^{4})))+
      theme_bw()+
      theme(legend.position = "none")+
      theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))
  
  cv_full_one[[t]]<-cv_fig
}

pdf(file = "scripts/users/hthomas/Output_Images/Traits_1/cv_figs_onebreak.pdf", height = 6,width = 10)
grid.arrange(arrangeGrob(cv_full_one[["LeafArea"]] + 
                           theme(legend.position="none")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         cv_full_one[["PlantHeight"]] + 
                           theme(legend.position="none")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         cv_full_one[["LMA"]] + 
                           theme(legend.position="none")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         cv_full_one[["LeafN"]] + 
                           theme(legend.position="none")+
                          theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         cv_full_one[["LDMC"]] + 
                           theme(legend.position="none")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         cv_full_one[["SeedMass"]] + 
                           theme(legend.position="none")+
                          theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm")),
                         nrow=2))
dev.off()

#---------------------------------------------------------------------------

##OVERALL FIGURE####

#---------------------------------------------------------------------------
pdf(file="scripts/users/hthomas/Output_Images/Traits_1/change_with_scale_log.pdf",width = 8, height = 5)
grid.arrange(arrangeGrob(scales[["PlantHeight"]] + 
                           theme(legend.position="none")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="PlantHeight",2]))),linetype="dashed")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="PlantHeight",3]))),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
                         scales[["LeafArea"]] + 
                           theme(legend.position="none")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="LeafArea",2]))),linetype="dashed")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="LeafArea",3]))),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
                         scales[["LMA"]] + 
                           theme(legend.position="none")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="LMA",2]))),linetype="dashed")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="LMA",3]))),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
                         no_species_full[["PlantHeight"]] + 
                           theme(legend.position="none")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="PlantHeight",2]))),linetype="dashed")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="PlantHeight",3]))),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
                         no_species_full[["LeafArea"]] + 
                           theme(legend.position="none")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafArea",2]))),linetype="dashed")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafArea",3]))),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
                         no_species_full[["LMA"]] + 
                           theme(legend.position="none")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LMA",2]))),linetype="dashed")+
                           geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LMA",3]))),linetype="dashed")+
                           theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
                         ncol=3))
dev.off()


#Supplementary Figure####
pdf(file="scripts/users/hthomas/Output_Images/Traits_1/change_with_scale_log_others.pdf",width = 8, height = 5)
grid.arrange(scales[["LeafN"]] + 
               theme(legend.position="none")+
               geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="LeafN",2]))),linetype="dashed")+
               geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="LeafN",3]))),linetype="dashed")+
               theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
             scales[["LDMC"]] + 
               theme(legend.position="none")+
               geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="LDMC",2]))),linetype="dashed")+
               geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="LDMC",3]))),linetype="dashed")+
               theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
             scales[["SeedMass"]] + 
               theme(legend.position="none")+
               #geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="SeedMass",2]))),linetype="dashed")+
               #geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks[all_mod_breaks$t=="SeedMass",3]))),linetype="dashed")+
               theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
             no_species_full[["LeafN"]] + 
               theme(legend.position="none")+
               geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafN",2]))),linetype="dashed")+
               geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LeafN",3]))),linetype="dashed")+
               theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
             no_species_full[["LDMC"]] + 
               theme(legend.position="none")+
               geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LDMC",2]))),linetype="dashed")+
               geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="LDMC",3]))),linetype="dashed")+
               theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
             no_species_full[["SeedMass"]] + 
               theme(legend.position="none")+
               #geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="SeedMass",2]))),linetype="dashed")+
               #geom_vline(aes(xintercept=as.numeric(as.character(all_mod_breaks_sp[all_mod_breaks_sp$t=="SeedMass",3]))),linetype="dashed")+
               theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")),
             ncol=3)
dev.off()

no_species_full[["PlantHeight"]] + 
  scale_x_log10()

#Sup Fig - balance of below BP and above####
load("scripts/users/hthomas/Traits_1/scale_change_dists_sponly.cleaned.larger_3.rds")

#Remove runs with fewer than 3 species (should be removed but for some reasons a few aren't)
all_traits<-subset(all_traits,no_species>2)
all_traits<-subset(all_traits,ndata>2)

#Remove zeros so logs to zero 
all_traits$scale<-all_traits$scale+1

#Log scale
all_traits$scale<-log10(all_traits$scale)
nrow(subset(all_traits,is.na(scale)))==0#check for nans

all_traits<-subset(all_traits,Taxonomic_hierarchy=="Within"|Taxonomic_hierarchy=="AccSpeciesName")
all_traits<-subset(all_traits,trait!="SeedMass"&trait!="SSD")

balances<-NULL
balances2<-NULL
all_site_balances<-NULL
for (s in unique(all_traits$trait)){ 
  #Subset to trait
  all_traits_within<-subset(all_traits,trait==s & Taxonomic_hierarchy=="Within")
  
  all_mod_breaks$V2<-as.numeric(as.character(all_mod_breaks$V2))
  all_mod_breaks$V3<-as.numeric(as.character(all_mod_breaks$V3))
  
  #Add side of breakpoint
  all_traits_within$BP<-ifelse(all_traits_within$scale<=all_mod_breaks[all_mod_breaks$t==s,]$V2,"Below","Above")
  
  #Calc if >=50%
  all_traits_within$dominant<-ifelse(all_traits_within$Value>=0.5,"Within","Among")
  #all_traits_within$dominant_3rd<-ifelse(all_traits_within$Value>=0.333333333,"Within","Among")
  
  #Build relevant object
  site_balance<- all_traits_within %>%
    group_by(BP) %>%
    summarise(Within = sum(dominant=="Within"), Among = sum(dominant=="Among"))
  
  site_balance <- gather(site_balance, type, length, Within:Among, factor_key=TRUE)
  site_balance$BP<-as.factor(site_balance$BP)
  levels(site_balance$BP)
  
  site_balance$BP = factor(site_balance$BP, levels = c("Below", "Above"))
  
  #Add trait
  site_balance$trait<-s
  
  #Combine
  all_site_balances<-rbind(all_site_balances,site_balance)
  
  site_balance$type <- factor(site_balance$type, levels = c("Among", "Within"))
  
  fig<-ggplot(site_balance,aes(x = BP, y = length,fill = type)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = percent_format())+
    scale_fill_manual(values = c("#F8766D00","#619CFF"),
                      breaks=c("AccSpeciesName", "Within"),
                      labels=c("Among species","Within species"))+
    theme_bw()+
    theme(legend.position="none")+
    labs(x="",y="Proportion of sites")+
    scale_x_discrete(labels=c("Below" = "Below breakpoint", "Above" = "Above breakpoint"))+
    theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    ggtitle(s)
  balances[[s]]<-fig
}

for (s in unique(all_traits$trait)){ 
  #Subset to trait
  all_traits_within<-subset(all_traits,trait==s & Taxonomic_hierarchy=="Within")
  
  all_mod_breaks$V2<-as.numeric(as.character(all_mod_breaks$V2))
  all_mod_breaks$V3<-as.numeric(as.character(all_mod_breaks$V3))
  
  #Add side of breakpoint
  all_traits_within$BP<-ifelse(all_traits_within$scale<=all_mod_breaks[all_mod_breaks$t==s,]$V2,"Below","Above")
  
  #Calc if >=33%
  all_traits_within$dominant<-ifelse(all_traits_within$Value>=0.333333333,"Within","Among")
  
  #Build relevant object
  site_balance<- all_traits_within %>%
    group_by(BP) %>%
    summarise(Within = sum(dominant=="Within"), Among = sum(dominant=="Among"))
  
  site_balance <- gather(site_balance, type, length, Within:Among, factor_key=TRUE)
  site_balance$BP<-as.factor(site_balance$BP)
  levels(site_balance$BP)
  
  site_balance$BP = factor(site_balance$BP, levels = c("Below", "Above"))
  
  #Add trait
  site_balance$trait<-s
  
  #Combine
  all_site_balances<-rbind(all_site_balances,site_balance)
  
  site_balance$type <- factor(site_balance$type, levels = c("Among", "Within"))
  
  fig<-ggplot(site_balance,aes(x = BP, y = length,fill = type)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = percent_format())+
    scale_fill_manual(values = c("#F8766D00","#619CFF"),
                      breaks=c("AccSpeciesName", "Within"),
                      labels=c("Among species","Within species"))+
    theme_bw()+
    theme(legend.position="none")+
    labs(x="",y="Proportion of sites")+
    scale_x_discrete(labels=c("Below" = "Below breakpoint", "Above" = "Above breakpoint"))+
    theme(plot.title = element_text(size=10, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.title=element_text(size=10),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    ggtitle(s)
  balances2[[s]]<-fig
}

pdf(file="...",width = 7.5, height = 10)
grid.arrange(balances[["PlantHeight"]],
             balances2[["PlantHeight"]],
             balances[["LeafArea"]],
             balances2[["LeafArea"]],
             balances[["LMA"]],
             balances2[["LMA"]],
             balances[["LeafN"]],
             balances2[["LeafN"]],
             balances[["LDMC"]],
             balances2[["LDMC"]],nrow=5)
dev.off()