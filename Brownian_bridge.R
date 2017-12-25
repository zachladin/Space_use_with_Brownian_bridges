#clear environment
rm(list=ls())

#load packages
library(ggmap)
library(cowplot)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(png)
library(reshape)
library(maptools)
library(maps)
library(raster)
library(rgdal)
library(sp)
library(gdalUtils)
library(ggplot2)
library(rasterVis)
library(rgeos)
library(Rmisc)
library(tiff)

#set wd
setwd("/Users/zach/Dropbox (ZachTeam)/Projects/WOTH_GRCA_PostFledging_Movement/")

#load functions
source(paste(getwd(),"TIBB_temporal","R_source","SummaryFunction.R",sep="/"))

###########################################################################
#read in clipped NLCD raster
nlcd.raster<-raster(paste(getwd(),"Results","NLCD_crop.tif",sep="/"))
plot(nlcd.raster)

#get projection
projection(nlcd.raster)
res(nlcd.raster)
nlcd.raster.proj<-projectRaster(nlcd.raster, res=30,crs=CRS("+proj=utm +zone=18 ellps=WGS84"),method="ngb")
############################################################################
#import tiger road shp fileas
roads.de<-readOGR(dsn=paste("/Users/zach/Dropbox (ZachTeam)/Projects/WOTH_GRCA_PostFledging_Movement/TIBB_temporal/FRAME_GIS","Roads","tl_2016_10003_roads",sep="/"), layer="tl_2016_10003_roads")
roads.md<-readOGR(dsn=paste("/Users/zach/Dropbox (ZachTeam)/Projects/WOTH_GRCA_PostFledging_Movement/TIBB_temporal/FRAME_GIS","Roads","tl_2016_24015_roads",sep="/"), layer="tl_2016_24015_roads")
roads.pa<-readOGR(dsn=paste("/Users/zach/Dropbox (ZachTeam)/Projects/WOTH_GRCA_PostFledging_Movement/TIBB_temporal/FRAME_GIS","Roads","tl_2016_42029_roads",sep="/"), layer="tl_2016_42029_roads")

roads.all<-rbind(roads.de, roads.md, roads.pa)

#crop roads to extent of ncld.raster
#roads.proj<-spTransform(roads.all, CRS(projection(nlcd.raster)))
roads.proj<-spTransform(roads.all, CRS("+proj=utm +zone=18 ellps=WGS84"))
roads.crop<-crop(roads.proj, extent(nlcd.raster.proj))

################################################################################
#read in roads raster
roads.raster<-raster(paste(getwd(),"Results","Roads_raster_export_R.tif",sep="/"))
roads.raster.proj<-projectRaster(roads.raster, res=30,crs=CRS("+proj=utm +zone=18 ellps=WGS84"),method="ngb")

###########################################################################
#read in bird location data
data<-read.csv(paste(getwd(),"Data","ALL_MOVEMENT_DATA_MASTER.csv",sep="/"),header=TRUE)

IndList<-sort(unique(data$Individual))

woth.data<-subset(data, Species=="WOTH")
wothIndList<-sort(unique(woth.data$Individual))

###########################################################################
#independent data
ind.data<-subset(data, Stage=="Post-Independence")

woth.ind<-subset(ind.data, Species=="WOTH")
grca.ind<-subset(ind.data, Species=="GRCA")

#remove any locations with NAs
data.1<-ind.data[is.na(ind.data$x)==FALSE,]

#convert UTM to lat/long
#make point data spatialPoints object
coords<-SpatialPoints(data.frame(x=data.1$x, y=data.1$y))

projection(coords)<-"+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000 +y_0=0 +ellps=GRS80 +units=m +no_defs"

#reproject point locations
coords.proj<-spTransform(coords, CRS(frame.projection))

xy<-as.data.frame(coordinates(coords.proj))
head(xy)
par(new=T)
plot(xy$x, xy$y)

#recombine projected points with data
data.1$Longitude<-xy$x
data.1$Latitude<-xy$y

#subset gps data
data.1$bird_num<-as.factor(as.character(data.1$bird_num))

#sort by Individual, then Fix
data.sort<-data.1[order(data.1$bird_num, data.1$day,decreasing = FALSE),]

#plot with ggplot
plot1<-ggplot(data=data.1)+
  aes(x=Longitude, y=Latitude)+
  geom_point(aes(color=as.factor(bird_num)))+
  geom_path(aes(color=as.factor(bird_num)))+
  theme(legend.position="none")
plot1 

#subset plots by species
plot1.facet<-plot1+facet_wrap(Independence~species)
plot1.facet    

#####################################################################################
#get bounding box for all birds

#get extent
minLon=min(data.sort$Longitude,na.rm=TRUE)-0.005
maxLon=max(data.sort$Longitude,na.rm=TRUE)+0.07
minLat=min(data.sort$Latitude,na.rm=TRUE)-0.03
maxLat=max(data.sort$Latitude,na.rm=TRUE)+0.03

bbox.1<-c(left=minLon, bottom=minLat,right=maxLon,top=maxLat)

#get base map
map.base <- get_map(location = bbox.1,maptype="terrain")
ggmap(map.base)
#######################################################
#overlay paths on map

data.sort$Independence<-factor(data.sort$Independence, levels=c("Pre-independence","Post-independence"))
plot2<-ggmap(map.base, darken = c(0.6, "white"))+
  geom_polygon(data=frame.poly.proj.fort,fill=alpha("darkgreen",0.8), aes(x=long, y=lat,group=id))+
  geom_point(data=data.sort, aes(x=Longitude, y=Latitude, color=bird_num),size=0.3)+
  geom_path(data=data.sort, aes(x=Longitude, y=Latitude,color=bird_num))+
  theme(legend.position="none")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  labs(x="Longitude",y="Latitude")

plot2
plot2.facet<-plot2+facet_grid(Independence~species)
plot2.facet

#######################################################################################################################
#load more libraries, but will cause problems with raster pkg, if loaded earlier
library(weights)
library(spatstat)

#####################################################################################################################
#now build funciton to do this for each individual

#define speciesList
speciesList<-c("WOTH","GRCA")

#define typeList (Pre- or Post-Independence)
typeList<-c("independent","dependent")

#test
numIn=15
species="WOTH"
type=typeList[2]

#####################################################################################################################

#extractLandcover function
#build funtion
extractLandcover<-function(numIn, species, type,r ){
  
  proj4string(r) <- CRS("+init=epsg:2776")
  r2 = projectRaster(r, crs=projection(nlcd.raster.proj))
  r3 <- r2 > -Inf
  pp <- rasterToPolygons(r3, dissolve=TRUE)
  
  # Crop data by extent of state subset
  crop.1 <- crop(nlcd.raster.proj, extent(pp))
  crop.2<- raster::mask(crop.1, pp)

  #resample raster r2, to get it in same resolution as crop.2
  s <- raster::resample(r2, crop.2, method='bilinear')
  unique(values(nlcd.raster.proj))
  Landcover <- extract(floor(nlcd.raster.proj), pp)
  landList<-as.numeric(unlist(Landcover))

  pdf.extract<-extract(s,pp)
  weightList<-as.numeric(unlist(pdf.extract))
  length(weightList)
  
  w.prop<-wpct(x=landList, weight=weightList, na.rm=TRUE)
  w.prop.df<-as.data.frame(w.prop)
  w.prop.df$num.codes<-row.names(w.prop.df)
  row.names(w.prop.df)<-NULL
  
  #get area of polygon (sq. km)
  bb.area<-gArea(pp)/1000000
  
  landcoverList<-c("Open Water","Developed Open Space","Developed Low Intensity","Developed Medium Intensity","Developed High Intensity","Barren Land","Deciduous Forest","Evergreen Forest","Mixed Forest", "Scrub/Shrub","Grassland/Herbaceous","Pasture/Hay","Cultivated Crops","Woody Wetlands","Emergent Herbaceous Wetlands")
  
  #convert raster to a data.frame
  raster2.pts <- rasterToPoints(nlcd.raster.proj)
  raster2.df <- data.frame(raster2.pts)
  valTab <- sort(unique(raster2.df[[3]]))
  
  #ratify raster
  rat.nlcd<-ratify(nlcd.raster.proj)
  rat <- levels(rat.nlcd)[[1]]
  rat$code <- c(valTab)
  rat$landcover <- landcoverList
  levels(rat.nlcd) <- rat
  
  # generate land cover number to name conversions
  num.codes <- floor(rat$ID)
  cover.names <- rat$landcover
  conversions <- data.frame(num.codes, cover.names)
  conversions <- na.omit(conversions)
  conversions <- conversions[order(conversions$num.codes),]
  conversions$num.codes<-as.factor(as.character(conversions$num.codes))
  
  #merge w.prop.df with coversions
  conversions.merge<-merge(conversions, w.prop.df, by="num.codes",all.x=TRUE)
  conversions.merge$w.prop[is.na(conversions.merge$w.prop)]<-0
  
  # convert to data frame
  prop.df <- data.frame(Species = species, Individual=paste(species, numIn,sep="_"), Type=type, Area=bb.area, num.codes=num.codes, proportion =  conversions.merge$w.prop)
  colnames(prop.df)<-c("Species","Individual","Type","Total_bridge_area_km2","num.codes","Proportion")
  
  #merge with land cover names
  prop.df.merge<-merge(prop.df,conversions, by=c("num.codes"),all.y=TRUE)
  prop.df.2<-prop.df.merge[,c("Species","Individual","Type","num.codes","cover.names","Total_bridge_area_km2","Proportion")]
  colnames(prop.df.2)<-c("Species","Individual","Type","LandCoverID","Land_Cover_Type","Total_bridge_area_km2","Proportion")
  
  #add Land_Cover_Area_km2
  prop.df.2$Land_Cover_area_km2<-prop.df.2$Total_bridge_area_km2*prop.df.2$Proportion
  
  #sort by land use
  prop.order<-prop.df.2[order(prop.df.2$Proportion, decreasing=FALSE),]
  prop.order$Land_Cover_Type<-factor(prop.order$Land_Cover_Type, levels=c(as.character(prop.order$Land_Cover_Type)))
  
  landcover.out<-prop.order
  return(landcover.out)
}

######################################################################################################################
#extract roads
extractRoads<-function(r){
  
  proj4string(r) <- CRS("+init=epsg:2776")
  r2 = projectRaster(r, crs=projection(nlcd.raster.proj))
  r3 <- r2 > -Inf
  pp <- rasterToPolygons(r3, dissolve=TRUE)
  
  #get roads (shp file)
  #subset roads within bb (pp)
  roads.sub<-roads.crop[pp,]
  rs <- raster(extent(roads.sub), crs=projection(roads.crop),resolution=30)
  rs[] <- 1:ncell(rs)
  # Intersect lines with raster "polygons" and add length to new lines segments
  rsp <- rasterToPolygons(rs)
  rp <- intersect(roads.sub, rsp)
  rp$length <- gLength(rp, byid=TRUE) / 1000
  x <- tapply(rp$length, rp$layer, sum)
  r.road <- raster(rs)
  r.road[as.integer(names(x))] <- x
  r.values<-values(r.road)
  r.crop<-crop(r.road, extent(pp))
  r.mask<-raster::mask(r.crop, pp)
  #resample raster r2, to get it in same resolution as crop.2
  resamp.temp <- raster::resample(r2, r.road, method='bilinear')
  road.den <- extract(r.mask, pp)
  roadList<-as.numeric(unlist(road.den))
  
  resamp.temp.values<-values(resamp.temp)
  
  pdf.extract<-extract(resamp.temp,pp)
  weightList<-as.numeric(unlist(pdf.extract))
  #length(pdf.extract)
  
  road.density<-weighted.mean(x=roadList, weight=weightList, na.rm=TRUE)
  road.density.df<-as.data.frame(road.density)
  colnames(road.density.df)<-"Road_density_km_per_km2"
  row.names(road.density.df)<-NULL
  
  #get roads raster

  # Crop data by extent of state subset
  crop.roads.1 <- crop(roads.raster.proj, extent(pp))
  crop.roads.2<- raster::mask(crop.roads.1, pp)

  #resample raster r2, to get it in same resolution as crop.2
  s <- raster::resample(r2, crop.roads.2, method='bilinear')
  RoadCells <- extract(floor(roads.raster.proj), pp)
  roadList<-as.numeric(unlist(RoadCells))
  
  pdf.roads.extract<-extract(s,pp)
  roads.weightList<-as.numeric(unlist(pdf.roads.extract))
  road.prop<-wpct(x=roadList, weight=roads.weightList, na.rm=TRUE)
  road.prop.df<-as.data.frame(road.prop)
  road.prop.df$num.codes<-row.names(road.prop.df)
  row.names(road.prop.df)<-NULL
  
  road.prop.out<-subset(road.prop.df, num.codes==1)
  
  #get area of polygon (sq. km)
  bb.area<-gArea(pp)/1000000
  
  roads.out.df<-data.frame(Road_proportion=road.prop.out[,1])
  roads.out.df$Road_area_km2<-roads.out.df$Road_proportion*bb.area

  #now combine with Road density
  
  roads.out.df$Road_density_km_per_km2<-road.density.df[,1]

 return(roads.out.df)
}

##############################################################################################################################
#function to get landcover and roads from bridges
extractBridgeInfo<-function(numIn, species, type, file){
  
  #First get bridge for a given individual
  new.file=NULL
  new.file<-file
  new.data<-read.table(new.file)

  out<-new.data[,c("V1","V2","V3")]
  colnames(out) <- c("x","y","average")
  
  #get 95%CI of bridge
  CI.out<-CI(out$average, ci = 0.95)
  lwr<-CI.out[3]
  upr=CI.out[1]
  
  new.data.1<-subset(out, c(average > lwr))
  
  r=rasterFromXYZ(new.data.1)
  
  # proj4string(r) <- CRS("+init=epsg:2776")
  # r2 = projectRaster(r, crs=projection(nlcd.raster))
  # r3 <- r2 > -Inf
  # pp <- rasterToPolygons(r3, dissolve=TRUE)

  #then extract landcover
  landcover=NULL
  landcover<-extractLandcover(numIn=numIn, species=speciesName, type=type, r=r)
  
  #get roads
  roads=NULL
  roads<-tryCatch({extractRoads(r=r)
    },error=function(cond2){
      cond2=data.frame(Road_proportion=0, Road_area_km2=0, Road_density_km_per_km2=0)
      cond2
    })
  
  bridge.out.df<-landcover
  bridge.out.df$Road_proportion<-roads$Road_proportion
  bridge.out.df$Road_area_km2<-roads$Road_area_km2
  bridge.out.df$Road_density_km_per_km2<-roads$Road_density_km_per_km2
  
  return(bridge.out.df)
}

###############################################################################################
#test missing road.density birds



###############################################################################################
#typeList
typeList<-c("dependent","independent")

#speciesList
speciesList<-c("GRCA","WOTH")

###############################################################################################

#for loop with function
all.results.out<-list()
for(k in 1:length(typeList)){
  type=typeList[k]
  
  dir.use<-ifelse(type=="dependent",
                  paste("/Users/zach/Dropbox (ZachTeam)/Projects/WOTH_GRCA_PostFledging_Movement/TIBB_temporal/TIBB_output","dependent",sep="/"),
                  paste("/Users/zach/Dropbox (ZachTeam)/Projects/WOTH_GRCA_PostFledging_Movement/TIBB_temporal/TIBB_output","independent",sep="/"))
  setwd(dir.use)
  
  prop.all.out<-list()
  for(j in 1:length(speciesList)){
  
    speciesName<-unique(as.character(speciesList[j]))

    fileList<-list.files(pattern="\\TIBB.txt$")
    
    ifelse(speciesName=="GRCA",
           fileList<-fileList[substr(fileList,1,1)=="G"],
           fileList<-fileList[substr(fileList,1,1)=="W"])
    
    species.out<-list()
    for(i in 1:length(fileList)){
    print(i)
      
    speciesName2<-substr(fileList[i],1,4)
    num.ind<-as.numeric(gsub("\\D", "", fileList[i])) 
    
    fileName<-paste("/Users/zach/Dropbox (ZachTeam)/Projects/WOTH_GRCA_PostFledging_Movement/TIBB_temporal/TIBB_output",type,fileList[i],sep="/")
    
    ind.results<-tryCatch({extractBridgeInfo(numIn=num.ind,species=speciesName2,type=type,file=fileName)
      
     },error=function(cond2){
       cond2=data.frame(
         Species=speciesName2,Individual=paste(speciesName2,num.ind,sep="_"),Type=type, LandCoverID=NA,Land_Cover_Type=NA,Total_bridge_area_km2=NA,Proportion=NA, Land_Cover_area_km2=NA, Road_proportion=NA, Road_area_km2=NA, Road_density_km_per_km2=NA)
       cond2
     })
   
      species.out<-rbind(species.out, ind.results)
    }
    prop.all.out<-rbind(prop.all.out, species.out)
  }
all.results.out<-rbind(all.results.out, prop.all.out)
}

write.csv(all.results.out, file=paste("/Users/zach/Dropbox (ZachTeam)/Projects/WOTH_GRCA_PostFledging_Movement/TIBB_temporal","Results","all_landcover_estimates_new_use_fileList_UTM.csv",sep="/"),row.names=FALSE)  

#########################################################################################
#########################################################################################
#########################################################################################
#summary function creates table of n, mean, var, SD, and SE
summaryFunction <- function(DataIn, factor, response){
  require(plyr)
  summaryOut <- ddply(DataIn, factor, .fun = function(xx){
    c(n = length(xx[,response]),
      mean = mean(xx[,response],na.rm=TRUE),
      var = var(xx[,response],na.rm=TRUE),
      SD = sd(xx[,response],na.rm=TRUE),
      SE = sqrt(var(xx[,response])/length(xx[,response])))
  })
  return(summaryOut)
  dev.off()
}
#########################################################################################
#read in results
all.results.in<-read.csv(paste("/Users/zach/Dropbox (ZachTeam)/Projects/WOTH_GRCA_PostFledging_Movement/TIBB_temporal","Results","all_landcover_estimates_new_use_fileList_UTM.csv",sep="/"),header=TRUE)

#remove any NAs (individuals w/o enough data)
all.results<-all.results.in[!is.na(all.results.in$LandCoverID),]

#create Landcover.Species.Individual column
all.results$Landcover.Species.Individual.Type<-paste(all.results$Land_Cover_Type, all.results$Individual,all.results$Type,sep="_")

#create species.individual columnn
all.results$species.individual.type<-paste( all.results$Individual,all.results$Type, sep="_")

#get unique Landcover.Species.Individual
individuals<-data.frame(individual=unique(all.results$species.individual.type))

#make dataframe with all posible landcover types for both species
land.cover.df<-data.frame(Land_Cover_Type=c("Open Water","Developed Open Space","Developed Low Intensity","Developed Medium Intensity","Developed High Intensity","Barren Land","Deciduous Forest","Evergreen Forest","Mixed Forest", "Scrub/Shrub","Grassland/Herbaceous","Pasture/Hay","Cultivated Crops","Woody Wetlands","Emergent Herbaceous Wetlands"))

#merge two dataframes
merge1<-merge(individuals, land.cover.df,all=TRUE)

#create Landcover.Species.Individual column
merge1$Landcover.Species.Individual.Type<-paste(merge1$Land_Cover_Type, merge1$individual,sep="_")


#now merge with all.results
merge2<-merge(merge1, all.results, by=c("Landcover.Species.Individual.Type"),all.x=TRUE)

drops <- c("Land_Cover_Type.y","Individual","LandCoverID","Species","species.individual.type")
merge.out<-merge2[ , !(names(merge2) %in% drops)]
merge.out$Landcover.Species.Individual.Type<-as.factor(as.character(merge.out$Landcover.Species.Individual.Type))

#change NAs to zeros
merge.out[is.na(merge.out)]<-0

#get separate Landcover.Species.Individual.Type
separate<-read.table(text=as.character(merge.out$Landcover.Species.Individual.Type), colClasses="character",sep="_")
colnames(separate)<-c("Landcover","Species","Indivdual","Type")

#add species column 
merge.out$Species<-separate$Species
#overwrite Type column
merge.out$Type<-separate$Type

#renames landcover column
colnames(merge.out)[3]<-"Land_Cover_Type"

merge.out$Species.Landcover.Type<-paste(merge.out$Species, merge.out$Land_Cover_Type,merge.out$Type,sep="_")
merge.out$Species.Landcover.Type<-as.factor(merge.out$Species.Landcover.Type)

#########################################################################
#use weighted.mean
#weighted.mean(x,y, pdf)

#use aggregate to get mean proportion
woth.test<-subset(merge.out, c(Species=="WOTH" & Type=="independent"))
mean.1<-aggregate(Proportion~Species.Landcover.Type, FUN="mean",data=woth.test)
sum(mean.1$Proportion)


#remove GRCA 50 and 16 and WOTH 35
merge.out.sub<-merge.out[c(merge.out$individual!="GRCA_50_dependent" & merge.out$individual !="GRCA_16_dependent" &
                             merge.out$individual !="WOTH_35_independent"),]


#summaryFunction
require(plyr)
summary.1<-summaryFunction(DataIn=merge.out, factor="Species.Landcover.Type",response="Land_Cover_area_km2")

#break apart Species.Landcover
factors<-read.table(text=as.character(summary.1$Species.Landcover), sep="_",colClasses="character")
summary.2<-cbind(factors, summary.1)
colnames(summary.2)<-c("Species","Land_Cover_Type","Type","Species.Landcover.Type","n","Mean","var","SD","SE")
head(summary.2)

#sort by GRCA highest to lowest Land cover area
summary.2.order<-summary.2[order(summary.2$Mean, summary.2$Species),]
summary.2.order$Species<-as.factor(summary.2.order$Species)

#remove NAs
summary.3<-summary.2.order

#set factor levels for habitat types
lc.list<-data.frame(Land_Cover_Type=c("Developed High Intensity","Developed Medium Intensity","Developed Low Intensity","Developed Open Space","Pasture/Hay","Cultivated Crops","Barren Land","Grassland/Herbaceous","Scrub/Shrub","Deciduous Forest","Evergreen Forest","Mixed Forest", "Woody Wetlands", "Emergent Herbaceous Wetlands", "Open Water"))


summary.3$Land_Cover_Type<-factor(summary.3$Land_Cover_Type, levels=rev(as.character(lc.list$Land_Cover_Type)))

summary.3$Type<-as.factor(as.character(summary.3$Type))

levels(summary.3$Land_Cover_Type)

#change factor levels
levels(summary.3$Type)
levels(summary.3$Type) <- c("Dependent", "Independent")


#change Species factor order
summary.3$Species<-as.factor(as.character(summary.3$Species))
levels(summary.3$Species)
summary.3$Species<-factor(summary.3$Species, levels=c("WOTH","GRCA"))
levels(summary.3$Species)



#mycolors
mycolors<-c("sienna2","royalblue2")

#now plot
#look at proportions plotted
area.plot2<-ggplot(data=summary.3, aes(x=Land_Cover_Type, y=Mean, ymin=Mean-SE, ymax=Mean+SE))+
  coord_flip()+
  geom_errorbar(width=0.1, size=0.5, position=position_dodge(width=0.95),aes(color=Species))+
  geom_bar(stat="identity",position="dodge",aes(fill=Species), alpha=0.7)+
  scale_color_manual(values=mycolors)+
  scale_fill_manual(values=mycolors,breaks=c("GRCA", "WOTH"), labels=c("GRCA", "WOTH"))+
  theme(panel.border = element_rect(color="black", size=1, linetype="solid"))+
  theme(legend.position = c(0.89, 0.99), 
        legend.justification = c(0, 1.2), 
        legend.background = element_rect(fill="transparent",color="black",size=0.3, linetype="solid"))+
  ylab(bquote('Area ('*km^2*')'))+
labs(x="Land cover type")+
  guides(color=FALSE)
  #geom_text(aes(label=Land_Cover_Type))
  #ggtitle("2011 NLCD underlying mean Wood Thrush\nBrownian bridge extent")
area.plot2

#now facet by pre- and post-independence
area.plot.facet<-area.plot2+facet_grid(.~Type)
area.plot.facet

#area.plot2.facet<-area.plot2+facet_grid(.~Species)

ggsave(area.plot.facet, file = paste(getwd(),"Figures","All_LandCover_areas.png",sep="/"), dpi=300, width=7,height=7, limitsize=FALSE)

#############################################################################################################################
#now look at relative proportions

#summaryFunction
summary.prop.1<-summaryFunction(DataIn=merge.out, factor="Species.Landcover.Type",response="Proportion")

#break apart Species.Landcover
factors<-read.table(text=as.character(summary.prop.1$Species.Landcover), sep="_",colClasses="character")
summary.prop.2<-cbind(factors, summary.prop.1)
colnames(summary.prop.2)<-c("Species","Land_Cover_Type","Type","Species.Landcover.Type","n","Mean","var","SD","SE")

#sort by GRCA highest to lowest Land cover area
summary.prop.2.order$Species<-as.factor(summary.prop.2.order$Species)
summary.prop.2.order<-summary.prop.2[order(summary.prop.2$Mean, summary.prop.2$Species),]
#summary.prop.2[names(summary.prop.2)=="Land_Cover_Type.x"]<-"Land_Cover_Type"
#remove NAs
summary.prop.3<-summary.prop.2.order

#set factor levels for habitat types
lc.list<-data.frame(Land_Cover_Type=c("Developed High Intensity","Developed Medium Intensity","Developed Low Intensity","Developed Open Space","Pasture/Hay","Cultivated Crops","Barren Land","Grassland/Herbaceous","Scrub/Shrub","Deciduous Forest","Evergreen Forest","Mixed Forest", "Woody Wetlands", "Emergent Herbaceous Wetlands", "Open Water"))

summary.prop.3$Land_Cover_Type<-factor(summary.prop.3$Land_Cover_Type, levels=rev(as.character(lc.list$Land_Cover_Type)))

summary.prop.3$Species<-as.factor(as.character(summary.prop.3$Species))
summary.prop.3$Type<-as.factor(as.character(summary.prop.3$Type))


#change factor levels
levels(summary.prop.3$Type)
levels(summary.prop.3$Type) <- c("Dependent", "Independent")

#change Species factor order
summary.prop.3$Species<-as.factor(as.character(summary.prop.3$Species))
levels(summary.prop.3$Species)
summary.prop.3$Species<-factor(summary.prop.3$Species, levels=c("WOTH","GRCA"))

BroadLand<-data.frame(Land_Cover_Type=lc.list,BroadLand=c("Developed","Developed","Developed","Developed","Ag & grassland","Ag & grassland","Ag & grassland","Ag & grassland","Forest","Forest","Forest","Forest","Forest","Wetland", "Wetland"))


#create some broad categories
summary.prop.4<-merge(summary.prop.3, BroadLand, by="Land_Cover_Type",all.x=TRUE)


#mycolors
mycolors<-c("sienna2", "royalblue2")

#now plot
#look at proportions plotted
prop.plot3<-ggplot(data=summary.prop.4, aes(x=Land_Cover_Type, y=Mean, ymin=Mean-SE, ymax=Mean+SE))+
  coord_flip()+
  geom_errorbar(width=0.1, size=0.5, position=position_dodge(width=0.95),aes(color=Species))+
  geom_bar(stat="identity",position="dodge",aes(fill=Species), alpha=0.7)+
  scale_color_manual(values=mycolors)+
  scale_fill_manual(values=mycolors,breaks=c("GRCA", "WOTH"), labels=c("GRCA", "WOTH"))+
  theme(panel.border = element_rect(color="black", size=1, linetype="solid"))+
  theme(legend.position = c(0.89, 0.99), 
        legend.justification = c(0, 1.2), 
        legend.background = element_rect(fill="transparent",color="black",size=0.3, linetype="solid"))+
  guides(color=FALSE)+
  ylab(bquote('Area ('*km^2*')'))+
  labs(x="Land cover type", y="Weighted relative proportion")+
  ylim(0,1)

#ggtitle("2011 NLCD underlying mean Wood Thrush\nBrownian bridge extent")
prop.plot3


prop.plot2.facet<-prop.plot3+facet_grid(.~Type)
prop.plot2.facet


ggsave(prop.plot3, file = paste(getwd(),"Figures","All_LandCover_Proportions.png",sep="/"), dpi=300, width=7,height=7, limitsize=FALSE)

###############################################################################################
#now look at road area

#make a Species.Type column
merge.out$Species.Type<-paste(merge.out$Species, merge.out$Type, sep=".")

road.df<-unique(merge.out[,c("Species.Type","Road_area_km2")])

#summaryFunction
summary.road.1<-summaryFunction(DataIn=road.df, factor="Species.Type",response="Road_area_km2")

#break apart Species.Landcover
factors<-read.table(text=as.character(summary.road.1$Species.Type), sep=".",colClasses="character")
summary.road.2<-cbind(factors, summary.road.1)
head(summary.road.2)
colnames(summary.road.2)<-c("Species","Type","Species.Type","n","Mean","var","SD","SE")

#sort by GRCA highest to lowest Land cover area
summary.road.2.order<-summary.road.2[order(summary.road.2$Mean, summary.road.2$Species),]
summary.road.2.order$Species<-as.factor(summary.road.2.order$Species)
summary.road.2[names(summary.road.2)=="Land_Cover_Type.x"]<-"Land_Cover_Type"
#remove NAs
summary.road.3<-summary.road.2.order


summary.road.3$Species<-as.factor(as.character(summary.road.3$Species))
summary.road.3$Type<-as.factor(as.character(summary.road.3$Type))


#change factor levels
levels(summary.road.3$Type)
levels(summary.road.3$Type) <- c("Dependent", "Independent")

#change Species factor order
levels(summary.road.3$Species)
summary.road.3$Species<-factor(summary.road.3$Species, levels=c("GRCA","WOTH"))
levels(summary.road.3$Species)


#mycolors
mycolors<-c("royalblue2","sienna2")

#now plot
#look at Road density plotted
road.plot3<-ggplot(data=summary.road.3, aes(x=Species, y=Mean, ymin=Mean-SE, ymax=Mean+SE))+
  #coord_flip()+
  geom_errorbar(width=0.1, size=0.5, position=position_dodge(width=0.95),aes(color=Species))+
  geom_bar(stat="identity",position="dodge",aes(fill=Species), alpha=0.8)+
  scale_color_manual(values=mycolors)+
  scale_fill_manual(values=mycolors,breaks=c("GRCA", "WOTH"), labels=c("GRCA", "WOTH"))+
  theme(panel.background = element_rect(fill="white"),panel.border = element_rect(color="black", size=1, linetype="solid"))+
  theme(legend.position = c(0.89, 0.99), 
        legend.justification = c(0, 1.2), 
        legend.background = element_rect(fill="transparent",color="black",size=0.3, linetype="solid"))+
  #guides(color=FALSE, alpha=FALSE)+
  ylab(bquote('Road area km2'))+
  labs(x="Species")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
  #ylim(0,1)

#ggtitle("2011 NLCD underlying mean Wood Thrush\nBrownian bridge extent")
road.plot3


road.plot2.facet<-road.plot3+facet_grid(.~Type)
road.plot2.facet



ggsave(road.plot2.facet, file = paste(getwd(),"Figures","All_Road_Area.png",sep="/"), dpi=300, width=7,height=7, limitsize=FALSE)

###############################################################################################
#now look at road proportions

#make a Species.Type column
merge.out$Species.Type<-paste(merge.out$Species, merge.out$Type, sep=".")

road.df<-unique(merge.out[,c("Species.Type","Road_proportion")])

#summaryFunction
summary.road.1<-summaryFunction(DataIn=road.df, factor="Species.Type",response="Road_proportion")

#break apart Species.Landcover
factors<-read.table(text=as.character(summary.road.1$Species.Type), sep=".",colClasses="character")
summary.road.2<-cbind(factors, summary.road.1)
head(summary.road.2)
colnames(summary.road.2)<-c("Species","Type","Species.Type","n","Mean","var","SD","SE")

#sort by GRCA highest to lowest Land cover area
summary.road.2.order<-summary.road.2[order(summary.road.2$Mean, summary.road.2$Species),]
summary.road.2.order$Species<-as.factor(summary.road.2.order$Species)
summary.road.2[names(summary.road.2)=="Land_Cover_Type.x"]<-"Land_Cover_Type"
#remove NAs
summary.road.3<-summary.road.2.order


summary.road.3$Species<-as.factor(as.character(summary.road.3$Species))
summary.road.3$Type<-as.factor(as.character(summary.road.3$Type))


#change factor levels
levels(summary.road.3$Type)
levels(summary.road.3$Type) <- c("Dependent", "Independent")

#change Species factor order
levels(summary.road.3$Species)
summary.road.3$Species<-factor(summary.road.3$Species, levels=c("GRCA","WOTH"))
levels(summary.road.3$Species)


#mycolors
mycolors<-c("royalblue2","sienna2")

#now plot
#look at Road density plotted
road.plot3<-ggplot(data=summary.road.3, aes(x=Species, y=Mean, ymin=Mean-SE, ymax=Mean+SE))+
  #coord_flip()+
  geom_errorbar(width=0.1, size=0.5, position=position_dodge(width=0.95),aes(color=Species))+
  geom_bar(stat="identity",position="dodge",aes(fill=Species), alpha=0.7)+
  scale_color_manual(values=mycolors)+
  scale_fill_manual(values=mycolors,breaks=c("GRCA", "WOTH"), labels=c("GRCA", "WOTH"))+
  theme(panel.background = element_rect(fill="white"),panel.border = element_rect(color="black", size=1, linetype="solid"))+
  theme(legend.position = c(0.89, 0.99), 
        legend.justification = c(0, 1.2), 
        legend.background = element_rect(fill="transparent",color="black",size=0.3, linetype="solid"))+
  #guides(color=FALSE, alpha=FALSE)+
  ylab(bquote('Weighted relative proportion of roads'))+
  labs(x="Species")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
#ylim(0,1)

#ggtitle("2011 NLCD underlying mean Wood Thrush\nBrownian bridge extent")
road.plot3


road.plot2.facet<-road.plot3+facet_grid(.~Type)
road.plot2.facet

ggsave(road.plot2.facet, file = paste(getwd(),"Figures","All_Road_Proportions.png",sep="/"), dpi=300, width=7,height=7, limitsize=FALSE)


###############################################################################################
#now look at road density (km per km2)

#make a Species.Type column
merge.out$Species.Type<-paste(merge.out$Species, merge.out$Type, sep=".")

road.df<-unique(merge.out[,c("Species.Type","Road_density_km_per_km2")])

#summaryFunction
summary.road.1<-summaryFunction(DataIn=road.df, factor="Species.Type",response="Road_density_km_per_km2")

#break apart Species.Landcover
factors<-read.table(text=as.character(summary.road.1$Species.Type), sep=".",colClasses="character")
summary.road.2<-cbind(factors, summary.road.1)
head(summary.road.2)
colnames(summary.road.2)<-c("Species","Type","Species.Type","n","Mean","var","SD","SE")

#sort by GRCA highest to lowest Land cover area
summary.road.2.order<-summary.road.2[order(summary.road.2$Mean, summary.road.2$Species),]
summary.road.2.order$Species<-as.factor(summary.road.2.order$Species)
summary.road.2[names(summary.road.2)=="Land_Cover_Type.x"]<-"Land_Cover_Type"
#remove NAs
summary.road.3<-summary.road.2.order


summary.road.3$Species<-as.factor(as.character(summary.road.3$Species))
summary.road.3$Type<-as.factor(as.character(summary.road.3$Type))


#change factor levels
levels(summary.road.3$Type)
levels(summary.road.3$Type) <- c("Dependent", "Independent")

#change Species factor order
levels(summary.road.3$Species)
summary.road.3$Species<-factor(summary.road.3$Species, levels=c("GRCA","WOTH"))
levels(summary.road.3$Species)


#mycolors
mycolors<-c("royalblue2","sienna2")

#now plot
#look at Road density plotted
road.plot3<-ggplot(data=summary.road.3, aes(x=Species, y=Mean, ymin=Mean-SE, ymax=Mean+SE))+
  #coord_flip()+
  geom_errorbar(width=0.1, size=0.5, position=position_dodge(width=0.95),aes(color=Species))+
  geom_bar(stat="identity",position="dodge",aes(fill=Species), alpha=0.8)+
  scale_color_manual(values=mycolors)+
  scale_fill_manual(values=mycolors,breaks=c("GRCA", "WOTH"), labels=c("GRCA", "WOTH"))+
  theme(panel.background = element_rect(fill="white"),panel.border = element_rect(color="black", size=1, linetype="solid"))+
  theme(legend.position = c(0.89, 0.99), 
        legend.justification = c(0, 1.2), 
        legend.background = element_rect(fill="transparent",color="black",size=0.3, linetype="solid"))+
  #guides(color=FALSE, alpha=FALSE)+
  ylab(bquote('Road density (km per km2)'))+
  labs(x="Species")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
#ylim(0,1)

#ggtitle("2011 NLCD underlying mean Wood Thrush\nBrownian bridge extent")
road.plot3


road.plot2.facet<-road.plot3+facet_grid(.~Type)
road.plot2.facet

ggsave(road.plot2.facet, file = paste(getwd(),"Figures","All_LandCover_Road_Density.png",sep="/"), dpi=300, width=7,height=7, limitsize=FALSE)


#################################################################################################
