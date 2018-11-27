#Forest RedList species
require(rgbif)
require(raster)
require(data.table)
require(rasterVis)
require(dismo)

#Norway polygon
norway<-raster::getData('GADM', country='NOR', level=0)
norwayP<-spTransform(norway,"+proj=utm +zone=32")
plot(norwayP)

# #Redlist -----------------------------------------------------------------
# #Make species lists with those linked to grazing changes
# #Data from Artsdatabanken Rødlist 2015. All redlisted spp, påvirkning factor: changed habitat use
# #Only species (not subspecies)
# #Filter by 'beite' 
vasc<-read.table('Rodlista2015_vascplant.csv'
                  ,header=T,sep=";",fileEncoding="UTF-16LE")
 
herbvasc<-droplevels(vasc[grep('*beit*',vasc$Påvirkningsfaktorer,ignore.case=T),])
dim(herbvasc)
herbvasc$Vitenskapelig.navn

moss<-read.table('Rodlista2015_mosses.csv'
                 ,header=T,sep=";",fileEncoding="UTF-16LE")

herbmoss<-droplevels(moss[grep('*beit*',moss$Påvirkningsfaktorer,ignore.case=T),])
dim(herbmoss)
herbmoss$Vitenskapelig.navn

lichen<-read.table('Rodlista2015_lichens.csv'
                 ,header=T,sep=";",fileEncoding="UTF-16LE")

herblichen<-droplevels(lichen[grep('*beit*',lichen$Påvirkningsfaktorer,ignore.case=T),])
dim(herblichen)
herblichen$Vitenskapelig.navn


redlistsp_all<-rbind(vasc,moss,lichen)
redlistsp_all$PlantGroup<-c(rep('Vascular',times=nrow(vasc)),rep('Bryphyte',times=nrow(moss)),rep('Lichen',times=nrow(lichen)))

View(redlistsp_all)


# GBIF import -------------------------------------------------------------

#Keyed download with doi linkk
keysL<-sapply(as.character(herblichen$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)
paste(keysL,collapse=',')#Copy and paste to below
odLichen<-occ_download('taxonKey = 2609133,2607387,3408876,2603661,2602561,2609349,3422592,3433876,2601243,2608140,5260765,7083298,2605872,7425899,5260770,9198278,5258394,8704544,3397573,3429282,3390552,2609409,2609180,5260761,2607725,8335421,5260565,2599762,3429909,3389602'
                       ,'country = NO','hasCoordinate = TRUE',
                       user='jamesspeed',pwd='*****',email='*****')
occ_download_meta(odLichen)
gbif_citation(occ_download_meta(odLichen))# GBIF Occurrence Download https://doi.org/10.15468/dl.pl144i Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"

keysM<-sapply(as.character(herbmoss$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)
paste(keysM,collapse=',')
odMoss<-occ_download('taxonKey = 5280585,4276928,5283214,7800507,2689413'
                       ,'country = NO','hasCoordinate = TRUE',
                       user='jamesspeed',pwd='*****',email='*****')
occ_download_meta(odMoss)
gbif_citation(occ_download_meta(odMoss))# "GBIF Occurrence Download https://doi.org/10.15468/dl.ydfbtt Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"

keysV<-sapply(as.character(herbvasc$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)
paste(keysV,collapse=',')
odVasc<-occ_download('taxonKey = 9139754,8558004,2704845,3025813,7931051,2792588,9020552,5405976,2849252,8231647,3001509,8917443,2820517,3012376,2975152,8277403,3033129,5361866,5284517,8869754,7270427,5410857,8915737,2926086,5409958,2787993,3111049,9177060,2927078,3012509,2926055,7660935,2914396,2975380,2925944'
                     ,'country = NO','hasCoordinate = TRUE',
                     user='jamesspeed',pwd='*****',email='*****')
occ_download_meta(odVasc)
gbif_citation(occ_download_meta(odVasc))# "GBIF Occurrence Download https://doi.org/10.15468/dl.7eqijw Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"

odLichdat<-occ_download_get(odLichen,overwrite=T)
odMossdat<-occ_download_get(odMoss,overwrite=T)
odVascdat<-occ_download_get(odVasc,overwrite=T)

odLichendata<-occ_download_import(odLichdat)
odMossdata<-occ_download_import(odMossdat)
odVascdata<-occ_download_import(odVascdat)

summary(as.factor(odVascdata$species))
summary(as.factor(odMossdata$species))
summary(as.factor(odLichendata$species))

#Convert to SPDF
mossspdf<-SpatialPointsDataFrame(cbind(odMossdata$decimalLongitude,odMossdata$decimalLatitude),as.data.frame(odMossdata),proj4string = crs(norway))
lichenspdf<-SpatialPointsDataFrame(cbind(odLichendata$decimalLongitude,odLichendata$decimalLatitude),as.data.frame(odLichendata),proj4string = crs(norway))
vascspdf<-SpatialPointsDataFrame(cbind(odVascdata$decimalLongitude,odVascdata$decimalLatitude),as.data.frame(odVascdata),proj4string = crs(norway))

moss_utm<-spTransform(mossspdf,crs(norwayP))
lichen_utm<-spTransform(lichenspdf,crs(norwayP))
vasc_utm<-spTransform(vascspdf,crs(norwayP))

plot(norwayP)
points(vasc_utm,pch=16,cex=0.1)


# Clipping to Norway outline bufferd --------------------------------------

#Clip to remove points outside of Norway
#Polygon of 1km buffer around Norway
library(rgeos)
norway1km<-gBuffer(norwayP, width=1000)
#plot(norway1km)
#plot(norwayP,add=T)

#Clip
lichen_clip<-lichen_utm[which(!is.na(over(lichen_utm,norway1km))),]
moss_clip<-moss_utm[which(!is.na(over(moss_utm,norway1km))),]
vasc_clip<-vasc_utm[which(!is.na(over(vasc_utm,norway1km))),]
plot(norwayP)
points(lichen_utm,pch=16,col=2)
points(lichen_clip,pch=16,col=3)
plot(norwayP)
points(moss_utm,pch=16,col=2)
points(moss_clip,pch=16,col=3)
plot(norwayP)
points(vasc_utm,pch=16,col=2)
points(vasc_clip,pch=16,col=3)

AllForestRedList<-rbind(lichen_clip,moss_clip,vasc_clip)
AllForestRedList$PlantGroup<-c(rep('Lichen',times=nrow(lichen_clip)),rep('Bryophyte',times=nrow(moss_clip)),rep('Vascular',times=nrow(vasc_clip)))
plot(norwayP)
points(AllForestRedList[AllForestRedList$PlantGroup=='Lichen',])

#Merge with red list details
match(AllForestRedList$species,redlistsp_all$Vitenskapelig.navn)
ForestRedList_adb<-merge(AllForestRedList,redlistsp_all,by.x='species',by.y='Vitenskapelig.navn',all.x=T)


# Final spp data ----------------------------------------------------------

write.table(lichen_clip,'Lichens_forest_redlisted_herbivory.csv')
write.table(moss_clip,'Moss_forest_redlisted_herbivory.csv')
write.table(vasc_clip,'Vascular_forest_redlisted_herbivory.csv')

fwrite(ForestRedList_adb@data,file='RedListedForestSpeciesNorwayBeite.csv')



# Environmental data ------------------------------------------------------


#Elevation
Norelev<-getData('alt',country='NOR')
plot(Norelev)

# #Climate
# Norbioclim<-getData('worldclim',var='bio',res=0.5,lon=5,lat=60)
# Norbioclim1<-getData('worldclim',var='bio',res=0.5,lon=5,lat=70)
# Norbioclim2<-getData('worldclim',var='bio',res=0.5,lon=40,lat=70)
# plot(Norbioclim[[1]])
# plot(Norbioclim1[[1]])
# mergclim<-merge(Norbioclim,Norbioclim1)
# mergclim1<-merge(mergclim,Norbioclim2)
# cropclim<-crop(mergclim1,Norelev)
# Norclimdat<-mask(cropclim,Norelev)
# plot(Norclimdat[[1]])
# 
# NorClimElev<-stack(Norclimdat,Norelev)
# names(NorClimElev)<-c(names(Norbioclim),names(Norelev))
# writeRaster(NorClimElev,'NorClimElev')
# 
# 
# NorClimElev<-stack('NorClimElev')
# NorClimElev
#  
# 
# #Land cover data
# landcover<-stack('landcoverstack')
# 
# 
# #Climate and elevation in same stack
# NorClimElev_utm<-projectRaster(NorClimElev,landcover[[1]],method='bilinear')
# 
# #Soil data
# soilph<-raster('geonode_phihox_m_sl2_250m.tif')
# soilph_utm<-projectRaster(soilph,landcover[[1]],method='bilinear')
# Norsoilph<-mask(crop(soilph_utm,landcover[[1]]),landcover[[1]])
# plot(Norsoilph)
# 
# #Herbivore data
# herbdat<-readOGR('KommuneMetabolicBiomass')
# herbdat1<-spTransform(herbdat,crs(landcover))
# moose1949<-rasterize(herbdat1,field='moose.1949',landcover[[1]])
# moose1959<-rasterize(herbdat1,field='moose.1959',landcover[[1]])
# moose1969<-rasterize(herbdat1,field='moose.1969',landcover[[1]])
# moose1979<-rasterize(herbdat1,field='moose.1979',landcover[[1]])
# moose1989<-rasterize(herbdat1,field='moose.1989',landcover[[1]])
# moose1999<-rasterize(herbdat1,field='moose.1999',landcover[[1]])
# moose2009<-rasterize(herbdat1,field='moose.2009',landcover[[1]])
# moose2015<-rasterize(herbdat1,field='moose.2015',landcover[[1]])
# moosestack<-stack(moose1949,moose1959,moose1969,moose1979,moose1989,moose1999,moose2009,moose2015)
# names(moosestack)<-c('moose1949','moose1959','moose1969','moose1979','moose1989','moose1999','moose2009','moose2015')
# writeRaster(moosestack,'Moose_metbio',overwrite=T)
# 
# red_deer1949<-rasterize(herbdat1,field='red_deer.1949',landcover[[1]])
# red_deer1959<-rasterize(herbdat1,field='red_deer.1959',landcover[[1]])
# red_deer1969<-rasterize(herbdat1,field='red_deer.1969',landcover[[1]])
# red_deer1979<-rasterize(herbdat1,field='red_deer.1979',landcover[[1]])
# red_deer1989<-rasterize(herbdat1,field='red_deer.1989',landcover[[1]])
# red_deer1999<-rasterize(herbdat1,field='red_deer.1999',landcover[[1]])
# red_deer2009<-rasterize(herbdat1,field='red_deer.2009',landcover[[1]])
# red_deer2015<-rasterize(herbdat1,field='red_deer.2015',landcover[[1]])
# red_deerstack<-stack(red_deer1949,red_deer1959,red_deer1969,red_deer1979,red_deer1989,red_deer1999,red_deer2009,red_deer2015)
# names(red_deerstack)<-c('red_deer1949','red_deer1959','red_deer1969','red_deer1979','red_deer1989','red_deer1999','red_deer2009','red_deer2015')
# writeRaster(red_deerstack,'red_deer_metbio',overwrite=T)
# 
# roe_deer1949<-rasterize(herbdat1,field='roe_deer.1949',landcover[[1]])
# roe_deer1959<-rasterize(herbdat1,field='roe_deer.1959',landcover[[1]])
# roe_deer1969<-rasterize(herbdat1,field='roe_deer.1969',landcover[[1]])
# roe_deer1979<-rasterize(herbdat1,field='roe_deer.1979',landcover[[1]])
# roe_deer1989<-rasterize(herbdat1,field='roe_deer.1989',landcover[[1]])
# roe_deer1999<-rasterize(herbdat1,field='roe_deer.1999',landcover[[1]])
# roe_deer2009<-rasterize(herbdat1,field='roe_deer.2009',landcover[[1]])
# roe_deer2015<-rasterize(herbdat1,field='roe_deer.2015',landcover[[1]])
# roe_deerstack<-stack(roe_deer1949,roe_deer1959,roe_deer1969,roe_deer1979,roe_deer1989,roe_deer1999,roe_deer2009,roe_deer2015)
# names(roe_deerstack)<-c('roe_deer1949','roe_deer1959','roe_deer1969','roe_deer1979','roe_deer1989','roe_deer1999','roe_deer2009','roe_deer2015')
# writeRaster(roe_deerstack,'roe_deer_metbio',overwrite=T)
# 
# #Stack together
# PredVars<-stack(NorClimElev_utm,landcover,Norsoilph,moosestack,red_deerstack,roe_deerstack)
# writeRaster(PredVars,'PredictorVariables',overwrite=T)


# Setup -------------------------------------------------------------------


PredVars<-stack('PredictorVariables')
names(PredVars)[21:26]<-c('Land_Cover',"Land_Use",'Forest_Productivity','Forest_Type','Vegetation_Type','SoilpH')


#Ratify land cover
PredVars$Land_Cover<-ratify(PredVars$Land_Cover)
ratlc<- levels(PredVars$Land_Cover)[[1]]
ratlc[["Land_Cover"]] <- c("Built up","Agricultural","Forest","Natural vegetation","Mires","Glaciers/Ice/Snow","Freshwater")
levels(PredVars$Land_Cover)<-ratlc
levelplot(PredVars$Land_Cover)

#Ratify forest productivty
PredVars$Forest_Productivity[PredVars$Forest_Productivity>18]<-NA
PredVars$Forest_Productivity<-ratify(PredVars$Forest_Productivity)
ratlcp<-levels(PredVars$Forest_Productivity)[[1]]
ratlcp[['Forest_Productivity']]<-c('Unproductive','Low','Medium','High')
levels(PredVars$Forest_Productivity)<-ratlcp
levelplot(PredVars$Forest_Productivity)
pairs(PredVars)

#Ratify forest type
PredVars$ForestType[PredVars$ForestType>33]<-NA
PredVars$ForestType<-ratify(PredVars$ForestType)
ratlct<-levels(PredVars$ForestType)[[1]]
ratlct
ratlct[['ForestType']]<-c('Coniferous','Deciduous','Mixed')
levels(PredVars$ForestType)<-ratlct
levelplot(PredVars$ForestType)


