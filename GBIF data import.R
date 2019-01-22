#Forest RedList species
require(rgbif)#Spatial data
require(raster)#Spatial data
require(rgdal)#Spatial data
require(sp)#Spatial data
require(data.table)#Reading in big tables
require(rasterVis)#Visualising spatial data
require(dismo)#Distribution modelling
require(ENMeval)#Tuning maxent

#Norway polygon
norway<-raster::getData('GADM', country='NOR', level=0)
norwayP<-spTransform(norway,"+proj=utm +zone=32")
plot(norwayP)

# # #Redlist -----------------------------------------------------------------
# # #Make species lists with those linked to grazing changes
# # #Data from Artsdatabanken Rødlist 2015. All redlisted spp, påvirkning factor: changed habitat use
# # #Only species (not subspecies)
# # #Filter by 'beite' 
# vasc<-read.table('Rodlista2015_vascplant.csv'
#                   ,header=T,sep=";",fileEncoding="UTF-16LE")
#  
# herbvasc<-droplevels(vasc[grep('*beit*',vasc$Påvirkningsfaktorer,ignore.case=T),])
# dim(herbvasc)
# herbvasc$Vitenskapelig.navn
# 
# moss<-read.table('Rodlista2015_mosses.csv'
#                  ,header=T,sep=";",fileEncoding="UTF-16LE")
# 
# herbmoss<-droplevels(moss[grep('*beit*',moss$Påvirkningsfaktorer,ignore.case=T),])
# dim(herbmoss)
# herbmoss$Vitenskapelig.navn
# 
# lichen<-read.table('Rodlista2015_lichens.csv'
#                  ,header=T,sep=";",fileEncoding="UTF-16LE")
# 
# herblichen<-droplevels(lichen[grep('*beit*',lichen$Påvirkningsfaktorer,ignore.case=T),])
# dim(herblichen)
# herblichen$Vitenskapelig.navn
# 
# 
# redlistsp_all<-rbind(vasc,moss,lichen)
# redlistsp_all$PlantGroup<-c(rep('Vascular',times=nrow(vasc)),rep('Bryphyte',times=nrow(moss)),rep('Lichen',times=nrow(lichen)))
# 
# View(redlistsp_all)
# 
# 
# # GBIF import -------------------------------------------------------------
# 
# #Keyed download with doi linkk
# keysL<-sapply(as.character(herblichen$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)
# paste(keysL,collapse=',')#Copy and paste to below
# odLichen<-occ_download('taxonKey = 2609133,2607387,3408876,2603661,2602561,2609349,3422592,3433876,2601243,2608140,5260765,7083298,2605872,7425899,5260770,9198278,5258394,8704544,3397573,3429282,3390552,2609409,2609180,5260761,2607725,8335421,5260565,2599762,3429909,3389602'
#                        ,'country = NO','hasCoordinate = TRUE',
#                        user='jamesspeed',pwd='*****',email='*****')
# occ_download_meta(odLichen)
# gbif_citation(occ_download_meta(odLichen))# GBIF Occurrence Download https://doi.org/10.15468/dl.pl144i Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"
# 
# keysM<-sapply(as.character(herbmoss$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)
# paste(keysM,collapse=',')
# odMoss<-occ_download('taxonKey = 5280585,4276928,5283214,7800507,2689413'
#                        ,'country = NO','hasCoordinate = TRUE',
#                        user='jamesspeed',pwd='*****',email='*****')
# occ_download_meta(odMoss)
# gbif_citation(occ_download_meta(odMoss))# "GBIF Occurrence Download https://doi.org/10.15468/dl.ydfbtt Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"
# 
# keysV<-sapply(as.character(herbvasc$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)
# paste(keysV,collapse=',')
# odVasc<-occ_download('taxonKey = 9139754,8558004,2704845,3025813,7931051,2792588,9020552,5405976,2849252,8231647,3001509,8917443,2820517,3012376,2975152,8277403,3033129,5361866,5284517,8869754,7270427,5410857,8915737,2926086,5409958,2787993,3111049,9177060,2927078,3012509,2926055,7660935,2914396,2975380,2925944'
#                      ,'country = NO','hasCoordinate = TRUE',
#                      user='jamesspeed',pwd='*****',email='*****')
# occ_download_meta(odVasc)
# gbif_citation(occ_download_meta(odVasc))# "GBIF Occurrence Download https://doi.org/10.15468/dl.7eqijw Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"
# 
# odLichdat<-occ_download_get(odLichen,overwrite=T)
# odMossdat<-occ_download_get(odMoss,overwrite=T)
# odVascdat<-occ_download_get(odVasc,overwrite=T)
# 
# odLichendata<-occ_download_import(odLichdat)
# odMossdata<-occ_download_import(odMossdat)
# odVascdata<-occ_download_import(odVascdat)
# 
# summary(as.factor(odVascdata$species))
# summary(as.factor(odMossdata$species))
# summary(as.factor(odLichendata$species))
# 
# #Convert to SPDF
# mossspdf<-SpatialPointsDataFrame(cbind(odMossdata$decimalLongitude,odMossdata$decimalLatitude),as.data.frame(odMossdata),proj4string = crs(norway))
# lichenspdf<-SpatialPointsDataFrame(cbind(odLichendata$decimalLongitude,odLichendata$decimalLatitude),as.data.frame(odLichendata),proj4string = crs(norway))
# vascspdf<-SpatialPointsDataFrame(cbind(odVascdata$decimalLongitude,odVascdata$decimalLatitude),as.data.frame(odVascdata),proj4string = crs(norway))
# 
# moss_utm<-spTransform(mossspdf,crs(norwayP))
# lichen_utm<-spTransform(lichenspdf,crs(norwayP))
# vasc_utm<-spTransform(vascspdf,crs(norwayP))
# 
# plot(norwayP)
# points(vasc_utm,pch=16,cex=0.1)
# 
# 
# # Clipping to Norway outline bufferd --------------------------------------
# 
# #Clip to remove points outside of Norway
# #Polygon of 1km buffer around Norway
# library(rgeos)
# norway1km<-gBuffer(norwayP, width=1000)
# #plot(norway1km)
# #plot(norwayP,add=T)
# 
# #Clip
# lichen_clip<-lichen_utm[which(!is.na(over(lichen_utm,norway1km))),]
# moss_clip<-moss_utm[which(!is.na(over(moss_utm,norway1km))),]
# vasc_clip<-vasc_utm[which(!is.na(over(vasc_utm,norway1km))),]
# plot(norwayP)
# points(lichen_utm,pch=16,col=2)
# points(lichen_clip,pch=16,col=3)
# plot(norwayP)
# points(moss_utm,pch=16,col=2)
# points(moss_clip,pch=16,col=3)
# plot(norwayP)
# points(vasc_utm,pch=16,col=2)
# points(vasc_clip,pch=16,col=3)
# 
# AllForestRedList<-rbind(lichen_clip,moss_clip,vasc_clip)
# AllForestRedList$PlantGroup<-c(rep('Lichen',times=nrow(lichen_clip)),rep('Bryophyte',times=nrow(moss_clip)),rep('Vascular',times=nrow(vasc_clip)))
# plot(norwayP)
# points(AllForestRedList[AllForestRedList$PlantGroup=='Lichen',])
# 
# #Merge with red list details
# match(AllForestRedList$species,redlistsp_all$Vitenskapelig.navn)
# ForestRedList_adb<-merge(AllForestRedList,redlistsp_all,by.x='species',by.y='Vitenskapelig.navn',all.x=T)
# 
# 
# # Final spp data ----------------------------------------------------------
# 
# write.table(lichen_clip,'Lichens_forest_redlisted_herbivory.csv')
# write.table(moss_clip,'Moss_forest_redlisted_herbivory.csv')
# write.table(vasc_clip,'Vascular_forest_redlisted_herbivory.csv')
# 
# fwrite(ForestRedList_adb@data,file='RedListedForestSpeciesNorwayBeite.csv')


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
names(PredVars)[20:25]<-c('Elevation','Land_Cover','Forest_Type','Forest_Productivity','Vegetation_Type','SoilpH')


#Convert to factor variables
#Ratify land cover
PredVars$Land_Cover<-ratify(PredVars$Land_Cover)
ratlc<- levels(PredVars$Land_Cover)[[1]]
ratlc[["Land_Cover"]] <- c("Built-up","Agricultural","Forest","Open-natural vegetation","Mires","Glaciers/Ice/Snow","Freshwater","Sea","NA")
levels(PredVars$Land_Cover)<-ratlc
levelplot(PredVars$Land_Cover)

#Ratify forest productivty
PredVars$Forest_Productivity[PredVars$Forest_Productivity>18]<-NA #Class 99 ikke registrert. Gjelder skogområder som ligger utenfor AR5 kartleggingsområder
PredVars$Forest_Productivity<-ratify(PredVars$Forest_Productivity)
ratlcp<-levels(PredVars$Forest_Productivity)[[1]]
ratlcp[['Forest_Productivity']]<-c('Unproductive','Low','Medium','High')
levels(PredVars$Forest_Productivity)<-ratlcp
levelplot(PredVars$Forest_Productivity)

#Ratify forest type
PredVars$Forest_Type[PredVars$Forest_Type>33]<-NA #	99 Ikke registrert. Gjelder skogområder som ligger utenfor AR5 kartleggingsområder
PredVars$Forest_Type<-ratify(PredVars$Forest_Type)
ratlct<-levels(PredVars$Forest_Type)[[1]]
ratlct
ratlct[['ForestType']]<-c('Coniferous','Deciduous','Mixed')
levels(PredVars$Forest_Type)<-ratlct
levelplot(PredVars$Forest_Type)



#Correlation plots
pairs(PredVars)


#Reimport Species data
redlistforest<-fread('RedListedForestSpeciesNorwayBeite.csv')

#Convert to spdf
rlforsp<-SpatialPointsDataFrame(cbind(redlistforest$decimalLongitude,redlistforest$decimalLatitude),redlistforest,proj4string = crs(norway))
rlforsp_utm<-spTransform(rlforsp,crs(norwayP))

plot(norwayP)
points(rlforsp_utm[rlforsp_utm$PlantGroup.x=='Vascular',],col='green',pch=16,cex=0.2)
points(rlforsp_utm[rlforsp_utm$PlantGroup.x=='Lichen',],col='brown',pch=16,cex=0.2)
points(rlforsp_utm[rlforsp_utm$PlantGroup.x=='Bryophyte',],col='blue',pch=16,cex=0.2)

#Select points only with good geographic precision (coordinateuncertainty <1415m (sqrt(1000^2+1000^2)))
rlfor_use<-rlforsp_utm[!is.na(rlforsp_utm$coordinateUncertaintyInMeters) & rlforsp_utm$coordinateUncertaintyInMeters<=1415,]

#Select points Only in forest land-cover
extforest<-extract(PredVars$Land_Cover,rlfor_use)
extforest[is.na(extforest)]<-0
forestonly<-rlfor_use[extforest==30,]

#Points per species
spforocc<-with(forestonly@data,tapply(species,species,length))
spforocc
hist(spforocc)

#Points per higher taxa
tapply(forestonly$PlantGroup.x,forestonly$PlantGroup.x,length)

#Plot all species
col<-colorRampPalette('white')
#Make a raster to plot
noralt<-getData('alt',country='NOR')
nornull<-crop(projectRaster(noralt,crs=crs(norwayP)),norwayP)
values(nornull)<-0
levelplot(nornull,margin=F,colorkey=F
          ,col.regions=col,main='Red listed forest species', xlab=NULL, ylab=NULL, scales=list(draw=FALSE))+
  layer(sp.polygons(norwayP))+
  layer(sp.points(forestonly[forestonly$PlantGroup.x=='Vascular',],pch=16,cex=0.2,col='green'))+
  layer(sp.points(forestonly[forestonly$PlantGroup.x=='Lichen',],pch=16,cex=0.2,col='blue'))+
  layer(sp.points(forestonly[forestonly$PlantGroup.x=='Bryophyte',],pch=16,cex=0.5,col='tan4'))


# Plot species ------------------------------------------------------------
#Lichens
lichen_use<-forestonly[forestonly$PlantGroup.x=='Lichen',]
p<-list()
for(i in 1:length(levels(as.factor(lichen_use$species)))){
  p[[i]]<-levelplot(nornull,margin=F,colorkey=F,col.regions=col, xlab=NULL, ylab=NULL, scales=list(draw=FALSE)
                    ,main=list(label=levels(as.factor(lichen_use$species))[i],cex=0.7,fontface=3))+
    layer(sp.polygons(norwayP))+
    layer(sp.points(lichen_use[lichen_use$species==levels(as.factor(lichen_use$species))[i],],pch=16,cex=0.5))
}

#Make locations to plot (28 lichen species)
rowi<-c(rep(1:5,times=6))
coli<-c(rep(1:6,each=5))
tiff(width=210,height=297,units='mm',res=300,'Lichen_distributions.tif')
lattice.options(  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
                  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0)))
for (i in 1:length(p)){
  print(p[[i]],split=c(rowi[i],coli[i],5,6),more=T)}
dev.off()

#Mosses
moss_use<-forestonly[forestonly$PlantGroup.x=='Bryophyte',]
m<-list()
for(i in 1:length(levels(as.factor(moss_use$species)))){
  m[[i]]<-levelplot(nornull,margin=F,colorkey=F,col.regions=col, xlab=NULL, ylab=NULL, scales=list(draw=FALSE)
                    ,main=list(label=levels(as.factor(moss_use$species))[i],cex=0.7,fontface=3))+
    layer(sp.polygons(norwayP))+
    layer(sp.points(moss_use[moss_use$species==levels(as.factor(moss_use$species))[i],],pch=16,cex=0.5))
}

#Make locations to plot (5 moss species)
rowi<-c(rep(1:5,times=1))
coli<-c(rep(1,each=5))
tiff(width=210,height=60,units='mm',res=300,'Bryophyte_distributions.tif')
lattice.options(  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
                  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0)))
for (i in 1:length(m)){
  print(m[[i]],split=c(rowi[i],1,5,1),more=T)}
dev.off()


#Vascular
vasc_use<-forestonly[forestonly$PlantGroup.x=='Vascular',]
v<-list()
for(i in 1:length(levels(as.factor(vasc_use$species)))){
  v[[i]]<-levelplot(nornull,margin=F,colorkey=F,col.regions=col, xlab=NULL, ylab=NULL, scales=list(draw=FALSE)
                    ,main=list(label=levels(as.factor(vasc_use$species))[i],cex=0.7,fontface=3))+
    layer(sp.polygons(norwayP))+
    layer(sp.points(vasc_use[vasc_use$species==levels(as.factor(vasc_use$species))[i],],pch=16,cex=0.5))
}

#Make locations to plot (31 vasc species)
rowi<-c(rep(1:5,times=7))
coli<-c(rep(1:7,each=5))
tiff(width=210,height=297,units='mm',res=150,'Vascular distributions.tif')
lattice.options(  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
                  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0)))
for (i in 1:length(v)){
  print(v[[i]],split=c(rowi[i],coli[i],5,7),more=T)}
dev.off()


# Distribution modelling --------------------------------------------------

#Background data
background<-read.table('BackgroundBiasCorrected.txt',header=T)#Vascular plants bias file used for Speed & Austrheim et al. 2017
backsp<-SpatialPoints(background,proj4string = crs(norwayP))
plot(norwayP)
points(backsp,cex=0.1,pch=16)

#MaxEnt Tuning
tuneparameters<-ENMevaluate(occ=forestonly@coords[forestonly$species==levels(as.factor(forestonly$species))[[1]],],
                            env=PredVars[[c(10,12,15,33,41,49,22,23)]],
                            categoricals=c("Forest_Type","Forest_Productivity"),
                            method="block",
                            bg.coords=backsp)
tuneparameters@results[which.min(tuneparameters@results$AICc),]

#MaxEnt modelling
me1<-maxent(p=forestonly@coords[forestonly$species==levels(as.factor(forestonly$species))[[1]],],
            x=PredVars[[c(10,12,15,33,41,49,22,23)]],
            factors=c("Forest_Type","Forest_Productivity"),
            a=backsp,
            args=c('betamultiplier=2.0','threshold=TRUE','product=TRUE',"-P","-J"))
me1
