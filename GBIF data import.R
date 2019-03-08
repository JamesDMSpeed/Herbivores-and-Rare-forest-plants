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

PredVars<-stack('PredictorVariables')#Available here https://ntnu.box.com/s/wcmr0dgoyz2yu6ielw6er1pm7h0gaisa 
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

#Mask out points outside of Norway coverage (worldclim data)
rlfor_useNor<-rlfor_use[!is.na(extract(PredVars[[10]],rlfor_use)),]

#Points per species
spNorwayGoodPrec<-with(rlfor_useNor@data,tapply(species,species,length))

#Select points Only in forest land-cover
extforest<-extract(PredVars$Land_Cover,rlfor_useNor)
extforest[is.na(extforest)]<-0
forestonly<-rlfor_use[extforest==30,]

#Points per species
spforocc<-with(forestonly@data,tapply(species,species,length))
spforocc
hist(spforocc)

#Only in forest with productivity
extforestprod<-extract(PredVars$Forest_Productivity,forestonly)
extforestprod[is.na(extforestprod)]<-0
forestprodonly<-forestonly[extforestprod>0,]
#Only in forest with forest type too
extforestprodtype<-extract(PredVars$Forest_Type,forestprodonly)
extforestprodtype[is.na(extforestprodtype)]<-0
forestprodonlytype<-forestprodonly[extforestprodtype>0,]

#Points per species - prod
spforprodocc<-with(forestprodonly@data,tapply(species,species,length))
spforprodocc
hist(spforprodocc)

#Points per species - prod and type
spforprodtypeocc<-with(forestprodonlytype@data,tapply(species,species,length))
spforprodtypeocc
hist(spforprodtypeocc)

#Counts of all poitns per species in different type
reccts<-as.data.frame(t(Reduce(function(...) merge(..., all=TRUE), list(t(as.data.frame(spNorwayGoodPrec)),
                                                t(as.data.frame(spforocc)),
                                                t(as.data.frame(spforprodocc)),
                                                t(as.data.frame(spforprodtypeocc))))))
                      
colnames(reccts)<-c('Records with forest type and productivity','Records with forest productivity','Records in forest','Records with good precisicion within worldclim')
reccts$Group<-rlfor_use$PlantGroup.x[match(rownames(reccts),rlfor_use$species)]
write.csv(reccts[,5:1],'RecordCountsperSpecies.csv')


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
#Currently only points with forest productivity and type data 


#Background data
background<-read.table('BackgroundBiasCorrected.txt',header=T)#Vascular plants bias file used for Speed & Austrheim et al. 2017
backsp<-SpatialPoints(background,proj4string = crs(norwayP))
plot(norwayP)
points(backsp,cex=0.1,pch=16)

#Background in forest
backforestext<-extract(PredVars$Land_Cover,backsp)
backforestext[is.na(backforestext)]<-0
backforest<-backsp[backforestext==30,]
  
#MaxEnt Tuning ####
#Always allow linear features

 args_sp<-vector("list", length(levels(as.factor(forestprodonlytype$species))))
 names(args_sp)<-levels(as.factor(forestprodonlytype$species))
 nrecs<-summary(as.factor(forestprodonlytype$species))

 # for(i in 1:length(levels(as.factor(forestprodonlytype$species)))){
#   print(i)
#   print(levels(as.factor(forestprodonlytype$species))[[i]])
#   print(nrecs[i])  
#   if(nrecs[i]>=15){
#   tuneparameters<-ENMevaluate(occ=forestprodonlytype@coords[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[[i]],],
#                             env=PredVars[[c(10,12,15,33,41,49,22:23)]],#Many NA values for forest type and productivity
#                             RMvalues = c(0.5,1,1.5,2,2.5,3,3.5,4,6,8), 
#                             fc = c("L", "LQ","LQH", "LQHP", "LQHPT"),
#                             categoricals=c("Forest_Type","Forest_Productivity"),
#                             method="block",
#                             bg.coords=backforest)
# tuneparameters@results[which.min(tuneparameters@results$AICc),]
# 
# b<-tuneparameters@results$rm[which.min(tuneparameters@results$AICc)]
# lin <-grepl('L',tuneparameters@results$features[which.min(tuneparameters@results$AICc)])
# quad<-grepl('Q',tuneparameters@results$features[which.min(tuneparameters@results$AICc)])
# prod<-grepl('P',tuneparameters@results$features[which.min(tuneparameters@results$AICc)])
# hing<-grepl('H',tuneparameters@results$features[which.min(tuneparameters@results$AICc)]) 
# thres<-grepl('T',tuneparameters@results$features[which.min(tuneparameters@results$AICc)])
# 
# args_sp[[i]]<-as.character(c(paste0('betamultiplier=',b),
#                       paste0('linear=',lin),
#                       paste0('quadratic=',quad),
#                       paste0('product=',prod),
#                       paste0('hinge=',hing),
#                       paste0('threshold=',thres),
#                       "-P",
#                       "-J",
#                       'replicates=5'
# ))}
# }

#Write tuning as list
#saveRDS(args_sp,file='Tuning')
args_sp<-readRDS('Tuning')

#MaxEnt modelling ####
me_list<-list()
#for(i in 1:4){
for(i in 1:length(levels(as.factor(forestprodonlytype$species)))){
  print(i)
  print(levels(as.factor(forestprodonlytype$species))[[i]])
  print(args_sp[[i]])
  #nrecs<-nrow(forestprodonlytype@coords[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[[i]],])
  print(nrecs[i])  
  if(nrecs[i]>=25){
  me_list[[i]]<-maxent(p=forestprodonlytype@coords[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[[i]],],
            x=PredVars[[c(10,12,15,33,41,49,22,23)]],#Many NA values for forest type and productivity
            factors=c("Forest_Type","Forest_Productivity"),
            a=backforest,
            args=args_sp[[i]],
            path=paste0('MaxEnt/',levels(as.factor(forestprodonlytype$species))[[i]]))
}}
names(me_list)<-levels(as.factor(forestprodonlytype$species))[nrecs>=25]


#Predictions ####
predictions<-list()
#for(i in 1:4){
for(i in 1:length(levels(as.factor(forestprodonlytype$species)))){
  #nrecs<-nrow(forestprodonlytype@coords[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[[i]],])
  print(nrecs[i])  
  if(nrecs[i]>=25){ predictions[[i]]<-predict(me_list[[i]],PredVars)
} }
predictions_sp<-predictions[nrecs>=25]
names(predictions_sp)<-levels(as.factor(forestprodonlytype$species))[nrecs>25]

#Averaging predictions across multiple runs
meanpreds<-lapply(predictions_sp,function(x)calc(x,mean))

allspeciespredictions<-stack(meanpreds)
writeRaster(allspeciespredictions,'ModelPredictions/.tif',format='GTiff',bylayer=T,suffix=names(allspeciespredictions))
writeRaster(allspeciespredictions,'ModelPredictions/AllSpeciesPredictions.tif',format='GTiff')




# Using sdm package ------------------------------------------------------
library(sdm)
#https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.01881

sem<-function(x)sd(x,na.rm=T)/length(!is.na(x))

#InstallAll()#One time to install all dependent packages
#Testing

ulmgla<-forestprodonlytype[forestprodonlytype@data$species=='Ulmus glabra',]#@coords
names(ulmgla)[1]<-'ulmgla'
ulmgla<-cbind(ulmgla[,1],ulmgla@coords)
ulmgla$ulmgla<-'ulmgla'#Problems with spaces in species names...

#ulmgla<-data.frame(cbind(ulmgla=rep('ulmgla',times=nrow(ulmgla)),ulmgla))
sdmdataset<-sdmData(ulmgla~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16,train=ulmgla,
                    predictors=PredVars,bg=list(n=1000,method='gRandom',remove=TRUE))
sdmdataset
plot(sdmdataset)

mod1<-sdm(ulmgla~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16,data=sdmdataset,
          methods=c('glm','gam','gbm','cart','fda','rf'),
          replication=c('cv','boot'),cv.folds=5)

mod1<-sdm(ulmgla~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16,data=sdmdataset,
          methods=c('gbm','tree','mda','fda'),replication=c('cv','boot'),cv.folds=5,n=10)
roc(mod1)
rcurve(mod1)
getVarImp(mod1,1)# 1 model at a time


# Making good dataframe ---------------------------------------------------

listdf<-list()
for (i in 1: length(levels(as.factor(forestprodonlytype$species)))){
#  for (i in 1:3){
  a<-forestprodonlytype[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[i],]
  b<-cbind(a[,1],a@coords)
 # names(b)[1]<-levels(as.factor(forestprodonlytype$species))[i]
  names(b)[1]<-'species'
  listdf[[i]]<-b
  }

  
sdmdataset<-sdmData(ulmgla~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16,train=ulmgla,
                      predictors=PredVars,bg=list(n=1000,method='gRandom',remove=TRUE))
  
s1<-sdmData(Ajuga_reptans~roe_deer2015,train = listdf[[1]],predictors = PredVars,
            bg=list(n=1000,method='gRandom',remove=TRUE))

slist<-list()
for(i in 1:3){
  #s[[i]]<-sdmData(paste(levels(as.factor(forestprodonlytype$species))[i])~
  s[[1]]<-sdmData(species~
                    roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16,
    train=listdf[[1]],
    predictors=PredVars,bg=list(n=1000,method='gRandom',remove=TRUE))
}


AllSpp <- do.call("rbind", listdf)
#Remove space from species name to avoid errors
AllSpp$species<-sub(" ","_",AllSpp$species)

#Make the sdm dataset with all species and relevent environmental variables (specifiy factors)
#1000 random points in the forest area
bg<-sampleRandom(PredVars$Forest_Productivity,1000,sp=T)

sdmdataset<-sdmData(species~roe_deer2015+red_deer2015+moose2015
                    +bio10_16+bio12_16+bio15_16+SoilpH
                    +f(Forest_Type)+f(Forest_Productivity)
                    ,train=AllSpp,predictors=PredVars,bg=list(n=1000,method='gRandom',remove=TRUE))
sdmdataset

#Model with all species concurrently 
#All species with >=20 records in forest
sdm_Allspp<-sdm( Ajuga_reptans               +Anastrophyllum_donnianum   
                +Arnica_montana              +Asperugo_procumbens          
                +Campanula_barbata           +Campanula_cervicaria        +Cetrelia_olivetorum        
                +Cinna_latifolia             +Collema_curtisporum        
                +Collema_occultatum          +Crepis_praemorsa            +Cypripedium_calceolus      
                +Dactylorhiza_sambucina      +Epipogium_aphyllum          +Galium_sterneri            
                +Gentianella_campestris      +Gyalecta_flotowii           +Gyalecta_truncigena        
                +Gyalecta_ulmi               +Hackelia_deflexa            +Herbertus_stramineus        +Heterodermia_speciosa      
                +Lithospermum_officinale     +Malus_sylvestris            +Menegazzia_subsimilis       +Menegazzia_terebrata       
                +Opegrapha_vermicellifera    +Ophrys_insectifera          +Pectenia_cyanoloma         
                +Phaeophyscia_kairamoi       +Physconia_detersa           +Pseudorchis_albida          +Ramalina_dilacerata        
                +Ramalina_sinensis           +Ramboldia_subcinnabarina    +Rinodina_disjuncta          
                +Schismatomma_graphidioides  +Scorzonera_humilis          +Sorbus_lancifolia           +Sorbus_subpinnata           +Staurolemma_omphalarioides  +Taxus_baccata              
                +Thalictrum_minus            +Thalictrum_simplex          +Thelotrema_macrosporum                  
                +Ulmus_glabra                +Vicia_cassubica                     
          
                ~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16+Forest_Type+Forest_Productivity+SoilpH,
                data=sdmdataset,
          methods=c('glm','gam','rf','gbm','mda','fda','brt'),
          replication=c('cv'),cv.folds=5)

saveRDS(sdm_Allspp,'SDM package/SDMAllSpecies')

#Extract model evaluations
#modeval<-cbind(sdm_Allspp@run.info,getEvaluation(sdm_Allspp))
modeval<-merge(sdm_Allspp@run.info,getEvaluation(sdm_Allspp),by='modelID')
write.csv(modeval,'SDM package/AllModels_Evaluation')

with(modeval,tapply(AUC,list(species,method),mean))
with(modeval,tapply(AUC,list(species,method),sd))

#Extract variable importances
varimplist<-list()
for (i in 1:max(sdm_Allspp@run.info$modelID)){
  ifelse(sdm_Allspp@run.info$success[i]==TRUE,
         {varimplist[[i]]<-getVarImp(sdm_Allspp,id=i)@varImportance
         varimplist[[i]]$species<-sdm_Allspp@run.info$species[i]
         varimplist[[i]]$method<-sdm_Allspp@run.info$method[i]
         varimplist[[i]]$repid<-sdm_Allspp@run.info$replicationID[i]}
         ,print(paste('Model failiure run ',i)))
}
AllVarImp<-do.call('rbind',varimplist)
write.csv(AllVarImp,'SDM package/AllModelsVariableImportance')
#Plot
varimpmean<-with(AllVarImp,tapply(corTest,list(variables,species),mean))
varimpsem<-with(AllVarImp,tapply(corTest,list(variables,species),sem))
par(mar=c(5,12,1,1))
b1<-barplot(varimpmean,beside=T,horiz=T,las=1,legend.text=T)
arrows(varimpmean+varimpsem,b1,varimpmean-varimpsem,b1,code=3,angle=90,length=0.05)

#Response curves
responsecurvelist<-list()
for (i in 1:length(levels(as.factor(sdm_Allspp@run.info$species)))){
responsecurvelist[[i]]<-rcurve(sdm_Allspp,id=sdm_Allspp@run.info$modelID[sdm_Allspp@run.info$species==levels(as.factor(sdm_Allspp@run.info$species))[i]]
       ,mean=T,main=levels(as.factor(sdm_Allspp@run.info$species))[i])
}
responsecurvelist[[1]]
#Ensemble models for each species
ensemblelist<-list()
for (i in 1:length(levels(as.factor(sdm_Allspp@run.info$species)))){
  ensemblelist[[i]]<-ensemble(sdm_Allspp,newdata=PredVars,filename=paste('EnsemblePredictions/',levels(as.factor(sdm_Allspp@run.info$species))[i]),
                              setting=list(method='weighted',stat='AUC',id=sdm_Allspp@run.info$modelID[sdm_Allspp@run.info$species==levels(as.factor(sdm_Allspp@run.info$species))[i]]))
}

#Niches 
niche(PredVars,ensemblelist[[1]],n=c('bio16_16','moose2015'))