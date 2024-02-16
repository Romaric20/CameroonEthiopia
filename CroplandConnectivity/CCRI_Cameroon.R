
#devtools::install_github("GarrettLab/CroplandConnectivity", subdir = "geohabnet")
library(igraph)
library(sp)
library(maps)
library(geosphere)
library(RColorBrewer)
library("colorspace") 
data("countriesLow")
library(terra)
library(geohabnet)


#use geohabnet function to download potato raster

banana<-cropharvest_rast("banana","mapspam2017Africa") #may take a few minutes
plantain<-cropharvest_rast("plantain","mapspam2017Africa")

banplan <- mean(banana, plantain)

cassava<-cropharvest_rast("cassava","mapspam2017Africa")
potato<-cropharvest_rast("potato","mapspam2017Africa")
sweetpotato<-cropharvest_rast("sweet potato","mapspam2017Africa")



#----------- East hemisphere-------------------
latifrom <- 1 #latitude: 
latito <- 16
longifrom <- 8#longitude: 
longito<- 18
CMR_ext <- extent(8, 18, 1, 16)


Cameroon<- raster::getData("GADM",
                           
                           country = "Cameroon",
                           
                           level = 1)
#----------------------------------------------------------
#------------ Set palette

palette1 <- c( "#F4E156FF", "#F6D746FF", "#F8CD37FF", "#FAC329FF", "#FBB91EFF", "#FCAF13FF", 
               "#FCA50BFF", "#FB9C06FF", "#FA9207FF", "#F8890CFF", "#F68013FF", "#F37819FF",
               "#F06F20FF", "#EC6727FF", "#E85F2EFF", "#E25834FF", "#DD5139FF", "#D74B3FFF",
               "#D04545FF", "#CA404AFF", "#C33B4FFF", "#BC3754FF", "#B43359FF", "#AC305EFF",
               "#A42C60FF", "#9B2964FF", "#932667FF", "#922568FF", "#902568FF", "#8F2469FF",
               "#8D2369FF", "#8C2369FF", "#8A226AFF", "#88226AFF", "#87216BFF", "#85216BFF",
               "#84206BFF", "#82206CFF", "#801F6CFF", "#7F1E6CFF", "#7D1E6DFF", "#7C1D6DFF",
               "#7A1D6DFF", "#781C6DFF", "#771C6DFF", "#751B6EFF", "#741A6EFF", "#721A6EFF",
               "#71196EFF", "#6E196EFF", "#6D186EFF", "#6B186EFF", "#6A176EFF", "#68166EFF",
               "#66166EFF", "#65156EFF", "#63156EFF", "#61136EFF", "#60136EFF", "#5E126EFF",
               "#5C126EFF", "#5B126EFF", "#59106EFF", "#58106EFF", "#560F6DFF", "#540F6DFF",
               "#530E6DFF", "#510E6CFF", "#500D6CFF", "#4D0D6CFF", "#4C0C6BFF", "#4A0C6BFF",
               "#490B6AFF", "#470B6AFF", "#450A69FF", "#440A68FF", "#420A68FF", "#400A67FF",
               "#3E0966FF", "#3D0965FF", "#3B0964FF", "#390963FF", "#380962FF", "#360961FF",
               "#340A5FFF", "#320A5EFF", "#310A5CFF", "#2F0A5BFF", "#2D0B59FF", "#2B0B57FF",
               "#290B55FF", "#280B53FF", "#250C51FF", "#240C4EFF", "#230C4BFF", "#200C49FF",
               "#1F0C47FF", "#1D0C44FF", "#1C0C42FF", "#1A0C40FF", "#190C3DFF", "#170C3BFF",
               "#150B38FF", "#150B36FF", "#130A33FF", "#110A31FF", "#11092EFF", "#0F092CFF",
               "#0D082AFF", "#0C0827FF", "#0B0725FF", "#0A0723FF", "#090620FF", "#08051EFF",
               "#07051CFF", "#060419FF", "#050418FF", "#040315FF", "#040312FF", "#030210FF",
               "#02020EFF", "#02020CFF", "#02010AFF", "#010108FF", "#010106FF", "#010005FF",
               "#000004FF", "#000004FF", "#000004FF")


## 1.2 Customize crop and values of parameters
beta0<-0.5                                       ###
beta<-1                                          ###
beta1<-1.5                                       ###
gamma00<-0.05                                    ###
gamma0<-0.1                                      ###
gamma<-0.2                                       ###
gamma1<-0.3                                      ###
gamma2<-1                                        ###                                
cutoff1<- 0.001                                  ###
cutoff2 <- 0.01 # cutoff of adjancecy matrix     ###


#Select the crophavest

cropharvest <- banplan 

#---------------------------------------------------
# aggregated resolution
Resolution <- 2 # Set aggregated resolution, for example, assign 12 for 1 degree.

#----------- total mean aggregration -----------------------------
cropharvestAGG <- aggregate(cropharvest, fact = Resolution, fun=sum, na.action = na.omit)
cropharvestAGGTM <- cropharvestAGG / Resolution / Resolution #TOTAL MEAN


#----------- crop cropland area for the west hemisphere ----------
cropharvestAGGTM_crop <- crop(cropharvestAGGTM, CMR_ext)	
plot(cropharvestAGGTM_crop, col = palette1,zlim= c(0, 2))
plot(countriesLow, add=TRUE, border = "white")
#----------- Extract cropland density data -----------------------
CropValues <- values(cropharvestAGGTM_crop)
CropValuesAzero <- which(CropValues > cutoff) # find the cells with value > 0.0001
cropValue <- CropValues[CropValuesAzero]
#----------- Extract xy corrdination for "povalue" cells ---------
lon <- NULL 
lat <- NULL 

for(i in 1:length(CropValuesAzero)){
  temp <- extentFromCells(cropharvestAGGTM_crop, CropValuesAzero[i])
  AVxminO <- temp[1]
  lon <- c(lon, AVxminO)
  AVymaxO <- temp[4]
  lat <- c(lat, AVymaxO)
}
#---------------------------------------------------------------
#---------------------------------------------------------------

cropdata1 <- data.frame(lon, lat, cropValue)

latilongimatr <- cropdata1[ ,c(1:2)]
#---- use Geosphere package, function distVincentyEllipsoid() is used to calculate the distance, defult distance is meter
dvse <- distVincentyEllipsoid(c(0,0), cbind(1, 0)) 
latilongimatr <- as.matrix(latilongimatr)
TemMat <- matrix(-999, nrow( latilongimatr),nrow(latilongimatr))

for (i in 1:nrow(latilongimatr)) {
  TemMat[i, ] <- distVincentyEllipsoid(latilongimatr[i,], latilongimatr)/dvse
}
distance_matrix <- TemMat


map_grey_background <- rast("map_grey_background.tif")

# with normalization east hemisphere
#----------------------------------------------------------
#plot cropland density

map_grey_background_CMR <- crop(map_grey_background, CMR_ext)
plot(map_grey_background_CMR, col = "grey75",  xaxt='n',  yaxt='n', axes=F, box=F, legend = F, 
     main=paste('crop density: banana-plantain'), cex.main=0.7)
plot(cropharvestAGGTM_crop, col = palette1, xaxt='n',
     yaxt='n', axes=F, box=F, add = TRUE)
plot(countriesLow, add=TRUE, border = "darkblue")
plot(Cameroon, col = NA, border = "darkblue", add=TRUE)


#-------------------------------------------------------------------------
# CCRI calculated by Inverse power-law function  

CCRI_powerlaw_function <- function(beta, cutoffadja, distance_matrix, lon, lat, cropValue, cropRaster, CellNumber)   {
  ##############################################
  #### create adjacency matrix
  
  distancematr <- distance_matrix # pairwise distance matrix
  #---- end of code
  distancematrexp <- distancematr^(-beta) #use function C=AX^(-beta), here A=1, X=distancematr
  cropmatr <- cropValue # complete gravity model with crop data
  cropmatr1 <- matrix(cropmatr, , 1 )
  cropmatr2 <- matrix(cropmatr, 1, )
  
  cropmatrix <- cropmatr1 %*% cropmatr2
  cropmatrix <- as.matrix(cropmatrix)
  cropdistancematr <- distancematrexp * cropmatrix # adjacecy matrix
  logicalmatr <- cropdistancematr > cutoffadja # adjacency matrix after threshold
  stan <- cropdistancematr * logicalmatr
  stan <- round(stan, 6) # use round() because betweenness() may have problem when do the calculation
  cropdistancematrix <- graph.adjacency(stan,mode=c("undirected"),diag=F,weighted=T)#create adjacency matrix
  ##############################################
  ## sum of nearest neighbors degree
  knnpref0<-graph.knn(cropdistancematrix,weights=NA)$knn
  knnpref0[is.na(knnpref0)]<-0
  degreematr<-degree(cropdistancematrix)
  knnpref<-knnpref0*degreematr
  if(max(knnpref)==0){knnprefp=0}else
    if(max(knnpref)>0){knnprefp=knnpref/max(knnpref)/6}
  
  ##############################################
  #### node degree, node strengh 
  ####
  nodestrength<-graph.strength(cropdistancematrix) 
  nodestrength[is.na(nodestrength)]<-0
  if(max(nodestrength)==0){nodestr=0}else
    if(max(nodestrength)>0){nodestr=nodestrength/max(nodestrength)/6}
  ##############################################
  #### betweenness centrality
  #### 
  between<-betweenness(cropdistancematrix, weights = (1.0001*max(E(cropdistancematrix)$weight)) - E(cropdistancematrix)$weight)
  between[is.na(between)]<-0
  if(max(between)==0){betweenp=0}else
    if(max(between)>0){betweenp=between/max(between)/2}
  ##############################################
  #### eigenvector and eigenvalues
  #### 
  eigenvectorvalues<-evcent(cropdistancematrix)
  ev<-eigenvectorvalues$vector
  ev[is.na(ev)]<-0
  if(max(ev)==0){evp=0}else
    if(max(ev)!=0){evp=ev/max(ev)/6}
  ##############################################
  #### CCRI is a weighted mean of 4 network metric
  ####    
  index<-knnprefp+evp+betweenp+nodestr
  
  indexpre<-cropRaster
  indexpre[]<-0
  indexpre[CellNumber]<- index
  indexv<-indexpre
  return(indexv)
}

#---------------------------------------------------------------------------
# CCRI calculated by negative exponential function 

CCRI_negExponential_function <-function(gamma,cutoffadja, distance_matrix, lon, lat, cropValue, cropRaster, CellNumber)   {
  ##############################################
  #### create adjacency matrix
  ####
  distancematr <- distance_matrix
  #---- end of code
  
  eulernumber<-exp(1)
  distancematrexponential <- eulernumber ^ (-gamma * distancematr)# exponential model
  cropmatr <- cropValue # complete gravity model with crop data
  cropmatr1 <- matrix(cropmatr,,1) # complete gravity model with crop data
  cropmatr2 <- matrix(cropmatr,1,)
  cropmatrix <- cropmatr1 %*% cropmatr2
  cropmatrix <- as.matrix(cropmatrix)
  cropdistancematr <- distancematrexponential * cropmatrix
  logicalmatr <- cropdistancematr > cutoffadja
  stan <- cropdistancematr * logicalmatr
  stan <- round(stan, 6) # use round() because betweenness() may have problem when do the calculation
  cropdistancematrix<-graph.adjacency(stan,mode=c("undirected"),diag=F,weighted=T)#create adjacency matrix
  ##############################################
  #### create network for all the selected nodes
  ####
  #V(cropdistancematrix)$color=colororder
  V(cropdistancematrix)$label.cex=0.7
  edgeweight<-E(cropdistancematrix)$weight*4000
  E(cropdistancematrix)$color="red"
  
  knnpref0<-graph.knn(cropdistancematrix,weights=NA)$knn
  knnpref0[is.na(knnpref0)]<-0
  degreematr<-degree(cropdistancematrix)
  knnpref<-knnpref0*degreematr
  if(max(knnpref)==0){knnprefp=0}else
    if(max(knnpref)>0){knnprefp=knnpref/max(knnpref)/6}
  
  ##############################################
  #### node degree, node strengh 
  ####
  nodestrength<-graph.strength(cropdistancematrix) 
  nodestrength[is.na(nodestrength)]<-0
  if(max(nodestrength)==0){nodestr=0}else
    if(max(nodestrength)>0){nodestr=nodestrength/max(nodestrength)/6}
  ##############################################
  #### betweenness centrality
  #### 
  between<-betweenness(cropdistancematrix, weights = (1.0001*max(E(cropdistancematrix)$weight)) - E(cropdistancematrix)$weight)
  between[is.na(between)]<-0
  if(max(between)==0){betweenp=0}else
    if(max(between)>0){betweenp=between/max(between)/2}
  ##############################################
  #### eigenvector and eigenvalues
  #### 
  eigenvectorvalues<-evcent(cropdistancematrix)
  ev<-eigenvectorvalues$vector
  ev[is.na(ev)]<-0
  if(max(ev)==0){evp=0}else
    if(max(ev)!=0){evp=ev/max(ev)/6}
  ##############################################
  #### plot index layer
  ####    
  index<-knnprefp+evp+betweenp+nodestr
  
  indexpre<-cropRaster
  indexpre[]<-0
  indexpre[CellNumber] <- index
  indexv<-indexpre
  return(indexv)
  
}

#--------------------------------------------------------
## sensitivity analysis CCRI BY Inverse power-law function and negative exponential 

index1 <- CCRI_powerlaw_function(beta0, cutoffadja, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index2 <- CCRI_powerlaw_function(beta, cutoffadja, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index3 <- CCRI_powerlaw_function(beta1, cutoffadja, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index4 <- CCRI_negExponential_function(gamma00, cutoffadja, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index5 <- CCRI_negExponential_function(gamma0, cutoffadja, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index6 <- CCRI_negExponential_function(gamma, cutoffadja, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index7 <- CCRI_negExponential_function(gamma1, cutoffadja, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index8 <- CCRI_negExponential_function(gamma2, cutoffadja, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

#----------------------------------------------------------------
# Complete sensitivity analysis of CCRI 

mean_index_raster <- sum (index1, index2, index3, index4, index5, index6, index7, index8) / 8
plot(mean_index_raster,col=c("white",palette1))

mean_index_raster_val <- values(mean_index_raster)
zeroId <- which(mean_index_raster_val == 0)
mean_index_raster[zeroId] <- NaN

#--- remove pixels outside of boundary
ZeroRaster <- rast("ZeroRaster.tif")
CMR_Zero <- crop(ZeroRaster, CMR_ext)
mean_index_raster <- disagg(mean_index_raster, fact = c(Resolution, Resolution))
mean_index_raster_CMR <- mean_index_raster + CMR_Zero

#---------------------------------------------------------


map_grey_background_CMR <- crop(map_grey_background, CMR_ext)

plot(map_grey_background_CMR, col = "grey75",  xaxt='n',  yaxt='n', axes=F, box=F, legend = F, 
     main=paste('Mean in cropland connectivity risk index from sensitivity analysis: banana and plantain'), cex.main=0.7)

plot(mean_index_raster_CMR , col = palette1, xaxt='n',
     yaxt='n', axes=F, box=F, add = TRUE)
plot(countriesLow, add=TRUE, border = "darkblue")
plot(Cameroon, col = NA, border = "darkblue", add=TRUE)

