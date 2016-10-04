# Implementation of PRA algorithm in R.
# Jochen Veitinger, 20.11.2015

# forest module: if you wish to use a forest mask (ASCII file where forest = 1, no forest = 0), please uncomment the respective lines at the end of the # script



 

#### Enter input parameters ####

inputRas = "nameofDTM.asc"
outPRA = "nameofoutputPRA.asc"
HS = 2.3 				# snow depth in [m]
smooth = "Regular" 			# type Regular or Smooth to define degree of smoothing (see manual)
wind = 180 				# wind direction: N= 0, W =90, S=180
windTol = 30 				# wind tolerance (see manaual) 
work_dir = "H:/mydocuments"		# working directory to save intermediate results


 

#### Using R libraries ####
print("Loading Libraries....")

library(sp)
library(gstat)
library(shapefiles)
library(foreign)
library(methods)

library(plyr)
#library(rgdal)   
library(raster)
library(RSAGA)
#library(maptools)

####  preliminary calculations ####

myenv=rsaga.env(workspace=work_dir,
                path="H:/Tool/SAGA-GIS",
                modules="H:/Tool/SAGA-GIS/modules") 

setwd(work_dir)

print("Begin Calculations....")


asc<- raster(inputRas)
asc.extent  <- extent(asc)
head <- read.ascii.grid.header(inputRas)
writeRaster(asc, "outputRas", format="SAGA", overwrite=TRUE)


if (smooth == "Regular"){
  cv = 0.35
} else{
  cv = 0.2
}

####  experimental function to relate snow depth with scale ####


i_max = ceiling((3*HS * cv) - 0.1)


#### Computation of windshelter index ####


# DTM for windshelter calculation


windRaster<-aggregate(asc, fact= 5,expand=TRUE)


if (i_max > 2) { 
  windRaster_new<-aggregate(asc, fact= 2*i_max,expand=TRUE)
  windRaster<- resample(windRaster_new, windRaster, method= "bilinear")
  
}
writeRaster(windRaster, "windRaster", format="ascii", datatype = 'FLT4S', overwrite=TRUE)

#calculate windshelter

ctrl = wind.shelter.prep(5, (wind*pi)/180, (windTol*pi)/180 ,2*i_max)  #wind, windTOl*
focal.function("windRaster.asc",fun=wind.shelter,prob = 0.5, control=ctrl,
               radius=5,search.mode="circle")
f <- list.files(pattern='windshelter.asc$', full.names=TRUE)
windshelter <- raster(f)
windshelter <- resample(windshelter, asc, method= "bilinear")


#### calculate ruggedness at different scales #####


for (i in 1:i_max ) {
  
  
  #calculate slope and aspect ####
  
  slope_name <- paste("slope", i, sep="")
  aspect_name <- paste("aspect", i, sep="")
  rsaga.geoprocessor("ta_morphometry",23,env=myenv,list(             
    DEM ="outputRas.sgrd",
    SLOPE = slope_name,
    ASPECT = aspect_name,
    SIZE=i,
    TOL_SLOPE="1.00000",
    TOL_CURVE="0.000100",
    EXPONENT="0.00000",
    ZSCALE="1.000000",
    CONSTRAIN=FALSE))
  rsaga.sgrd.to.esri(slope_name, slope_name,
                     format = "ascii", georef = "corner", prec = 2)
  rsaga.sgrd.to.esri(aspect_name, aspect_name,
                     format = "ascii", georef = "corner", prec = 2)
  
  
  
  
  #create raster object of slope raster 
  f <- list.files(pattern=paste(slope_name, ".asc$", sep=""), full.names=TRUE)
  slope <- raster(f)
  f <- list.files(pattern=paste(aspect_name, ".asc$", sep=""), full.names=TRUE)
  aspect <- raster(f)
  
  
  
  
  #convert to radians
  slope_rad <- slope*pi/180
  aspect_rad <- aspect*pi/180
  
  #calculate xyz components
  xy_raster <- sin(slope_rad)
  z_raster <- cos(slope_rad)
  x_raster <- sin(aspect_rad) * xy_raster
  y_raster <- cos(aspect_rad) * xy_raster
  
  xsum_raster <- focal(x_raster, w=matrix(1,3,3), fun=sum)
  ysum_raster <- focal(y_raster, w=matrix(1,3,3), fun=sum)
  zsum_raster <- focal(z_raster, w=matrix(1,3,3), fun=sum)
  
  result_raster <- sqrt((xsum_raster)^2 + (ysum_raster)^2 + (zsum_raster)^2)
  
  ruggedness_raster <- (1- (result_raster/9))
  rugg_name <- paste("ruggedness", i, sep="")
  writeRaster(ruggedness_raster, rugg_name, format="ascii", overwrite=TRUE)
  
}



#### Correction of snow surface roughness with slope ####

f <- list.files(pattern=paste("ruggedness", i_max, ".asc$", sep=""), full.names=TRUE)
rugg<- raster(f)

if (i_max > 1) { 
  
  
  f <- list.files(pattern=paste("slope", i_max, ".asc$", sep=""), full.names=TRUE)
  slp_coef <- as.matrix(raster(f))
  slp_coef <- 1- ((slp_coef-30)/30)
  slp_coef[slp_coef < 0] <- 0
  slp_coef[slp_coef > 1] <- 1
  slp_coef <- 1 + (slp_coef * (i_max -1))
  slp_coef <- round(slp_coef, digits = 0)
  #write.ascii.grid(data = slp_coef, "coef", header = head, write.header = TRUE, digits = 2,
  #                hdr.digits = 3 , dec = ".", georef = "corner")
  
  
  for (i in 1:(i_max -1)) {
    f <- list.files(pattern=paste("ruggedness", i, ".asc$", sep=""), full.names=TRUE)
    rugg_i<- raster(f)
    rugg[which(slp_coef == i )] <- rugg_i[which(slp_coef == i)]
  }
}




#### Definition of membership functions #####


# define bell curve parameters for roughness
a <- 0.01
b <- 5
c <- -0.007


rugg1 <- 1/ (1+((rugg-c)/a)^(2*b))

rugg1[rugg > 0.01] <- 0



# define bell curve parameters for slope
a <- 11
b <- 4
c <- 43


slope1 <- 1/ (1+((slope-c)/a)^(2*b))
slope1[slope < 28] <- 0
slope1[slope > 60] <- 0


# define bell curve parametrs
a <- 2
b <- 5
c <- 2


windshelter <- 1/ (1+((windshelter-c)/a)^(2*b))


#### Fuzzy logic operator

minvar<- min(slope1, rugg1, windshelter)

PRA <- (1- minvar)*minvar + minvar*(slope1 + rugg1 + windshelter)/3


#forest mask
#forest <- raster(forest)
#PRA <- crop(PRA, forest)
#PRA <- PRA* (1-forest)

PRA.expand <- extend(PRA, asc.extent, value=NA)
writeRaster(PRA.expand, outPRA, format="ascii", overwrite=TRUE)

print("Calculations Complete...")

