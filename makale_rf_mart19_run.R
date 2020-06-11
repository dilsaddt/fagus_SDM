library(biomod2)
library(raster)
library(rgdal)

setwd("~/Desktop/fagus_project")
#load the data
fagus <- read.csv("presence25.csv", header = TRUE, sep = ",")
head(fagus)
plot(fagus)

#getting presence values from data
myRespName <- 'Species'
myResp <- as.numeric(fagus[,myRespName])
#getting presence coordinates fromd data
myRespXY <- fagus[,c("X", "Y")]

#environmental layers as explanatory variables

bio1 <- raster("~/Desktop/fagus_project/climate/current_2.5_tif/bio1.tif")
bio2 <- raster("~/Desktop/fagus_project/climate/current_2.5_tif/bio2.tif")
bio3 <- raster("~/Desktop/fagus_project/climate/current_2.5_tif/bio3.tif")
bio4 <- raster("~/Desktop/fagus_project/climate/current_2.5_tif/bio4.tif")
bio8 <- raster("~/Desktop/fagus_project/climate/current_2.5_tif/bio8.tif")
bio9 <- raster("~/Desktop/fagus_project/climate/current_2.5_tif/bio9.tif")
bio12 <- raster("~/Desktop/fagus_project/climate/current_2.5_tif/bio12.tif")
bio15 <- raster("~/Desktop/fagus_project/climate/current_2.5_tif/bio15.tif")
bio19 <- raster("~/Desktop/fagus_project/climate/current_2.5_tif/bio19.tif")

current_bios = stack(bio1, bio2, bio3, bio4, bio8, bio9, bio12, bio15, bio19)

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = current_bios,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 20796,
                                     PA.strategy = 'random',
                                     na.rm = TRUE)

plot(myBiomodData)

#modeling

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GAM', 'RF', 'SRE', 'MAXENT.Phillips'),
  models.options = myBiomodOption,
  NbRunEval=3,
  DataSplit=70,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste('Fagus orientalis',"FirstModeling",sep=""))

myBiomodModelOut

#get all models evaluation
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
myBiomodModelEval
# print the dimnames of this object
dimnames(myBiomodModelEval)
# let's print the TSS scores
myBiomodModelEval["TSS","Testing.data","RF",,]
# let's print the ROC scores of all selected models
myBiomodModelEval["ROC","Testing.data",,,]
# print variable importances
get_variables_importance(myBiomodModelOut)

#Model Projection
# First let’s project the individual models on our current conditions
# (the globe) to visualize them.

# projection over the globe under current conditions
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = current_bios,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')
# summary of crated oject
myBiomodProj
# files created on hard drive
list.files("Fagus.orientalis/proj_current/")
# make some plots sub-selected by str.grep argument
plot(myBiomodProj, str.grep = 'RUN1_RF')
# if you want to make custom plots, you can also get the projected map
myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj
names(myCurrentProj)
writeRaster(myCurrentProj$Species_PA1_RUN1_RF, filename = "current_RF", format="GTiff")


#2070 MIROC
#8.5
d_bio1 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp85_tif/bio1.tif")
d_bio2 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp85_tif/bio2.tif")
d_bio3 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp85_tif/bio3.tif")
d_bio4 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp85_tif/bio4.tif")
d_bio8 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp85_tif/bio8.tif")
d_bio9 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp85_tif/bio9.tif")
d_bio12 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp85_tif/bio12.tif")
d_bio15 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp85_tif/bio15.tif")
d_bio19 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp85_tif/bio19.tif")

d_bios = stack(d_bio1, d_bio2, d_bio3, d_bio4, d_bio8, d_bio9, 
               d_bio12, d_bio15, d_bio19)

# projection under 2070_miroc conditions
future2_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = d_bios,
  proj.name = '2070 MIROC',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future2_MIROC_Proj
# files created on hard drive
list.files("Fagus.orientalis/proj_2070 MIROC")
# make some plots sub-selected by str.grep argument

plot(future2_MIROC_Proj, str.grep = 'RUN1_RF')

# if you want to make custom plots, you can also get the projected map
my2070MProj <- get_predictions(future2_MIROC_Proj)
my2070Mproj
writeRaster(my2070MProj$Species_PA1_RUN1_RF, filename = "2070M_85_RF", format="GTiff")

#4.5
e_bio1 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio1.tif")
e_bio2 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio2.tif")
e_bio3 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio3.tif")
e_bio4 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio4.tif")
e_bio8 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio8.tif")
e_bio9 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio9.tif")
e_bio12 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio12.tif")
e_bio15 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio15.tif")
e_bio19 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio19.tif")

e_bios = stack(e_bio1, e_bio2, e_bio3, e_bio4, e_bio8, e_bio9, 
               e_bio12, e_bio15, e_bio19)

# projection under 2070_miroc conditions
future2_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = e_bios,
  proj.name = '2070 MIROC',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future2_MIROC_Proj
# files created on hard drive
list.files("Fagus.orientalis/proj_2070 MIROC")
# make some plots sub-selected by str.grep argument

plot(future2_MIROC_Proj, str.grep = 'RUN1_RF')

# if you want to make custom plots, you can also get the projected map
my2070MProj <- get_predictions(future2_MIROC_Proj)
my2070Mproj
writeRaster(my2070MProj$Species_PA1_RUN1_RF, filename = "2070M_45_RF", format="GTiff")

#2.6
f_bio1 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio1.tif")
f_bio2 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio2.tif")
f_bio3 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio3.tif")
f_bio4 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio4.tif")
f_bio8 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio8.tif")
f_bio9 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio9.tif")
f_bio12 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio12.tif")
f_bio15 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio15.tif")
f_bio19 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio19.tif")

f_bios = stack(f_bio1, f_bio2, f_bio3, f_bio4, f_bio8, f_bio9, 
               f_bio12, f_bio15, f_bio19)

# projection under 2070_miroc conditions
future2_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = f_bios,
  proj.name = '2070 MIROC',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future2_MIROC_Proj
# files created on hard drive
list.files("Fagus.orientalis/proj_2070 MIROC")
# make some plots sub-selected by str.grep argument

plot(future2_MIROC_Proj, str.grep = 'RUN1_RF')

# if you want to make custom plots, you can also get the projected map
my2070MProj <- get_predictions(future2_MIROC_Proj)
my2070Mproj
writeRaster(my2070MProj$Species_PA1_RUN1_RF, filename = "2070M_26_RF", format="GTiff")

