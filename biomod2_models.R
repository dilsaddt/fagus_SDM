library(raster)
library(rgdal)
library(biomod2)

setwd("~/Desktop/Fagus")
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
  models = c('GLM','GAM','RF','SRE','MAXENT.Phillips'),
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
myBiomodModelEval["TSS","Testing.data","GLM",,]
myBiomodModelEval["TSS","Testing.data","GAM",,]
myBiomodModelEval["TSS","Testing.data","RF",,]
myBiomodModelEval["TSS","Testing.data","SRE",,]
myBiomodModelEval["TSS","Testing.data","MAXENT.Phillips",,]


# let's print the ROC scores of all selected models
myBiomodModelEval["ROC","Testing.data",,,]
# print variable importances
get_variables_importance(myBiomodModelOut)

# Ensemble Modeling
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )
# print summary
myBiomodEM
# get evaluation scores
EM <- get_evaluations(myBiomodEM)
EM
plot(myBiomodEM)


# Model Projection
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
plot(myBiomodProj, str.grep = 'RUN1_GLM')
plot(myBiomodProj, str.grep = 'RUN1_GAM')
plot(myBiomodProj, str.grep = 'RUN1_RF')
plot(myBiomodProj, str.grep = 'RUN1_SRE')
plot(myBiomodProj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj
writeRaster(myCurrentProj$Species_PA1_RUN1_RF, filename = "current_RF_out", format="GTiff")
writeRaster(myCurrentProj$Species_PA1_RUN1_GLM, filename = "current_GLM_out", format="GTiff")
writeRaster(myCurrentProj$Species_PA1_RUN1_GAM, filename = "current_GAM_out", format="GTiff")
writeRaster(myCurrentProj$Species_PA1_RUN1_SRE, filename = "current_SRE_out", format="GTiff")
writeRaster(myCurrentProj$Species_PA1_RUN1_MAXENT.Phillips, filename = "current_MAXENT.Phillips_out", format="GTiff")



# ##CLIMATE CHANGE PROJECTIONS --------------------------------------------

# # MIROC # ---------------------------------------------------------------

# #Last Glacial Maximum MIROC ---------------------------------------------

a_bio1 <- raster("~/Desktop/fagus_project/climate/lgm_miroc_2.5_tif/bio1.tif")
a_bio2 <- raster("~/Desktop/fagus_project/climate/lgm_miroc_2.5_tif/bio2.tif")
a_bio3 <- raster("~/Desktop/fagus_project/climate/lgm_miroc_2.5_tif/bio3.tif")
a_bio4 <- raster("~/Desktop/fagus_project/climate/lgm_miroc_2.5_tif/bio4.tif")
a_bio8 <- raster("~/Desktop/fagus_project/climate/lgm_miroc_2.5_tif/bio8.tif")
a_bio9 <- raster("~/Desktop/fagus_project/climate/lgm_miroc_2.5_tif/bio9.tif")
a_bio12 <- raster("~/Desktop/fagus_project/climate/lgm_miroc_2.5_tif/bio12.tif")
a_bio15 <- raster("~/Desktop/fagus_project/climate/lgm_miroc_2.5_tif/bio15.tif")
a_bio19 <- raster("~/Desktop/fagus_project/climate/lgm_miroc_2.5_tif/bio19.tif")

a_bios = stack(a_bio1, a_bio2, a_bio3, a_bio4, a_bio8, a_bio9, 
               a_bio12, a_bio15, a_bio19)


# projection under lgm_miroc conditions
LGM_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = a_bios,
  proj.name = 'LGM_MIROC',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
LGM_MIROC_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_LGM_MIROC")
# make some plots sub-selected by str.grep argument
plot(LGM_MIROC_Proj, str.grep = 'RUN1_GLM')
plot(LGM_MIROC_Proj, str.grep = 'RUN1_GAM')
plot(LGM_MIROC_Proj, str.grep = 'RUN1_RF')
plot(LGM_MIROC_Proj, str.grep = 'RUN1_SRE')
plot(LGM_MIROC_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
myLGMMProj <- get_predictions(LGM_MIROC_Proj)
myLGMMProj

writeRaster(myLGMMProj$Species_PA1_RUN1_RF, filename = "LGMM_RF", format="GTiff")
writeRaster(myLGMMProj$Species_PA1_RUN1_GLM, filename = "LGMM_GLM", format="GTiff")
writeRaster(myLGMMProj$Species_PA1_RUN1_GAM, filename = "LGMM_GAM", format="GTiff")
writeRaster(myLGMMProj$Species_PA1_RUN1_SRE, filename = "LGMM_SRE", format="GTiff")
writeRaster(myLGMMProj$Species_PA1_RUN1_MAXENT.Phillips, filename = "LGMM_MAXENT.Phillips", format="GTiff")



# #Mid-Holocene MIROC -----------------------------------------------------

b_bio1 <- raster("~/Desktop/fagus_project/climate/mh_miroc_2.5_tif/bio1.tif")
b_bio2 <- raster("~/Desktop/fagus_project/climate/mh_miroc_2.5_tif/bio2.tif")
b_bio3 <- raster("~/Desktop/fagus_project/climate/mh_miroc_2.5_tif/bio3.tif")
b_bio4 <- raster("~/Desktop/fagus_project/climate/mh_miroc_2.5_tif/bio4.tif")
b_bio8 <- raster("~/Desktop/fagus_project/climate/mh_miroc_2.5_tif/bio8.tif")
b_bio9 <- raster("~/Desktop/fagus_project/climate/mh_miroc_2.5_tif/bio9.tif")
b_bio12 <- raster("~/Desktop/fagus_project/climate/mh_miroc_2.5_tif/bio12.tif")
b_bio15 <- raster("~/Desktop/fagus_project/climate/mh_miroc_2.5_tif/bio15.tif")
b_bio19 <- raster("~/Desktop/fagus_project/climate/mh_miroc_2.5_tif/bio19.tif")

b_bios = stack(b_bio1, b_bio2, b_bio3, b_bio4, b_bio8, b_bio9, 
               b_bio12, b_bio15, b_bio19)


# projection under mh_miroc conditions
MH_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = b_bios,
  proj.name = 'MH MIROC',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
MH_MIROC_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_MH MIROC")
# make some plots sub-selected by str.grep argument
plot(MH_MIROC_Proj, str.grep = 'RUN1_GLM')
plot(MH_MIROC_Proj, str.grep = 'RUN1_GAM')
plot(MH_MIROC_Proj, str.grep = 'RUN1_RF')
plot(MH_MIROC_Proj, str.grep = 'RUN1_SRE')
plot(MH_MIROC_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
myMHMProj <- get_predictions(MH_MIROC_Proj)
myMHMProj

writeRaster(myMHMProj$Species_PA1_RUN1_RF, filename = "MHM_RF", format="GTiff")
writeRaster(myMHMProj$Species_PA1_RUN1_GLM, filename = "MHM_GLM", format="GTiff")
writeRaster(myMHMProj$Species_PA1_RUN1_GAM, filename = "MHM_GAM", format="GTiff")
writeRaster(myMHMProj$Species_PA1_RUN1_SRE, filename = "MHM_SRE", format="GTiff")
writeRaster(myMHMProj$Species_PA1_RUN1_MAXENT.Phillips, filename = "MHM_MAXENT.Phillips", format="GTiff")



# #2050 MIROC 8.5 ---------------------------------------------------------

c_bio1 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp85_tif/bio1.tif")
c_bio2 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp85_tif/bio2.tif")
c_bio3 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp85_tif/bio3.tif")
c_bio4 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp85_tif/bio4.tif")
c_bio8 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp85_tif/bio8.tif")
c_bio9 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp85_tif/bio9.tif")
c_bio12 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp85_tif/bio12.tif")
c_bio15 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp85_tif/bio15.tif")
c_bio19 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp85_tif/bio19.tif")

c_bios = stack(c_bio1, c_bio2, c_bio3, c_bio4, c_bio8, c_bio9, 
               c_bio12, c_bio15, c_bio19)

# projection under 2050_miroc conditions
future1_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = c_bios,
  proj.name = '2050 MIROC RCP8.5',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future1_MIROC_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2050 MIROC RCP8.5")
# make some plots sub-selected by str.grep argument
plot(future1_MIROC_Proj, str.grep = 'RUN1_GLM')
plot(future1_MIROC_Proj, str.grep = 'RUN1_GAM')
plot(future1_MIROC_Proj, str.grep = 'RUN1_RF')
plot(future1_MIROC_Proj, str.grep = 'RUN1_SRE')
plot(future1_MIROC_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2050MProj <- get_predictions(future1_MIROC_Proj)
my2050MProj

writeRaster(my2050MProj$Species_PA1_RUN1_RF, filename = "2050M_85_RF", format="GTiff")
writeRaster(my2050MProj$Species_PA1_RUN1_GLM, filename = "2050M_85_GLM", format="GTiff")
writeRaster(my2050MProj$Species_PA1_RUN1_GAM, filename = "2050M_85_GAM", format="GTiff")
writeRaster(my2050MProj$Species_PA1_RUN1_SRE, filename = "2050M_85_SRE", format="GTiff")
writeRaster(my2050MProj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2050M_85_MAXENT.Phillips", format="GTiff")



# #2050 MIROC RCP2.6 ------------------------------------------------------

i_bio1 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp26_tif/bio1.tif")
i_bio2 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp26_tif/bio2.tif")
i_bio3 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp26_tif/bio3.tif")
i_bio4 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp26_tif/bio4.tif")
i_bio8 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp26_tif/bio8.tif")
i_bio9 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp26_tif/bio9.tif")
i_bio12 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp26_tif/bio12.tif")
i_bio15 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp26_tif/bio15.tif")
i_bio19 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp26_tif/bio19.tif")

i_bios = stack(i_bio1, i_bio2, i_bio3, i_bio4, i_bio8, i_bio9, 
               i_bio12, i_bio15, i_bio19)

# projection under 2050_miroc conditions
future3_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = i_bios,
  proj.name = '2050 MIROC RCP2.6',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future3_MIROC_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2050 MIROC RCP2.6")
# make some plots sub-selected by str.grep argument
plot(future3_MIROC_Proj, str.grep = 'RUN1_GLM')
plot(future3_MIROC_Proj, str.grep = 'RUN1_GAM')
plot(future3_MIROC_Proj, str.grep = 'RUN1_RF')
plot(future3_MIROC_Proj, str.grep = 'RUN1_SRE')
plot(future3_MIROC_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2050M26Proj <- get_predictions(future3_MIROC_Proj)
my2050M26Proj

writeRaster(my2050M26Proj$Species_PA1_RUN1_RF, filename = "2050M_26_RF", format="GTiff")
writeRaster(my2050M26Proj$Species_PA1_RUN1_GLM, filename = "2050M_26_GLM", format="GTiff")
writeRaster(my2050M26Proj$Species_PA1_RUN1_GAM, filename = "2050M_26_GAM", format="GTiff")
writeRaster(my2050M26Proj$Species_PA1_RUN1_SRE, filename = "2050M_26_SRE", format="GTiff")
writeRaster(my2050M26Proj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2050M_26_MAXENT.Phillips", format="GTiff")



# #2050 MIROC RCP4.5 ------------------------------------------------------

m_bio1 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp45_tif/bio1.tif")
m_bio2 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp45_tif/bio2.tif")
m_bio3 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp45_tif/bio3.tif")
m_bio4 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp45_tif/bio4.tif")
m_bio8 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp45_tif/bio8.tif")
m_bio9 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp45_tif/bio9.tif")
m_bio12 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp45_tif/bio12.tif")
m_bio15 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp45_tif/bio15.tif")
m_bio19 <- raster("~/Desktop/fagus_project/climate/2050_miroc_2.5_rcp45_tif/bio19.tif")

m_bios = stack(m_bio1, m_bio2, m_bio3, m_bio4, m_bio8, m_bio9, 
               m_bio12, m_bio15, m_bio19)

# projection under 2050_miroc conditions
future5_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = m_bios,
  proj.name = '2050 MIROC RCP4.5',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future5_MIROC_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2050 MIROC RCP4.5")
# make some plots sub-selected by str.grep argument
plot(future5_MIROC_Proj, str.grep = 'RUN1_GLM')
plot(future5_MIROC_Proj, str.grep = 'RUN1_GAM')
plot(future5_MIROC_Proj, str.grep = 'RUN1_RF')
plot(future5_MIROC_Proj, str.grep = 'RUN1_SRE')
plot(future5_MIROC_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2050M45Proj <- get_predictions(future5_MIROC_Proj)
my2050M45Proj

writeRaster(my2050M45Proj$Species_PA1_RUN1_RF, filename = "2050M_45_RF", format="GTiff")
writeRaster(my2050M45Proj$Species_PA1_RUN1_GLM, filename = "2050M_45_GLM", format="GTiff")
writeRaster(my2050M45Proj$Species_PA1_RUN1_GAM, filename = "2050M_45_GAM", format="GTiff")
writeRaster(my2050M45Proj$Species_PA1_RUN1_SRE, filename = "2050M_45_SRE", format="GTiff")
writeRaster(my2050M45Proj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2050M_45_MAXENT.Phillips", format="GTiff")



# #2070 MIROC 8.5 ---------------------------------------------------------

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
  proj.name = '2070 MIROC RCP8.5',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future2_MIROC_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2070 MIROC RCP8.5")
# make some plots sub-selected by str.grep argument
plot(future2_MIROC_Proj, str.grep = 'RUN1_GLM')
plot(future2_MIROC_Proj, str.grep = 'RUN1_GAM')
plot(future2_MIROC_Proj, str.grep = 'RUN1_RF')
plot(future2_MIROC_Proj, str.grep = 'RUN1_SRE')
plot(future2_MIROC_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2070MProj <- get_predictions(future2_MIROC_Proj)
my2070MProj
writeRaster(my2070MProj$Species_PA1_RUN1_RF, filename = "2070M_85_RF", format="GTiff")
writeRaster(my2070MProj$Species_PA1_RUN1_GLM, filename = "2070M_85_GLM", format="GTiff")
writeRaster(my2070MProj$Species_PA1_RUN1_GAM, filename = "2070M_85_GAM", format="GTiff")
writeRaster(my2070MProj$Species_PA1_RUN1_SRE, filename = "2070M_85_SRE", format="GTiff")
writeRaster(my2070MProj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2070M_85_MAXENT.Phillips", format="GTiff")


# #2070 MIROC RCP2.6 ------------------------------------------------------

j_bio1 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio1.tif")
j_bio2 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio2.tif")
j_bio3 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio3.tif")
j_bio4 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio4.tif")
j_bio8 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio8.tif")
j_bio9 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio9.tif")
j_bio12 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio12.tif")
j_bio15 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio15.tif")
j_bio19 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp26_tif/bio19.tif")

j_bios = stack(j_bio1, j_bio2, j_bio3, j_bio4, j_bio8, j_bio9, 
               j_bio12, j_bio15, j_bio19)

# projection under 2070_miroc conditions
future4_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = j_bios,
  proj.name = '2070 MIROC RCP2.6',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future4_MIROC_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2070 MIROC RCP2.6")
# make some plots sub-selected by str.grep argument
plot(future4_MIROC_Proj, str.grep = 'RUN1_GLM')
plot(future4_MIROC_Proj, str.grep = 'RUN1_GAM')
plot(future4_MIROC_Proj, str.grep = 'RUN1_RF')
plot(future4_MIROC_Proj, str.grep = 'RUN1_SRE')
plot(future4_MIROC_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2070M26Proj <- get_predictions(future4_MIROC_Proj)
my2070M26Proj

writeRaster(my2070M26Proj$Species_PA1_RUN1_RF, filename = "2070M_26_RF", format="GTiff")
writeRaster(my2070M26Proj$Species_PA1_RUN1_GLM, filename = "2070M_26_GLM", format="GTiff")
writeRaster(my2070M26Proj$Species_PA1_RUN1_GAM, filename = "2070M_26_GAM", format="GTiff")
writeRaster(my2070M26Proj$Species_PA1_RUN1_SRE, filename = "2070M_26_SRE", format="GTiff")
writeRaster(my2070M26Proj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2070M_26_MAXENT.Phillips", format="GTiff")


# #2070 MIROC RCP4.5 ------------------------------------------------------

n_bio1 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio1.tif")
n_bio2 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio2.tif")
n_bio3 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio3.tif")
n_bio4 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio4.tif")
n_bio8 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio8.tif")
n_bio9 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio9.tif")
n_bio12 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio12.tif")
n_bio15 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio15.tif")
n_bio19 <- raster("~/Desktop/fagus_project/climate/2070_miroc_2.5_rcp45_tif/bio19.tif")

n_bios = stack(n_bio1, n_bio2, n_bio3, n_bio4, n_bio8, n_bio9, 
               n_bio12, n_bio15, n_bio19)

# projection under 2070_miroc conditions
future6_MIROC_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = n_bios,
  proj.name = '2070 MIROC RCP4.5',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future6_MIROC_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2070 MIROC RCP4.5")
# make some plots sub-selected by str.grep argument
plot(future6_MIROC_Proj, str.grep = 'RUN1_GLM')
plot(future6_MIROC_Proj, str.grep = 'RUN1_GAM')
plot(future6_MIROC_Proj, str.grep = 'RUN1_RF')
plot(future6_MIROC_Proj, str.grep = 'RUN1_SRE')
plot(future6_MIROC_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2070M45Proj <- get_predictions(future6_MIROC_Proj)
my2070M45Proj

writeRaster(my2070M45Proj$Species_PA1_RUN1_RF, filename = "2070M_45_RF", format="GTiff")
writeRaster(my2070M45Proj$Species_PA1_RUN1_GLM, filename = "2070M_45_GLM", format="GTiff")
writeRaster(my2070M45Proj$Species_PA1_RUN1_GAM, filename = "2070M_45_GAM", format="GTiff")
writeRaster(my2070M45Proj$Species_PA1_RUN1_SRE, filename = "2070M_45_SRE", format="GTiff")
writeRaster(my2070M45Proj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2070M_45_MAXENT.Phillips", format="GTiff")



# # CCSM4 # ---------------------------------------------------------------


# #Last Glacial Maximum CCSM4 ---------------------------------------------

e_bio1 <- raster("~/Desktop/fagus_project/climate/lgm_ccsm4_2.5_tif/bio1.tif")
e_bio2 <- raster("~/Desktop/fagus_project/climate/lgm_ccsm4_2.5_tif/bio2.tif")
e_bio3 <- raster("~/Desktop/fagus_project/climate/lgm_ccsm4_2.5_tif/bio3.tif")
e_bio4 <- raster("~/Desktop/fagus_project/climate/lgm_ccsm4_2.5_tif/bio4.tif")
e_bio8 <- raster("~/Desktop/fagus_project/climate/lgm_ccsm4_2.5_tif/bio8.tif")
e_bio9 <- raster("~/Desktop/fagus_project/climate/lgm_ccsm4_2.5_tif/bio9.tif")
e_bio12 <- raster("~/Desktop/fagus_project/climate/lgm_ccsm4_2.5_tif/bio12.tif")
e_bio15 <- raster("~/Desktop/fagus_project/climate/lgm_ccsm4_2.5_tif/bio15.tif")
e_bio19 <- raster("~/Desktop/fagus_project/climate/lgm_ccsm4_2.5_tif/bio19.tif")

e_bios = stack(e_bio1, e_bio2, e_bio3, e_bio4, e_bio8, e_bio9, 
               e_bio12, e_bio15, e_bio19)


# projection under lgm_miroc conditions
LGM_CCSM4_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = e_bios,
  proj.name = 'LGM_CCSM4',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
LGM_CCSM4_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_LGM_CCSM4")
# make some plots sub-selected by str.grep argument
plot(LGM_CCSM4_Proj, str.grep = 'RUN1_GLM')
plot(LGM_CCSM4_Proj, str.grep = 'RUN1_GAM')
plot(LGM_CCSM4_Proj, str.grep = 'RUN1_RF')
plot(LGM_CCSM4_Proj, str.grep = 'RUN1_SRE')
plot(LGM_CCSM4_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
myLGMCProj <- get_predictions(LGM_CCSM4_Proj)
myLGMCProj

writeRaster(myLGMCProj$Species_PA1_RUN1_RF, filename = "LGMC_RF", format="GTiff")
writeRaster(myLGMCProj$Species_PA1_RUN1_GLM, filename = "LGMC_GLM", format="GTiff")
writeRaster(myLGMCProj$Species_PA1_RUN1_GAM, filename = "LGMC_GAM", format="GTiff")
writeRaster(myLGMCProj$Species_PA1_RUN1_SRE, filename = "LGMC_SRE", format="GTiff")
writeRaster(myLGMCProj$Species_PA1_RUN1_MAXENT.Phillips, filename = "LGMC_MAXENT.Phillips", format="GTiff")


# #Mid-Holocene CCSM4 -----------------------------------------------------

f_bio1 <- raster("~/Desktop/fagus_project/climate/mh_ccsm4_2.5_tif/bio1.tif")
f_bio2 <- raster("~/Desktop/fagus_project/climate/mh_ccsm4_2.5_tif/bio2.tif")
f_bio3 <- raster("~/Desktop/fagus_project/climate/mh_ccsm4_2.5_tif/bio3.tif")
f_bio4 <- raster("~/Desktop/fagus_project/climate/mh_ccsm4_2.5_tif/bio4.tif")
f_bio8 <- raster("~/Desktop/fagus_project/climate/mh_ccsm4_2.5_tif/bio8.tif")
f_bio9 <- raster("~/Desktop/fagus_project/climate/mh_ccsm4_2.5_tif/bio9.tif")
f_bio12 <- raster("~/Desktop/fagus_project/climate/mh_ccsm4_2.5_tif/bio12.tif")
f_bio15 <- raster("~/Desktop/fagus_project/climate/mh_ccsm4_2.5_tif/bio15.tif")
f_bio19 <- raster("~/Desktop/fagus_project/climate/mh_ccsm4_2.5_tif/bio19.tif")

f_bios = stack(f_bio1, f_bio2, f_bio3, f_bio4, f_bio8, f_bio9, 
               f_bio12, f_bio15, f_bio19)


# projection under mh_miroc conditions
MH_CCSM4_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = f_bios,
  proj.name = 'MH CCSM4',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
MH_CCSM4_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_MH MIROC")
# make some plots sub-selected by str.grep argument
plot(MH_CCSM4_Proj, str.grep = 'RUN1_GLM')
plot(MH_CCSM4_Proj, str.grep = 'RUN1_GAM')
plot(MH_CCSM4_Proj, str.grep = 'RUN1_RF')
plot(MH_CCSM4_Proj, str.grep = 'RUN1_SRE')
plot(MH_CCSM4_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
myMHCProj <- get_predictions(MH_CCSM4_Proj)
myMHCProj

writeRaster(myMHCProj$Species_PA1_RUN1_RF, filename = "MHC_RF", format="GTiff")
writeRaster(myMHCProj$Species_PA1_RUN1_GLM, filename = "MHC_GLM", format="GTiff")
writeRaster(myMHCProj$Species_PA1_RUN1_GAM, filename = "MHC_GAM", format="GTiff")
writeRaster(myMHCProj$Species_PA1_RUN1_SRE, filename = "MHC_SRE", format="GTiff")
writeRaster(myMHCProj$Species_PA1_RUN1_MAXENT.Phillips, filename = "MHC_MAXENT.Phillips", format="GTiff")


# #2050 CCSM4 RCP 8.5 -----------------------------------------------------

g_bio1 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp85_tif/bio1.tif")
g_bio2 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp85_tif/bio2.tif")
g_bio3 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp85_tif/bio3.tif")
g_bio4 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp85_tif/bio4.tif")
g_bio8 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp85_tif/bio8.tif")
g_bio9 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp85_tif/bio9.tif")
g_bio12 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp85_tif/bio12.tif")
g_bio15 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp85_tif/bio15.tif")
g_bio19 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp85_tif/bio19.tif")

g_bios = stack(g_bio1, g_bio2, g_bio3, g_bio4, g_bio8, g_bio9, 
               g_bio12, g_bio15, g_bio19)

# projection under 2050_ccsm4 conditions
future1_CCSM4_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = g_bios,
  proj.name = '2050 CCSM4 RCP8.5',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future1_CCSM4_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2050 CCSM4 RCP8.5")
# make some plots sub-selected by str.grep argument
plot(future1_CCSM4_Proj, str.grep = 'RUN1_GLM')
plot(future1_CCSM4_Proj, str.grep = 'RUN1_GAM')
plot(future1_CCSM4_Proj, str.grep = 'RUN1_RF')
plot(future1_CCSM4_Proj, str.grep = 'RUN1_SRE')
plot(future1_CCSM4_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2050CProj <- get_predictions(future1_CCSM4_Proj)
my2050CProj

writeRaster(my2050CProj$Species_PA1_RUN1_RF, filename = "2050C_85_RF", format="GTiff")
writeRaster(my2050CProj$Species_PA1_RUN1_GLM, filename = "2050C_85_GLM", format="GTiff")
writeRaster(my2050CProj$Species_PA1_RUN1_GAM, filename = "2050C_85_GAM", format="GTiff")
writeRaster(my2050CProj$Species_PA1_RUN1_SRE, filename = "2050C_85_SRE", format="GTiff")
writeRaster(my2050CProj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2050C_85_MAXENT.Phillips", format="GTiff")


# #2050 CCSM4 RCP2.6 ------------------------------------------------------

k_bio1 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp26_tif/bio1.tif")
k_bio2 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp26_tif/bio2.tif")
k_bio3 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp26_tif/bio3.tif")
k_bio4 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp26_tif/bio4.tif")
k_bio8 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp26_tif/bio8.tif")
k_bio9 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp26_tif/bio9.tif")
k_bio12 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp26_tif/bio12.tif")
k_bio15 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp26_tif/bio15.tif")
k_bio19 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp26_tif/bio19.tif")

k_bios = stack(k_bio1, k_bio2, k_bio3, k_bio4, k_bio8, k_bio9, 
               k_bio12, k_bio15, k_bio19)

# projection under 2050_miroc conditions
future3_CCSM4_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = k_bios,
  proj.name = '2050 CCSM4 RCP2.6',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future3_CCSM4_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2050 CCSM4 RCP2.6")
# make some plots sub-selected by str.grep argument
plot(future3_CCSM4_Proj, str.grep = 'RUN1_GLM')
plot(future3_CCSM4_Proj, str.grep = 'RUN1_GAM')
plot(future3_CCSM4_Proj, str.grep = 'RUN1_RF')
plot(future3_CCSM4_Proj, str.grep = 'RUN1_SRE')
plot(future3_CCSM4_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2050C26Proj <- get_predictions(future3_CCSM4_Proj)
my2050C26Proj

writeRaster(my2050C26Proj$Species_PA1_RUN1_RF, filename = "2050C_26_RF", format="GTiff")
writeRaster(my2050C26Proj$Species_PA1_RUN1_GLM, filename = "2050C_26_GLM", format="GTiff")
writeRaster(my2050C26Proj$Species_PA1_RUN1_GAM, filename = "2050C_26_GAM", format="GTiff")
writeRaster(my2050C26Proj$Species_PA1_RUN1_SRE, filename = "2050C_26_SRE", format="GTiff")
writeRaster(my2050C26Proj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2050C_26_MAXENT.Phillips", format="GTiff")


# #2050 CCSM4 RCP4.5 ------------------------------------------------------

o_bio1 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp45_tif/bio1.tif")
o_bio2 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp45_tif/bio2.tif")
o_bio3 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp45_tif/bio3.tif")
o_bio4 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp45_tif/bio4.tif")
o_bio8 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp45_tif/bio8.tif")
o_bio9 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp45_tif/bio9.tif")
o_bio12 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp45_tif/bio12.tif")
o_bio15 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp45_tif/bio15.tif")
o_bio19 <- raster("~/Desktop/fagus_project/climate/2050_ccsm4_2.5_rcp45_tif/bio19.tif")

o_bios = stack(o_bio1, o_bio2, o_bio3, o_bio4, o_bio8, o_bio9, 
               o_bio12, o_bio15, o_bio19)

# projection under 2050_miroc conditions
future5_CCSM4_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = o_bios,
  proj.name = '2050 MIROC RCP4.5',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future5_CCSM4_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2050 CCSM4 RCP4.5")
# make some plots sub-selected by str.grep argument
plot(future5_CCSM4_Proj, str.grep = 'RUN1_GLM')
plot(future5_CCSM4_Proj, str.grep = 'RUN1_GAM')
plot(future5_CCSM4_Proj, str.grep = 'RUN1_RF')
plot(future5_CCSM4_Proj, str.grep = 'RUN1_SRE')
plot(future5_CCSM4_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2050C45Proj <- get_predictions(future5_CCSM4_Proj)
my2050C45Proj

writeRaster(my2050C45Proj$Species_PA1_RUN1_RF, filename = "2050C_45_RF", format="GTiff")
writeRaster(my2050C45Proj$Species_PA1_RUN1_GLM, filename = "2050C_45_GLM", format="GTiff")
writeRaster(my2050C45Proj$Species_PA1_RUN1_GAM, filename = "2050C_45_GAM", format="GTiff")
writeRaster(my2050C45Proj$Species_PA1_RUN1_SRE, filename = "2050C_45_SRE", format="GTiff")
writeRaster(my2050C45Proj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2050C_45_MAXENT.Phillips", format="GTiff")

# #2070 CCSM4 RCP8.5 ------------------------------------------------------

h_bio1 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp85_tif/bio1.tif")
h_bio2 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp85_tif/bio2.tif")
h_bio3 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp85_tif/bio3.tif")
h_bio4 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp85_tif/bio4.tif")
h_bio8 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp85_tif/bio8.tif")
h_bio9 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp85_tif/bio9.tif")
h_bio12 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp85_tif/bio12.tif")
h_bio15 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp85_tif/bio15.tif")
h_bio19 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp85_tif/bio19.tif")

h_bios = stack(h_bio1, h_bio2, h_bio3, h_bio4, h_bio8, h_bio9, 
               h_bio12, h_bio15, h_bio19)

# projection under 2070_miroc conditions
future2_CCSM4_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = h_bios,
  proj.name = '2070 CCSM4 RCP8.5',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future2_CCSM4_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2070 CCSM4 RCP8.5")
# make some plots sub-selected by str.grep argument
plot(future2_CCSM4_Proj, str.grep = 'RUN1_GLM')
plot(future2_CCSM4_Proj, str.grep = 'RUN1_GAM')
plot(future2_CCSM4_Proj, str.grep = 'RUN1_RF')
plot(future2_CCSM4_Proj, str.grep = 'RUN1_SRE')
plot(future2_CCSM4_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2070CProj <- get_predictions(future2_CCSM4_Proj)
my2070CProj

writeRaster(my2070CProj$Species_PA1_RUN1_RF, filename = "2070C_85_RF", format="GTiff")
writeRaster(my2070CProj$Species_PA1_RUN1_GLM, filename = "2070C_85_GLM", format="GTiff")
writeRaster(my2070CProj$Species_PA1_RUN1_GAM, filename = "2070C_85_GAM", format="GTiff")
writeRaster(my2070CProj$Species_PA1_RUN1_SRE, filename = "2070C_85_SRE", format="GTiff")
writeRaster(my2070CProj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2070C_85_MAXENT.Phillips", format="GTiff")


# #2070 CCSM4 RCP2.6 ------------------------------------------------------

l_bio1 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp26_tif/bio1.tif")
l_bio2 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp26_tif/bio2.tif")
l_bio3 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp26_tif/bio3.tif")
l_bio4 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp26_tif/bio4.tif")
l_bio8 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp26_tif/bio8.tif")
l_bio9 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp26_tif/bio9.tif")
l_bio12 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp26_tif/bio12.tif")
l_bio15 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp26_tif/bio15.tif")
l_bio19 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp26_tif/bio19.tif")

l_bios = stack(l_bio1, l_bio2, l_bio3, l_bio4, l_bio8, l_bio9, 
               l_bio12, l_bio15, l_bio19)

# projection under 2070_miroc conditions
future4_CCSM4_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = l_bios,
  proj.name = '2070 CCSM4 RCP2.6',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future4_CCSM4_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2070 CCSM4 RCP2.6")
# make some plots sub-selected by str.grep argument
plot(future4_CCSM4_Proj, str.grep = 'RUN1_GLM')
plot(future4_CCSM4_Proj, str.grep = 'RUN1_GAM')
plot(future4_CCSM4_Proj, str.grep = 'RUN1_RF')
plot(future4_CCSM4_Proj, str.grep = 'RUN1_SRE')
plot(future4_CCSM4_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2070C26Proj <- get_predictions(future4_CCSM4_Proj)
my2070C26Proj

writeRaster(my2070C26Proj$Species_PA1_RUN1_RF, filename = "2070C_26_RF", format="GTiff")
writeRaster(my2070C26Proj$Species_PA1_RUN1_GLM, filename = "2070C_26_GLM", format="GTiff")
writeRaster(my2070C26Proj$Species_PA1_RUN1_GAM, filename = "2070C_26_GAM", format="GTiff")
writeRaster(my2070C26Proj$Species_PA1_RUN1_SRE, filename = "2070C_26_SRE", format="GTiff")
writeRaster(my2070C26Proj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2070C_26_MAXENT.Phillips", format="GTiff")


# #2070 CCSM4 RCP4.5 ------------------------------------------------------

p_bio1 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp45_tif/bio1.tif")
p_bio2 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp45_tif/bio2.tif")
p_bio3 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp45_tif/bio3.tif")
p_bio4 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp45_tif/bio4.tif")
p_bio8 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp45_tif/bio8.tif")
p_bio9 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp45_tif/bio9.tif")
p_bio12 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp45_tif/bio12.tif")
p_bio15 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp45_tif/bio15.tif")
p_bio19 <- raster("~/Desktop/fagus_project/climate/2070_ccsm4_2.5_rcp45_tif/bio19.tif")

p_bios = stack(p_bio1, p_bio2, p_bio3, p_bio4, p_bio8, p_bio9, 
               p_bio12, p_bio15, p_bio19)

# projection under 2070_miroc conditions
future6_CCSM4_Proj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = p_bios,
  proj.name = '2070 CCSM4 RCP4.5',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of crated object
future6_CCSM4_Proj

# files created on hard drive
list.files("Fagus.orientalis/proj_2070 CCSM4 RCP4.5")
# make some plots sub-selected by str.grep argument
plot(future6_CCSM4_Proj, str.grep = 'RUN1_GLM')
plot(future6_CCSM4_Proj, str.grep = 'RUN1_GAM')
plot(future6_CCSM4_Proj, str.grep = 'RUN1_RF')
plot(future6_CCSM4_Proj, str.grep = 'RUN1_SRE')
plot(future6_CCSM4_Proj, str.grep = 'RUN1_MAXENT.Phillips')
# if you want to make custom plots, you can also get the projected map
my2070C45Proj <- get_predictions(future6_CCSM4_Proj)
my2070C45Proj

writeRaster(my2070C45Proj$Species_PA1_RUN1_RF, filename = "2070C_45_RF", format="GTiff")
writeRaster(my2070C45Proj$Species_PA1_RUN1_GLM, filename = "2070C_45_GLM", format="GTiff")
writeRaster(my2070C45Proj$Species_PA1_RUN1_GAM, filename = "2070C_45_GAM", format="GTiff")
writeRaster(my2070C45Proj$Species_PA1_RUN1_SRE, filename = "2070C_45_SRE", format="GTiff")
writeRaster(my2070C45Proj$Species_PA1_RUN1_MAXENT.Phillips, filename = "2070C_45_MAXENT.Phillips", format="GTiff")


