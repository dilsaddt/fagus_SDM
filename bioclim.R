# cropping climatet data again
library(raster)
library(sp)
setwd("~/Desktop/fagus_project")


# get data
#CURRENT
bio1 <- raster("~/Desktop/fagus_project/bio_2-5m_bil/bio1.bil")
bio2 <- raster("~/Desktop/fagus_project/bio_2-5m_bil/bio2.bil")
bio3 <- raster("~/Desktop/fagus_project/bio_2-5m_bil/bio3.bil")
bio4 <- raster("~/Desktop/fagus_project/bio_2-5m_bil/bio4.bil")
bio8 <- raster("~/Desktop/fagus_project/bio_2-5m_bil/bio8.bil")
bio9 <- raster("~/Desktop/fagus_project/bio_2-5m_bil/bio9.bil")
bio12 <- raster("~/Desktop/fagus_project/bio_2-5m_bil/bio12.bil")
bio15 <- raster("~/Desktop/fagus_project/bio_2-5m_bil/bio15.bil")
bio19 <- raster("~/Desktop/fagus_project/bio_2-5m_bil/bio19.bil")

bios <- stack(bio1, bio2, bio3, bio4, bio8, bio9, bio12, bio15, bio19)

# crop
region <- extent(c(18, 62, 33, 51))
bio <- crop(bios, region)
plot(bio$bio1) # plotting one by one
plot(bio) # plotting each bio in the same frame