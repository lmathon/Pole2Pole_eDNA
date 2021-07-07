library(tidyverse)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(raster)


metadata_sampling <- read.csv("metadata/Metadata_eDNA_Pole2Pole.csv", header = T, sep = ";", stringsAsFactors = F, na.strings=c("","NA"))

## select samples with GPS data
metadata_dist  <- metadata_sampling %>%
  filter(!is.na(longitude_start))

metadata_na  <- metadata_sampling %>%
  filter(is.na(longitude_start))

b1 <- raster("c:/Users/mathon/Downloads/gebco_2020_geotiff/gebco_2020_n0.0_s-90.0_w-180.0_e-90.0.tif")
b2 <- raster("c:/Users/mathon/Downloads/gebco_2020_geotiff/gebco_2020_n0.0_s-90.0_w-90.0_e0.0.tif")
b3 <- raster("c:/Users/mathon/Downloads/gebco_2020_geotiff/gebco_2020_n0.0_s-90.0_w0.0_e90.0.tif")
b4 <- raster("c:/Users/mathon/Downloads/gebco_2020_geotiff/gebco_2020_n0.0_s-90.0_w90.0_e180.0.tif")
b5 <- raster("c:/Users/mathon/Downloads/gebco_2020_geotiff/gebco_2020_n90.0_s0.0_w-180.0_e-90.0.tif")
b6 <- raster("c:/Users/mathon/Downloads/gebco_2020_geotiff/gebco_2020_n90.0_s0.0_w-90.0_e0.0.tif")
b7 <- raster("c:/Users/mathon/Downloads/gebco_2020_geotiff/gebco_2020_n90.0_s0.0_w0.0_e90.0.tif")
b8 <- raster("c:/Users/mathon/Downloads/gebco_2020_geotiff/gebco_2020_n90.0_s0.0_w90.0_e180.0.tif")


bathy <- list(b1,b2,b3,b4,b5,b6,b7,b8)
bathy$filename <- 'bathy.tif'
bathy$overwrite <- TRUE

bathy_tot <- do.call(merge, bathy)


crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
pts = metadata_dist[,c("longitude_start", "latitude_start")]
pts$latitude_start <- as.numeric(pts$latitude_start)
pts$longitude_start <- as.numeric(pts$longitude_start)
pts_sp <- SpatialPoints(pts,proj4string = crswgs84)
pts_sp <- data.frame(pts_sp)
names(pts_sp) <- c("long", "lat")

prof <- cbind(extract(bathy_tot, pts_sp, df=T), pts_sp)

