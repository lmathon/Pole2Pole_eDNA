library(tidyverse)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(gissr)
library(move)


## open metadata uptodate
metadata_sampling <- read.csv("metadata/Metadata_eDNA_Pole2Pole.csv", header = T, sep = ";", stringsAsFactors = F, na.strings=c("","NA"))

## select samples with GPS data
metadata_dist  <- metadata_sampling %>%
  filter(!is.na(longitude_start)) %>%
  distinct(station, .keep_all=T)

## distance for whole world (EPSG:3857)

# Define projections
crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

crsmerc=CRS("+proj=merc +a=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs")

# Formate GPS points and project
pts = metadata_dist[,c("longitude_start", "latitude_start")]
pts$latitude_start <- as.numeric(pts$latitude_start)
pts$longitude_start <- as.numeric(pts$longitude_start)
pts_sp <- SpatialPoints(pts,proj4string = crswgs84)
pts_sp2 <- spTransform(pts_sp, crsmerc)

# calculate and formate distance
dist <- gDistance(pts_sp, pts_sp, byid = T)
dist <- dist*100
dist_km <- round(dist, 3) 
rownames(dist) <- metadata_dist$station
colnames(dist) <- metadata_dist$station

save(dist_km, file="Rdata/geographic_distance_stations.rdata")
