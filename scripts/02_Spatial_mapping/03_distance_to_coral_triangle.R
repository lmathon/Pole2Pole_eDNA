library(tidyverse)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(gissr)
library(raster)


## open metadata uptodate
metadata_sampling <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v2.csv", header = T, sep = ",", stringsAsFactors = F, na.strings=c("","NA"))

## select samples with GPS data
metadata_dist  <- metadata_sampling %>%
  filter(!is.na(longitude_start))
metadata_na  <- metadata_sampling %>%
  filter(is.na(longitude_start))


# coordinates of the center of the Coral Triangle
center_CT <- data.frame(longitude=133.679826, latitude=-1.307436)


# calculate distance
metadata_dist$dist_to_CT <- pointDistance(metadata_dist[,c("longitude_start", "latitude_start")], center_CT, lonlat=TRUE)


# transform in KM with 2 decimals
metadata_dist$dist_to_CT <- metadata_dist$dist_to_CT/1000

metadata_dist$dist_to_CT <- round(metadata_dist$dist_to_CT, 2) 


metadata_na$dist_to_CT <- NA
metadata <- rbind(metadata_dist, metadata_na)

write.csv(metadata, "metadata/Metadata_eDNA_Pole2Pole_v3.csv", row.names = F)
