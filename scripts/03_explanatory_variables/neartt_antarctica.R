library(tidyverse)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(gissr)
library(raster)
'%ni%' <- Negate("%in%")

## open metadata uptodate
met_socio <- read.csv("metadata/eDNA_all_socio_variables.csv", header = T, sep = ";", stringsAsFactors = F, na.strings=c("","NA"))

antarctic <- met_socio %>%
  filter(TERRITORY1 %in% c("Antarctica", "Argentina"))%>%
  distinct(station, .keep_all=T) %>%
  dplyr::select(-neartt)

no_antarctic <- met_socio %>%
  filter(TERRITORY1 %ni% c("Antarctica", "Argentina"))

met_sampling <- read.csv("metadata/Metadata_eDNA_Pole2Pole.csv", header=T, sep=";")

met_ant <- met_sampling %>%
  filter(station %in% antarctic$station) %>%
  select(station, latitude_start, longitude_start)


# coordinates of the center of the Coral Triangle
ushuaia <- data.frame(longitude=-67.596442, latitude=-54.932768)


# calculate distance
met_ant$dist_pop <- pointDistance(met_ant[,c("longitude_start", "latitude_start")], ushuaia, lonlat=TRUE)


# transform in KM with 2 decimals
met_ant$dist_pop <- met_ant$dist_pop/1000

met_ant$dist_pop<- round(met_ant$dist_pop, 2) 

met_ant$neartt <- met_ant$dist_pop * 19

antarctic <- left_join(antarctic, met_ant[,c("station", "neartt")])
antarctic <- antarctic %>%
  distinct(station, .keep_all=T)
met_socio <- bind_rows(antarctic, no_antarctic)

write.csv(met_socio, "metadata/eDNA_all_socio_variables.csv", row.names = F)
