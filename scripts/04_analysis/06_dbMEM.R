library(codep)
library(adespatial)
library(adegraphics)
library(vegan)
library(car)
library(dplyr)
library(data.table)
library(ggplot2)
library(sf)
library(tidyr)
load("Rdata/MNTD_pairwise_station.rdata")

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(mntd),]
coor <- meta[,c("longitude_start", "latitude_start")]

plot(coor, asp=1)

DistSpatial=gcd.hf(coor)

dbmem = dbmem(DistSpatial)
rownames(dbmem) <- labels(DistSpatial)

summary(dbmem)

adegraphics::s.label(coor, nb = attr(dbmem, "listw"))

ade4::s.value(coor, dbmem[,1])

dbmem_coor <- cbind(coor, dbmem)

save(dbmem_coor, file="Rdata/db_mem.rdata")

