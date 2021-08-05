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
library(stringr)
library(maps)
library(mapdata)


load("Rdata/MNTD_pairwise_station.rdata")

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(mntd),]
coor <- meta[,c("latitude_start", "longitude_start")]


DistSpatial=gcd.hf(coor)

dbmem = dbmem(DistSpatial)

summary(dbmem)

adegraphics::s.label(coor, nb = attr(dbmem, "listw"))

ade4::s.value(coor, dbmem[,13])

attributes(dbmem)$values


# check with dbrda wich MEM explain the most spatial autocorrelation
dbrda0 <- capscale(mntd ~ 1, dbmem[,1:15])
dbrdaG <- capscale(mntd ~ ., dbmem[,1:15])
mem_sel <- ordiR2step(dbrda0, scope = formula(dbrdaG), direction="both")

RsquareAdj(mem_sel)

anova(mem_sel)
mem_sel$anova
anova(mem_sel, by = "axis",  permutations = 99)
anova(mem_sel, by = "margin", permutations = 99)


dbmem <- dbmem[,1:5]


save(dbmem, file="Rdata/db_mem.rdata")


