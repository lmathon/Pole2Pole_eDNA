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

# change from wide to long format
dbmem_gps_long <- gather(dbmem_coor, MEM, Value, MEM1:ncol(dbmem_coor))
# Calculate an MEM average value for each GPS point.
library(stringr)
dbmem_gps <- dbmem_gps_long %>% group_by(latitude_start, longitude_start, MEM)%>%
  summarise(mem_mean <- mean(Value))
setnames(dbmem_gps, "mem_mean <- mean(Value)", "Average_MEM")
dbmem_wide <- spread(dbmem_gps, MEM, Average_MEM)


#  Download a high resolution map with the sf package
library(maps)
library(mapdata)

wH <- map_data("worldHires",  xlim=c(-179,179), ylim=c(-90,90)) # subset polygons surrounding med sea

ggplot() +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = NA) +
  coord_fixed(xlim=c(-179,179), ylim=c(-90,90), ratio = 1.3) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())

dbmem_gps_long_mem14 <- subset(dbmem_gps, subset=MEM=="MEM1"| MEM=="MEM4")


pdf("MEM1_4_fas.pdf", width=5, height=5)
x_title="Longitude"
y_title="Latitude"
graph1 <- ggplot() +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = NA) +
  coord_fixed(xlim = c(-179,179), ylim=c(-90,90))+
  facet_wrap(~MEM)+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())+
  geom_point(aes(x = longitude_start, y = latitude_start,fill=Average_MEM), data=dbmem_gps_long_mem14,shape = 21, size=1.5)+
  theme_bw()+theme(legend.position = "none",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y=y_title)+  
  labs(x=x_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))+
  scale_fill_continuous(high="yellow",low="red")
graph1
dev.off()
