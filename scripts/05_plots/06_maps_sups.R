library(vegan)
library(tidyverse)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggthemes)
library(ggpubr)

# load richness data
load("Rdata/richness_station.rdata")
rownames(rich_station) <- rich_station$station

#load untransformed variables
load("Rdata/geographic_variables.rdata")
load("Rdata/socioeconomic_variables.rdata")
load("Rdata/environmental_variables.rdata")
load("Rdata/sampling_variables.rdata")

# load transformed variables
load("Rdata/all_explanatory_variables_numeric.rdata")
# Formate data
rownames(env_var) <- env_var$station
env_var <- env_var[rownames(rich_station),]
env_var <- env_var[, colSums(env_var != 0) > 0]

rownames(geo_var) <- geo_var$station
geo_var <- geo_var[rownames(rich_station),]

rownames(samp_var) <- samp_var$station
samp_var <- samp_var[rownames(rich_station),]

rownames(socio_var) <- socio_var$station
socio_var <- socio_var[rownames(rich_station),]

exp_var <- cbind(env_var, socio_var[,-1], geo_var[,-1], samp_var[,-1])


# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
data <- exp_var

coor <- meta[, c("longitude_start", "latitude_start")]
identical(as.character(rownames(coor)), as.character(rownames(exp_var_num)))
data2 <- cbind(exp_var_num, coor)

# select variables

data <- left_join(data, rich_station, by="station")
data2 <- left_join(data2, rich_station, by="station")
data2$MOTUs <- log10(data2$MOTUs +1)
data2$crypto_MOTUs <- log10(data2$crypto_MOTUs +1)
data2$largefish_MOTUs <- log10(data2$largefish_MOTUs +1)

#### Maps ####
world <- ne_countries(type="map_units",returnclass = 'sf')
world2 <- st_read("c://Users/mathon/Desktop/PhD/Projets/Megafauna/Carto_megafauna/GSHHS_f_L1.shp")
antarctica <- st_read("c://Users/mathon/Desktop/PhD/Projets/Megafauna/Carto_megafauna/GSHHS_f_L6.shp")
antarctica2 <- st_read("c://Users/mathon/Desktop/PhD/Projets/Megafauna/Carto_megafauna/Couche_ile_Antarctique.shp")

# Convert point in dataset to sf object
data$metadata_map_sf = st_as_sf(data[,c("longitude_start", "latitude_start")], coords = c("longitude_start", "latitude_start"), 
                           crs = 4326)

data2$metadata_map_sf = st_as_sf(data2[,c("longitude_start", "latitude_start")], coords = c("longitude_start", "latitude_start"), 
                                crs = 4326)

# MOTUs richness
data2 <- data2[order(data2$MOTUs),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = data2$MOTUs), size=2, data= data2$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(xlim = c(-180, 180), ylim = c(-80, 90)) +
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="log(MOTUs richness) ")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(3, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = "grey50"),
        plot.title = element_text(lineheight=.8, face="bold"))

ggsave("outputs/maps & plot sup/log(MOTUs_richness).png")

# Crypto richness
data2 <- data2[order(data2$crypto_MOTUs),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = data2$crypto_MOTUs), size=2, data= data2$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(xlim = c(-180, 180), ylim = c(-80, 90)) + 
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="log(Cryptobenthics richness) ")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(3, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = "grey50"),
        plot.title = element_text(lineheight=.8, face="bold"))

ggsave("outputs/maps & plot sup/log(crypto_richness).png")


# large fish richness
data2 <- data2[order(data2$largefish_MOTUs),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = data2$largefish_MOTUs), size=2, data= data2$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(xlim = c(-180, 180), ylim = c(-80, 90)) + 
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="log(large fish richness) ")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(3, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = "grey50"),
        plot.title = element_text(lineheight=.8, face="bold"))

ggsave("outputs/maps & plot sup/log(large_richness).png")


# ses.MNTD
load("Rdata/MNTD_station.rdata")
data <- left_join(data, mntd_stations)
data <- data[order(data$MNTD),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = data$MNTD), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(xlim = c(-180, 180), ylim = c(-80, 90)) + 
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="ses.MNTD")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(3, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = "grey50"),
        plot.title = element_text(lineheight=.8, face="bold"))

ggsave("outputs/maps & plot sup/ses.MNTD.png")


#### map for explanatory variables (change variable and title) by region ####
data <- data[order(data$MarineEcosystemDependency),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") +
  geom_sf(aes(fill = data$MarineEcosystemDependency), size=3, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(xlim = c(-180, 180), ylim = c(-80, 90)) +
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="Marine Ecosystem Dependency ")+
  #ggtitle("Tropical Southwestern Pacific")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(3, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = "grey95"),
        plot.title = element_text(lineheight=.8, face="bold"))

ggsave("outputs/maps & plot sup/MED.png")


#### map for explanatory variables (change variable and title) by region ####
data2 <- data2[order(data2$Gravity),]

arctic_grav <- ggplot() + 
  geom_sf(aes(), data = world2, fill = "grey80", col="grey80") +
  geom_sf(aes(fill = data2$Gravity), size=4, data= data2$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(xlim = c(3, 34), ylim = c(75, 82)) +
  scale_fill_gradient2(low="yellow", high="red", mid="orange", midpoint = 0.2, aesthetics = "fill")+
  labs(fill="log(Gravity+1) ")+
  ggtitle("Arctic")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(3, "mm"),
        #axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = "grey95"),
        plot.title = element_text(lineheight=.8, face="bold"))


ggarrange(arctic_MED, Atl_Med_MED, china_MED, leng_MED, nc_MED, eparse_MED, cari_MED, ant_MED,
          ncol=2, nrow=4, common.legend = T, legend = c("bottom"))
ggsave("outputs/maps & plot sup/map_MED_region.png", height = 16, width = 10)


ggarrange(arctic_grav, Atl_Med_grav, china_grav, leng_grav, nc_grav, eparse_grav, cari_grav, ant_grav,
          ncol=2, nrow=4, common.legend = T, legend = c("bottom"))
ggsave("outputs/maps & plot sup/map_grav_region.png", height = 16, width = 10)
