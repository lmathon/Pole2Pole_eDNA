library(vegan)
library(tidyverse)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggpubr)

# load richness data
load("Rdata/richness_station.rdata")
rownames(rich_station) <- rich_station$station


# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- meta[, c("longitude_start", "latitude_start")]
data <- cbind(rich_station, coor)

#### Maps ####
world <- ne_countries(type="map_units",returnclass = 'sf')

# Convert point in dataset to sf object
data$metadata_map_sf = st_as_sf(data[,c("longitude_start", "latitude_start")], coords = c("longitude_start", "latitude_start"), 
                                crs = 4326)

# MOTUs richness
data <- data[order(data$MOTUs),]

MOTUs <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = data$MOTUs), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="MOTUs richness")+
  ggtitle("36% of stations with < 20 MOTUs")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=25),
        legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(4, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(size = 10, face="bold", hjust = 0.5))

# Chondri richness
data <- data[order(data$chondri_MOTUs),]
data$chondri_MOTUs[data$chondri_MOTUs == 0] <- NA

chondri <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = data$chondri_MOTUs), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="yellow", high="red", na.value = "black", aesthetics = "fill")+
  labs(fill="Chondrichthyan richness")+
  ggtitle("46% of stations with 0 Chondrichthyan MOTUs")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=25),
        legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(4, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(size = 10, face="bold", hjust = 0.5))

# Crypto richness
data <- data[order(data$crypto_MOTUs),]

crypto <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = data$crypto_MOTUs), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="Cryptobenthics richness")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=25),
        legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(4, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(lineheight=.8, face="bold"))

# Large fish richness
data <- data[order(data$largefish_MOTUs),]
data$largefish_MOTUs[data$largefish_MOTUs == 0] <- NA

large <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = data$largefish_MOTUs), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="yellow", high="red", na.value = "black", aesthetics = "fill")+
  labs(fill="Large fish richness")+
  ggtitle("16% of stations with 0 large fish MOTUs")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=25),
        legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(4, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(size = 10, face="bold", hjust = 0.5))


ggarrange(MOTUs, crypto, chondri, large, nrow=2, ncol=2, labels = c("A", "B", "C", "D"))
ggsave("outputs/maps & plot sup/Figure1_maps_richness.png", width = 8, height = 6)
