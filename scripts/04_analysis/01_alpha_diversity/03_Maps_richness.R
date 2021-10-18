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
df1 <- data[order(data$MOTUs),]

MOTUs <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = df1$MOTUs), size=2, data= df1$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="MOTUs richness")+
  #ggtitle("36% of stations with < 20 MOTUs")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.key.height = unit(3, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(size = 10, face="bold", hjust = 0.5))

# Chondri richness
df2 <- data[order(data$chondri_MOTUs),]
df2 <- df2 %>%
  filter(chondri_MOTUs >0) 

chondri <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = df2$chondri_MOTUs), size=2, data= df2$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="yellow", high="red", na.value = "black", aesthetics = "fill")+
  labs(fill="Chondrichthyan richness")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.key.height = unit(3, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(size = 10, face="bold", hjust = 0.5))

# 0 Chondri 
df3 <- data[order(data$chondri_MOTUs),]
df3 <- df3 %>%
  filter(chondri_MOTUs == 0) 

chondri0 <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(col="grey20", fill = "grey20", size=2, data= df3$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  xlab("46% of stations with 0 Chondrichthyan MOTUs")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(size=12, face = "bold"),
        panel.grid.major = element_line(colour = "gray70"))

# Crypto richness
df4 <- data[order(data$crypto_MOTUs),]

crypto <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = df4$crypto_MOTUs), size=2, data= df4$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="Cryptobenthics richness")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.key.height = unit(3, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(lineheight=.8, face="bold"))

# Large fish richness
df5 <- data[order(data$largefish_MOTUs),]
df5 <- df5 %>%
  filter(largefish_MOTUs> 0)

large <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(aes(fill = df5$largefish_MOTUs), size=2, data= df5$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="yellow", high="red", aesthetics = "fill")+
  labs(fill="Large fish richness")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.key.height = unit(3, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(size = 10, face="bold", hjust = 0.5))

# 0 Large fish 
df6 <- data[order(data$largefish_MOTUs),]
df6 <- df6 %>%
  filter(largefish_MOTUs == 0)

large0 <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") + 
  geom_sf(col="grey20", fill = "grey20", size=2, data= df6$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  xlab("16% of stations with 0 large fish MOTUs")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = "gray70"),
        axis.title.x = element_text(size = 12, face="bold"))


ggarrange(MOTUs, crypto, large, large0, nrow=2, ncol=2, labels = c("A", "B", "C", "D"), font.label = list(size=12))
ggsave("outputs/maps & plot sup/Figure1_maps_richness.png", width = 8, height = 6)

ggarrange(chondri, chondri0, nrow = 1, ncol=2,labels = c("A", "B"), font.label = list(size=12))
ggsave("outputs/maps & plot sup/chondri_richness_null.png", width = 8, height = 3)
