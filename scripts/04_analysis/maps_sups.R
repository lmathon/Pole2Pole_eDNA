library(vegan)
library(tidyverse)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# load richness data
load("Rdata/richness_station.rdata")
rownames(rich_station) <- rich_station$station

#load untransformed variables
load("Rdata/geographic_variables.rdata")
load("Rdata/socioeconomic_variables.rdata")
load("Rdata/environmental_variables.rdata")
load("Rdata/sampling_variables.rdata")

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
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- meta[, c("longitude_start", "latitude_start")]
data <- cbind(exp_var, coor)


# select variables

data <- left_join(data, rich_station, by="station")



# plot dist_to_CT ~ Gravity
ggplot(data, aes(Gravity, dist_to_CT))+
  geom_bin2d()+
  labs(fill="Nb of stations")+
  theme(panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())



#### Maps ####
world <- ne_countries(type="map_units",returnclass = 'sf')

# Convert point in dataset to sf object
data$metadata_map_sf = st_as_sf(data[,c("longitude_start", "latitude_start")], coords = c("longitude_start", "latitude_start"), 
                           crs = 4326)

# MOTUs richness
data <- data[order(data$MOTUs),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(aes(fill = data$MOTUs), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="blue", high="red", aesthetics = "fill")+
  labs(fill="MOTUs richness")+
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

# Family richness
data <- data[order(data$Family),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(aes(fill = data$Family), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="blue", high="red", aesthetics = "fill")+
  labs(fill="Family richness")+
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

# Chondri richness
data <- data[order(data$chondri_MOTUs),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(aes(fill = data$chondri_MOTUs), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="blue", high="red", aesthetics = "fill")+
  labs(fill="Chondrichthyans richness")+
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

# Crypto richness
data <- data[order(data$crypto_MOTUs),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(aes(fill = data$crypto_MOTUs), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="blue", high="red", aesthetics = "fill")+
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

# MPD
load("Rdata/MPD_station.rdata")
data <- left_join(data, mpd)
data <- data[order(data$mpd),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(aes(fill = data$mpd), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="blue", high="red", aesthetics = "fill")+
  labs(fill="Mean Pairwise genetic distance")+
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


# Less than 20 MOTUs
data_20 <- data %>%
  filter(MOTUs < 21)

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(fill="blue", size=2, data= data_20$metadata_map_sf, shape=21) + 
  coord_sf(crs = "+proj=robin") + 
  ggtitle("Stations with less than 20 MOTUs")+
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
        plot.title = element_text(size=12, face="bold", hjust = 0.5))

# Zero chondrichthyans
data_zero_chondri <- data %>%
  filter(chondri_MOTUs == 0)

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(fill="blue", size=2, data= data_zero_chondri$metadata_map_sf, shape=21) + 
  coord_sf(crs = "+proj=robin") + 
  ggtitle("Stations with no Chondrichthyans")+
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
        plot.title = element_text(size=12, face="bold", hjust = 0.5))

#### map for explanatory variables (change variable and title) ####
data <- data[order(data$volume),]

ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(aes(fill = data$volume), size=2, data= data$metadata_map_sf, shape=21, show.legend = T) + 
  coord_sf(crs = "+proj=robin") + 
  scale_fill_gradient(low="blue", high="red", aesthetics = "fill")+
  labs(fill="Volume (L)")+
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
        legend.key.width = unit(8, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(lineheight=.8, face="bold"))
