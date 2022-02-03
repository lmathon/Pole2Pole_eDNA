library(vegan)
library(tidyverse)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggthemes)
library(ggrepel)
library(RColorBrewer)

load("Rdata/03-filter-data.Rdata")

Sites <- unique(df_filtered$site_name)

data <- data.frame(province=character(), site=character(), longitude=numeric(), latitude=numeric(), n_stations=numeric())

for (i in 1:length(Sites)) {
  df <- df_filtered %>%
    filter(site_name==Sites[i])
  data[i, "province"] <- unique(df$province)
  data[i, "site"] <- Sites[i]
  data[i, "longitude"] <- df[1,"longitude_start"]
  data[i, "latitude"] <- df[1,"latitude_start"]
  data[i, "n_stations"] <- n_distinct(df$station)
  
}

world <- ne_countries(type="map_units",returnclass = 'sf')

data$metadata_map_sf = st_as_sf(data[,c("longitude", "latitude")], coords = c("longitude", "latitude"), 
                                crs = 4326)
map_global <- ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80", col="grey80") +
  geom_text_repel(data = data, aes(x=longitude, y=latitude, label = n_stations),size=3, min.segment.length = 0.2, force = 2, max.overlaps=20) +
  geom_sf(aes(fill = data$province), size=2, data= data$metadata_map_sf, shape=21, show.legend = F) + 
  coord_sf(xlim = c(-180, 180), ylim = c(-80, 90)) +
  scale_fill_manual(values=c("#A6CEE3","#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#B15928", "#FF7F00", "#CAB2D6","#6A3D9A","#FFD92F"), aesthetics = "fill")+
  labs(x="", y="")+
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.height = unit(3, "mm"),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1, fill=NA))

ggsave("outputs/maps & plot sup/map_global_label.png")
save(map_global, file="Rdata/map_global.rdata")
