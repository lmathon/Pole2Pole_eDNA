library(ggplot2)
library(tidyverse)
library(ggpubr)

load("Rdata/richness_station.rdata")
load("Rdata/MPD_station.rdata")

mpd_rich_station <- left_join(rich_station, mpd_stations)

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")

mpd_rich_station <- left_join(mpd_rich_station, meta[,c("station", "province")])


SR_MPD <- ggplot(mpd_rich_station, aes(x=MPD, y=MOTUs)) +
  geom_point(size=2, aes(col=province), show.legend = F)+
  scale_fill_brewer(palette="Paired", direction = 1, aesthetics = "col")+
  labs(x="Mean Pairwise Distance within station", y="MOTUs richness")+
  theme_bw()

ggsave("outputs/SR~MPD.png", width = 8, height = 5)


### assemble with map
load("Rdata/map_global.rdata")


ggarrange(map_global, SR_MPD, nrow=2, ncol=1, labels = c("A", "B"))
