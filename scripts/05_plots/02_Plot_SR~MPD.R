library(ggplot2)
library(tidyverse)
library(ggpubr)

load("Rdata/richness_station.rdata")
load("Rdata/MNTD_station.rdata")

mntd_rich_station <- left_join(rich_station, mntd_stations)

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")

mntd_rich_station <- mntd_rich_station %>%
  left_join(., meta[,c("station", "province")]) %>%
  distinct(station, .keep_all=T)


SR_MNTD <- ggplot(mntd_rich_station, aes(x=MNTD, y=MOTUs)) +
  geom_point(size=2, aes(col=province), show.legend = T)+
  scale_fill_brewer(palette="Paired", direction = 1, aesthetics = "col")+
  labs(x="ses.MNTD within station", y="MOTUs richness")+
  theme_bw()

ggsave("outputs/SR~ses.MNTD.png", width = 8, height = 5)


### assemble with map
load("Rdata/map_global.rdata")


ggarrange(map_global, SR_MNTD, nrow=2, ncol=1, labels = c("A", "B"))
ggsave("outputs/fig1.png", width = 7, height = 7)
