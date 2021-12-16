library(ggplot2)
library(tidyverse)
library(ggpubr)

load("Rdata/richness_station.rdata")
load("Rdata/FD_Hill_alpha.rdata")

FD_rich_station <- left_join(rich_station, FD_Hill)

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")

FD_rich_station <- FD_rich_station %>%
  left_join(., meta[,c("station", "province")]) %>%
  distinct(station, .keep_all=T)


SR_FDq2 <- ggplot(FD_rich_station, aes(x=MOTUs, y=FD_q2)) +
  geom_point(size=2, aes(col=province), show.legend = T)+
  scale_fill_brewer(palette="Paired", direction = 1, aesthetics = "col")+
  labs(x="MOTUs richness", y= expression(paste("Sequence ", alpha,"-diversity")))+
  theme_bw()

ggsave("outputs/SR~FDq2.png", width = 8, height = 5)


### assemble with map
load("Rdata/map_global.rdata")


ggarrange(map_global, SR_FDq2, nrow=2, ncol=1, labels = c("A", "B"))
ggsave("outputs/fig1.png", width = 7, height = 7)

cor.test(FD_rich_station$FD_q2, FD_rich_station$MOTUs, method = "spearman")
