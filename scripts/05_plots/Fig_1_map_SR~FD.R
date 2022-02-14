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


SR_FDq0 <- ggplot(FD_rich_station, aes(x=MOTUs, y=FD_q0)) +
  geom_point(size=2, aes(col=province), show.legend = T)+
  labs(color = "Region")+
  scale_fill_manual(values=c("#A6CEE3","#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#B15928", "#FF7F00", "#CAB2D6","#6A3D9A","#FFD92F"), aesthetics = "col")+
  labs(x=expression(paste("Taxonomic ", alpha,"-diversity")), y= expression(paste("Sequence ", alpha,"-diversity")))+
  theme_bw()

ggsave("outputs/SR~FDq0.png", width = 8, height = 5)


### assemble with map
load("Rdata/map_global.rdata")


ggarrange(map_global, SR_FDq0, nrow=2, ncol=1, labels = c("A", "B"))
ggsave("outputs/Figures_papier/Fig1.pdf", width = 7.2, height = 7.2, dpi = 600)

cor.test(FD_rich_station$FD_q0, FD_rich_station$MOTUs, method = "pearson")
