library(vegan)
library(tidyverse)
library(ggplot2)
library(ade4)

load("Rdata/MNTD_pairwise_station.rdata")
mntd <- as.dist(mntd)


# run pcoa and select axis
pcoa <- dudi.pco(mntd)

pcoa_axis <- pcoa$li

# join metadata
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(pcoa_axis),]
identical(as.character(rownames(meta)), rownames(pcoa_axis))
met <- meta[, c("province", "longitude_start", "latitude_start")]
data <- cbind(pcoa_axis, met)

ggplot(data, aes(x=A1, y=A2))+
  geom_point(aes(x=A1, y=A2))+
  geom_encircle(aes(group = province, fill= province), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = TRUE) +
  scale_fill_brewer(palette="Paired", direction = 1, aesthetics = "fill") +
  xlab("PCoA1")+
  ylab("PCoA2")+
  #xlim(-0.6, 0.4)+
  #ylim(-0.5, 0.5)+
  theme(legend.position = c("right"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_text(size=11),
        axis.title.x = element_text(size=11),
        panel.border = element_rect(fill = NA))
