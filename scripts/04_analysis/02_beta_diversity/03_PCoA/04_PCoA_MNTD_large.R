library(dplyr)
library(stringr)
library(rsq)
library(margins)
library(reshape)
library(tidyverse)
library(vegan)
library(ggplot2)
library(patchwork)
library(ggalt)
library(ggrepel)
library(mctest)
library(ggpubr)
library(ade4)
library(ape)

load("Rdata/MNTD_largefish_pairwise_station.rdata")

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(mntd_large),]
identical(as.character(rownames(meta)), rownames(mntd_large))
province <- meta$province



#### PCoA ####

mntd_large <- as.dist(mntd_large)
pcoa <- dudi.pco(mntd_large)
summary(pcoa)
CAP1=23.6
CAP2=18.6

plot <- pcoa$li 
plot$province <- province


ggplot(plot, aes(x = A1, y=A2))+
  geom_point()+
  geom_hline(yintercept = 0, lty = 2, col = "grey", show.legend = F) +
  geom_vline(xintercept = 0, lty = 2, col = "grey", show.legend = F) +
  geom_encircle(aes(group=province, fill=province),
                s_shape = 1, expand = 0, alpha = 0.4, show.legend = TRUE)+
  scale_fill_brewer(palette="Paired", direction = 1, aesthetics = "fill") +
  labs(x = paste0("PCoA1 (", CAP1, "%)"), y = paste0("PCoA2 (", CAP2, "%)"))+
  ggtitle("PCoA - MNTD large fish")+
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1, fill="white")) 
