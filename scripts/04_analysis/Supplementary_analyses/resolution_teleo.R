library(tidyverse)
library(ggplot2)
library(forcats)
library(ggpubr)
'%ni%' <- Negate("%in%")

load("Rdata/03-filter-data.Rdata")

BE <- read.csv("data/resolution_Exclusive_taxonomic_resolution_BE_by_marker_bdr.csv", header = T)
BS <- read.csv("data/resolution_Taxonomic_resolution_BS_by_marker_bdr.csv", header = T)
all <- read.csv("data/sequencing_by_marker_bdr.csv")

family <- df_filtered %>%
  distinct(family_name_corrected)
family <- as.character(family$family_name_corrected)

BS <- BS %>%
  select(family_name, teleo_mitofish)
colnames(BS) <- c("Family", "resolution")

all <- all %>%
  select(family_name, n_sp_total, teleo_mitofish)
colnames(all) <- c("Family", "n_sp_total", "sequenced")

BS_seq <- left_join(BS, all)

BS_seq <- BS_seq %>%
  filter(Family %in% family) %>%
  filter(Family %ni% c("Cichlidae", "Cyprinidae", "Nemacheilidae", "Poeciliidae", "Salmonidae"))

BS_seq <- BS_seq %>%
  filter(!is.na(n_sp_total)) %>%
  filter(!is.na(resolution))

BS_seq$resolution <- BS_seq$resolution*100
BS_seq$sequenced <- BS_seq$sequenced*100

BS_seq1 <- BS_seq[1:60,]
BS_seq2 <- BS_seq[61:120,]
BS_seq3 <- BS_seq[121:163,]

seg1 <- data.frame(yend=numeric(60))
for (i in 1:60) {
  seg1[i,"yend"] <- max(BS_seq1[i, "resolution"], BS_seq1[i, "sequenced"])
}

seg2 <- data.frame(yend=numeric(60))
for (i in 1:60) {
  seg2[i,"yend"] <- max(BS_seq2[i, "resolution"], BS_seq2[i, "sequenced"])
}

seg3 <- data.frame(yend=numeric(43))
for (i in 1:43) {
  seg3[i,"yend"] <- max(BS_seq3[i, "resolution"], BS_seq3[i, "sequenced"])
}


plot1 <- ggplot(BS_seq1)+
  geom_segment(aes(x=fct_rev(reorder(Family,Family)), y=0, xend=fct_rev(reorder(Family,Family)), yend=seg1$yend), color="grey50")+
  geom_point(aes(x=fct_rev(reorder(Family,Family)), y=sequenced), col="#4587A1", size=2)+
  geom_point(aes(x=fct_rev(reorder(Family,Family)), y=resolution), col="#FFA500", size=2)+
  annotate(geom = "text", x=BS_seq1$Family, y=115, label=BS_seq1$n_sp_total, hjust=1, size=3)+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100))+
  xlab("")+
  ylab("% of species sequenced and resolution")+
  coord_flip()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9))

plot2 <- ggplot(BS_seq2)+
  geom_segment(aes(x=fct_rev(reorder(Family,Family)), y=0, xend=fct_rev(reorder(Family,Family)), yend=seg2$yend), color="grey50")+
  geom_point(aes(x=fct_rev(reorder(Family,Family)), y=sequenced), col="#4587A1", size=2)+
  geom_point(aes(x=fct_rev(reorder(Family,Family)), y=resolution), col="#FFA500", size=2)+
  annotate(geom = "text", x=BS_seq2$Family, y=115, label=BS_seq2$n_sp_total, hjust=1, size=3)+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100))+
  xlab("")+
  ylab("% of species sequenced and resolution")+
  coord_flip()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9))

plot3 <- ggplot(BS_seq3)+
  geom_segment(aes(x=fct_rev(reorder(Family,Family)), y=0, xend=fct_rev(reorder(Family,Family)), yend=seg3$yend), color="grey50")+
  geom_point(aes(x=fct_rev(reorder(Family,Family)), y=sequenced), col="#4587A1", size=2)+
  geom_point(aes(x=fct_rev(reorder(Family,Family)), y=resolution), col="#FFA500", size=2)+
  annotate(geom = "text", x=BS_seq3$Family, y=115, label=BS_seq3$n_sp_total, hjust=1, size=3)+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100))+
  xlab("")+
  ylab("% of species sequenced and resolution")+
  coord_flip()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9))


resolution1 <- ggarrange(plot1, plot2, ncol=2)
ggsave(resolution1, file="outputs/maps & plot sup/resolution1.png", width=7.5, height = 7.5)

resolution2 <- plot3
ggsave(resolution2, file="outputs/maps & plot sup/resolution2.png", width=3.8, height = 6)
