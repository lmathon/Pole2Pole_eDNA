library(dplyr)
library(forcats)
library(stringr)
library(rsq)
library(margins)
library(betapart)
library(reshape)
library(tidyverse)
library(tidyselect)
library(vegan)
library(ggplot2)
library(patchwork)
library(ggalt)
library(ggrepel)
library(grid)
source("scripts/04_analysis/00_functions.R")

load("Rdata/all_explanatory_variables.rdata")
load("Rdata/MNTD_pairwise_station.rdata")
load("Rdata/geographic_distance_stations.rdata")

# transform data
mntd <- as.matrix(mntd)
dist_km <- dist_km[rownames(mntd), colnames(mntd)]

df <- exp_var %>% filter(station %in% rownames(mntd))
rownames(df) <- df$station
df <- df[rownames(mntd), ]
df <- df[,-62]
df <- df[, colSums(df != 0) > 0]

#-----------------------------------------------------------------------------------------------------------------------------
# total dbrda
#-----------------------------------------------------------------------------------------------------------------------------

dbrda_tot <- capscale(mntd ~ mean_DHW_1year+mean_DHW_5year+mean_SST_1year+mean_SST_5year+mean_sss_1year+mean_sss_5year+mean_npp_1year+mean_npp_5year+pH_mean+dist_to_coast+dist_to_CT+province+depth_fin+depth_sampling+latitude_start+sample_type+sample_method+filter+sequencing_depth+volume, df)
RsquareAdj(dbrda_tot)
anova(dbrda_tot)
anova(dbrda_tot, by = "axis",  permutations = 99)
anova(dbrda_tot, by = "term", permutations = 99)

station_scores <- scores(dbrda_tot)$sites
var_scores <- dbrda_tot[["CCA"]][["biplot"]][, 1:2] %>% as.data.frame()

# extract the percentage variability explained by axes
sumdbrda <- summary(dbrda_tot)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(df$station), rownames(station_scores)) # verify that data in same order
station_scores_met <- cbind(station_scores, df)

grda_station <- ggplot(station_scores_met, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_point(col = "black", cex = 1) +
  geom_segment(data= var_scores, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "blue",
               arrow=arrow(length=unit(0.01,"npc")))+
  geom_label_repel(data= var_scores, 
                   aes(x= CAP1, y=CAP2), 
                   fontface=3, size = 3,
                   label = rownames(var_scores),
                   show.legend = F, max.overlaps = 18) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = c(0, 1),             # position in top left corner
        legend.justification = c(0, 1),        # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_text(size=11),
        legend.text = element_text(size=11),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_station

# var part (need to remove almost constant variables i.e. filter, sequencing_depth)
env_var <- as.data.frame(df[,2:47])
geo_var <- df[, 48:53]
samp_var <- df[, c(54,55,58)]
mntd <- as.dist(mntd)

varpart_tot <- VarPart(mntd, env_var, geo_var, samp_var)
varpart_tot

plot(varpart_tot, digits = 2, Xnames = c('env_var', 'geo_var', 'samp_var'), bg = c('navy', 'tomato', 'yellow'))
#--------------------------------------------------------------------------------------------------------------------------------------
# partial dbrda : correcting for sampling variables
#--------------------------------------------------------------------------------------------------------------------------------------

dbrda_part <- capscale(mntd ~ mean_DHW_1year+mean_DHW_5year+mean_SST_1year+mean_SST_5year+mean_sss_1year+mean_sss_5year+mean_npp_1year+mean_npp_5year+pH_mean+dist_to_coast+dist_to_CT+province+depth_fin+depth_sampling+latitude_start + Condition(sample_type,sample_method,filter,sequencing_depth,volume), df)
RsquareAdj(dbrda_part)
anova(dbrda_part)
anova(dbrda_part, by = "axis",  permutations = 99)
anova(dbrda_part, by = "term", permutations = 99)


station_scores <- scores(dbrda_part)$sites
var_scores <- dbrda_part[["CCA"]][["biplot"]][, 1:2] %>% as.data.frame()

# extract the percentage variability explained by axes
sumdbrda <- summary(dbrda_part)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(df$station), rownames(station_scores)) # verify that data in same order
station_scores_met <- cbind(station_scores, df)

grda_station <- ggplot(station_scores_met, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_point(col = "black", cex = 1) +
  geom_segment(data= var_scores, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "blue",
               arrow=arrow(length=unit(0.01,"npc")))+
  geom_label_repel(data= var_scores, 
                   aes(x= CAP1, y=CAP2), 
                       fontface=3, size = 3,
                   label = rownames(var_scores),
                   show.legend = F, max.overlaps = 18) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = c(0, 1),             # position in top left corner
        legend.justification = c(0, 1),        # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_text(size=11),
        legend.text = element_text(size=11),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_station

