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
library(grid)
library(spdep)
library(mctest)
library(ggpubr)
source("scripts/04_analysis/00_functions.R")

load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
load("Rdata/MNTD_pairwise_station.rdata")
load("Rdata/db_mem.rdata")

# transform data

data <- exp_var
df <- exp_var_num %>%
  select(-c("station")) # remove station
df_mem <- cbind(df, dbmem)

#---------------------------------------------------------------------------------------------------------------------------
#### Full model ####

dbrda_full <- capscale(mntd ~ .,df)

RsquareAdj(dbrda_full)
anova(dbrda_full)
anova(dbrda_full, by = "term", permutations = 99)



# check for colinearity and select variables
mctest::imcdiag(dbrda_full, method="VIF")

df_sel <- df%>%
  select(-c("pH_mean", "NGO"))

#### dbRDA with selected variables ####
dbrda_sel <- capscale(mntd ~ .,df_sel)
                       
RsquareAdj(dbrda_sel)
anova(dbrda_sel)
anova(dbrda_sel, by = "term", permutations = 99)

mctest::imcdiag(dbrda_sel, method="VIF")

#---------------------------------------------------------------------------------------------------------------------------
#### dbrda total with all selected variables, correcting for spatial ####

dbrda_fin <- capscale(mntd ~ mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+neartt+Gravity+MarineEcosystemDependency+Naturalresourcesrents+dist_to_CT+bathy+depth_sampling+latitude+distCoast+sample_method+sequencer+volume+mean_DHW_1year + Condition(MEM1+MEM2+MEM3+MEM4+MEM5), df_sel)
RsquareAdj(dbrda_fin)
anova(dbrda_fin)
anova(dbrda_fin, by = "term", permutations = 99)
AIC(dbrda_fin)


# var part (need to remove almost constant variables i.e. filter, sequencing_depth)
env_var <- df_sel[,c("mean_sss_1year", "mean_npp_1year", "mean_SST_1year", "mean_DHW_1year", "mean_DHW_5year")]
geo_var <- df_sel[, c("distCoast", "latitude", "bathy", "dist_to_CT", "depth_sampling")]
socio_var <- df_sel[,c("HDI2019", "neartt", "Gravity", "MarineEcosystemDependency", "Naturalresourcesrents")]
samp_var <- df[, c("volume")]
mntd <- as.dist(mntd)

varpart_tot <- varpart(mntd, env_var, geo_var, socio_var, samp_var)
varpart_tot

varplot <- plot(varpart_tot, digits = 2, Xnames = c('env_var', 'geo_var', 'socio_var', 'samp_var'), bg = c('navy', 'tomato', 'yellow', 'green'))

# get scores
station_scores <- scores(dbrda_sel)$sites
var_scores <- dbrda_sel[["CCA"]][["biplot"]][, 1:2] %>% as.data.frame()

# get most differentiated species along first axis
quant75_cap1 <- quantile(abs(var_scores$CAP1), probs = c(0.75))
quant75_cap2 <- quantile(abs(var_scores$CAP2), probs = c(0.75))
quant75 <- rbind(quant75_cap1, quant75_cap2 )
var_scores_diff75_cap1 <- var_scores[which(abs(var_scores$CAP1) > quant75_cap1["75%"]),]
var_scores_diff75_cap2 <- var_scores[which(abs(var_scores$CAP2) > quant75_cap2["75%"]),]
var_scores_diff75 <- rbind(var_scores_diff75_cap1, var_scores_diff75_cap2)
var_scores_diff75 <- unique(var_scores_diff75)

# extract the percentage variability explained by axes
sumdbrda <- summary(dbrda_sel)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(rownames(df)), rownames(station_scores)) # verify that data in same order
station_scores_met <- cbind(station_scores, data)

grda_station <- ggplot(station_scores_met, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_point(col = "black", cex = 1) +
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

grda_variables <- ggplot() + 
  geom_segment(data= var_scores, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "grey",
               arrow=arrow(length=unit(0.01,"npc"))) + # all species
  geom_segment(data= var_scores_diff75, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc"))) + # most differentiated species
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_label_repel(data= var_scores_diff75, 
                   aes(x= CAP1, y=CAP2, #hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2)),
                       fontface=3, size = 3),
                   label = rownames(var_scores_diff75),
                   show.legend = F) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position = c(0, 1),              # position in top left corner
        legend.justification = c(0, 1),         # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),   # add margin as to not overlap with axis box
        legend.background = element_rect(fill =  alpha("white", 0.0)),
        legend.title = element_text(size=11),
        legend.text = element_text(size=11)) +
  # re-add legend and change text legend key by making invisible points and overriding its key shape
  geom_point(data= var_scores_diff75, 
             aes(x=CAP1, y=CAP2),
             size = 0, stroke = 0) + 
  guides(colour = guide_legend(override.aes = list(size = 5, shape = c(utf8ToInt("C"), utf8ToInt("B"), utf8ToInt("D"), utf8ToInt("P")))))
grda_variables

ggarrange(grda_station, grda_variables, nrow=2)
ggsave("outputs/dbRDA/MNTD/dbrda_tot.png", width = 5, height = 8)
#---------------------------------------------------------------------------------------------------------------------------------
#### partial dbrda correcting for sampling and spatial ####

dbrda_part <- capscale(mntd ~ mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+neartt+Gravity+MarineEcosystemDependency+Naturalresourcesrents+dist_to_CT+bathy+depth_sampling+latitude+distCoast+mean_DHW_1year + Condition(volume), df_sel) #MEM1+MEM2+MEM3+MEM4+MEM5+
RsquareAdj(dbrda_part)
anova(dbrda_part)
anova(dbrda_part, by = "term", permutations = 99)
anova(dbrda_part, by = "margin", permutations = 99)


# variation partitioning
varpart_part <- varpart(mntd, env_var, geo_var, socio_var, samp_var)
varpart_part

plot(varpart_part, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# get scores
station_scores <- scores(dbrda_part)$sites
var_scores <- dbrda_part[["CCA"]][["biplot"]][, 1:2] %>% as.data.frame()

# get most differentiated species along first axis
quant75_cap1 <- quantile(abs(var_scores$CAP1), probs = c(0.75))
quant75_cap2 <- quantile(abs(var_scores$CAP2), probs = c(0.75))
quant75 <- rbind(quant75_cap1, quant75_cap2 )
var_scores_diff75_cap1 <- var_scores[which(abs(var_scores$CAP1) > quant75_cap1["75%"]),]
var_scores_diff75_cap2 <- var_scores[which(abs(var_scores$CAP2) > quant75_cap2["75%"]),]
var_scores_diff75 <- rbind(var_scores_diff75_cap1, var_scores_diff75_cap2)
var_scores_diff75 <- unique(var_scores_diff75)

# extract the percentage variability explained by axes
sumdbrda <- summary(dbrda_part)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(rownames(df)), rownames(station_scores)) # verify that data in same order
station_scores_met <- cbind(station_scores, data)

grda_station <- ggplot(station_scores_met, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_encircle(aes(group = province, fill= province), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = TRUE) + # hull area 
  geom_point(col = "black", cex = 1) +
  scale_fill_brewer(palette="Paired", direction = 1, aesthetics = "fill") +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "none",             # position in top left corner
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_text(size=11),
        legend.text = element_text(size=11),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_station

grda_variables <- ggplot() + 
  geom_segment(data= var_scores, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "grey",
               arrow=arrow(length=unit(0.01,"npc"))) + # all species
  geom_segment(data= var_scores_diff75, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc"))) + # most differentiated species
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_label_repel(data= var_scores_diff75, 
                   aes(x= CAP1, y=CAP2, #hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2)),
                       fontface=3, size = 3),
                   label = rownames(var_scores_diff75),
                   show.legend = F) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position = c(0, 1),              # position in top left corner
        legend.justification = c(0, 1),         # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),   # add margin as to not overlap with axis box
        legend.background = element_rect(fill =  alpha("white", 0.0)),
        legend.title = element_text(size=11),
        legend.text = element_text(size=11)) +
  # re-add legend and change text legend key by making invisible points and overriding its key shape
  geom_point(data= var_scores_diff75, 
             aes(x=CAP1, y=CAP2),
             size = 0, stroke = 0) + 
  guides(colour = guide_legend(override.aes = list(size = 5, shape = c(utf8ToInt("C"), utf8ToInt("B"), utf8ToInt("D"), utf8ToInt("P")))))
grda_variables

ggarrange(grda_station, grda_variables, nrow=2, common.legend = TRUE, legend = "right")
ggsave("outputs/dbRDA/MNTD/dbrda_part.png", width = 8, height = 8)


