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
library(spdep)
source("scripts/04_analysis/00_functions.R")

load("Rdata/all_explanatory_variables_numeric.rdata")
load("Rdata/Jaccard_MOTU_dissimilarity.rdata")
load("Rdata/db_mem.rdata")

# transform data

df <- exp_var_num
df <- df[,-16] # remove temporarily incomplete variables (sequencer)
df_mem <- cbind(df[,-1], dbmem)

#---------------------------------------------------------------------------------------------------------------------------
#### Model selection ####

# dbrda on null model, and full model
dbrda_0 <- capscale(dist_jac_mo ~ 1, df_mem)
dbrda_tot <- capscale(dist_jac_mo ~ .,df_mem)

vif(dbrda_tot)

# run selection of best model
dbrda_sel <- ordiR2step(dbrda_0, scope = formula(dbrda_tot), direction="both") 

# check out selected variables 
dbrda_sel$anova

vif(dbrda_sel) # -> select variables : filter, sequencing_depth & latitude cause problems

RsquareAdj(dbrda_sel)
anova(dbrda_sel)
anova(dbrda_sel, by = "term", permutations = 99)


#---------------------------------------------------------------------------------------------------------------------------
#### dbrda final with all selected variables, correcting for spatial ####

sel_vars <- rownames(dbrda_sel$anova)[-nrow(dbrda_sel$anova)] %>% sub("..", "",.)
df_sel <- df_mem %>% select(all_of(sel_vars))

dbrda_fin <- capscale(dist_jac_mo ~ dist_to_coast+mean_sss_1year+depth_fin+sample_method+mean_npp_1year+mean_SST_1year+dist_to_CT+province+volume+depth_sampling+pH_mean+mean_DHW_1year + Condition(MEM1+MEM2+MEM3+MEM4+MEM5), df_sel)
RsquareAdj(dbrda_fin)
anova(dbrda_fin)
anova(dbrda_fin, by = "term", permutations = 99)


# var part (need to remove almost constant variables i.e. filter, sequencing_depth)
env_var <- df_sel[,c("mean_sss_1year", "mean_npp_1year", "mean_SST_1year", "pH_mean", "mean_DHW_1year")]
geo_var <- df_sel[, c("dist_to_coast", "latitude_start", "depth_fin", "dist_to_CT", "province", "depth_sampling")]
samp_var <- df[, c("sample_method", "volume")]
dist_jac_mo <- as.dist(dist_jac_mo)

varpart_tot <- varpart(dist_jac_mo, env_var, geo_var, samp_var)
varpart_tot

plot(varpart_tot, digits = 2, Xnames = c('env_var', 'geo_var', 'samp_var'), bg = c('navy', 'tomato', 'yellow'))

# get scores
station_scores <- scores(dbrda_fin)$sites
var_scores <- dbrda_fin[["CCA"]][["biplot"]][, 1:2] %>% as.data.frame()

# extract the percentage variability explained by axes
sumdbrda <- summary(dbrda_fin)
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

#---------------------------------------------------------------------------------------------------------------------------------
#### partial dbrda correcting for sampling and spatial ####

dbrda_part <- capscale(dist_jac_mo ~ dist_to_coast+mean_sss_1year+latitude_start+depth_fin+mean_npp_1year+mean_SST_1year+dist_to_CT+province+depth_sampling+pH_mean+mean_DHW_1year + Condition(MEM1+MEM2+MEM3+MEM4+MEM5+filter+sequencing_depth+volume+sample_method), df_sel)
RsquareAdj(dbrda_part)
anova(dbrda_part)
anova(dbrda_part, by = "term", permutations = 99)

vif(dbrda_part)
# variation partitioning
varpart_part <- varpart(dist_jac_mo, env_var, geo_var)
varpart_part

plot(varpart_part, digits = 2, Xnames = c('env_var', 'geo_var'), bg = c('navy', 'tomato'))


# get scores
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
