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
df <- data %>%
  select(-c("station", "province")) # remove station
df_mem <- cbind(df, dbmem)

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(data),]
identical(as.character(rownames(meta)), rownames(data))
coor <- meta[, c("longitude_start", "latitude_start")]
data <- cbind(data, coor)

#---------------------------------------------------------------------------------------------------------------------------
#### Full model ####

dbrda_full <- capscale(mntd ~ .,df_mem)

RsquareAdj(dbrda_full)
anova(dbrda_full)
anova(dbrda_full, by = "margin", permutations = 99)



# check for colinearity and select variables
mctest::imcdiag(dbrda_full, method="VIF")

# selection variables
dbrda0 <- capscale(mntd ~ 1, df_mem)
dbrdaG <- capscale(mntd ~ ., df_mem)
mem_sel <- ordiR2step(dbrda0, scope = formula(dbrdaG), direction="both")


#### partial dbrda correcting for sampling ####

dbrda_part <- capscale(mntd ~ mean_sss_1year+MEM4+mean_SST_1year+MEM3+NoViolence_mean+pH_mean+neartt+HDI2019+MEM5+bathy+conflicts+dist_to_CT+MEM2+MEM1+mean_npp_1year+MarineEcosystemDependency+depth_sampling+Corruption_mean+NGO+distCoast+mean_DHW_1year+Gravity +Condition(volume), df_mem) 

mctest::imcdiag(dbrda_part, method="VIF")

RsquareAdj(dbrda_part)
anova(dbrda_part)
anova(dbrda_part, by = "term", permutations = 99)
anova(dbrda_part, by = "margin", permutations = 99)


# variation partitioning
#
env_var <- df_mem[,c("mean_sss_1year", "mean_npp_1year", "mean_SST_1year", "mean_DHW_1year", "pH_mean")]
geo_var <- df_mem[, c("distCoast", "bathy", "dist_to_CT", "depth_sampling","MEM1", "MEM2", "MEM4", "MEM3", "MEM5")]
socio_var <- df_mem[,c("HDI2019", "neartt", "Gravity", "MarineEcosystemDependency", "Corruption_mean", "conflicts", "NoViolence_mean", "NGO")]
mntd <- as.dist(mntd)


varpart_part <- varpart(mntd, env_var, geo_var, socio_var)
varpart_part

plot(varpart_part, digits = 2, Xnames = c('environment', 'geography', 'socio-economy'), bg = c('navy', 'tomato', 'yellow'))


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
identical(as.character(rownames(data)), rownames(station_scores)) # verify that data in same order
station_scores_met <- cbind(station_scores, data)

grda_station <- ggplot(station_scores_met, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey", show.legend = F) +
  geom_vline(xintercept = 0, lty = 2, col = "grey", show.legend = F) +
  geom_encircle(aes(group = province, fill= province), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = TRUE) + # hull area 
  geom_point(col = "black", cex = 1, show.legend = F) +
  scale_fill_brewer(palette="Paired", direction = 1, aesthetics = "fill") +
  geom_segment(data= var_scores, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "grey",
               arrow=arrow(length=unit(0.01,"npc")), show.legend = F) + # all variables
  geom_segment(data= var_scores_diff75, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc")), show.legend = F) + # most differentiated variables
  geom_label_repel(data= var_scores_diff75, 
                   aes(x= CAP1, y=CAP2, 
                       fontface=3),
                   label = rownames(var_scores_diff75),
                   label.size = NA, 
                   size = 4,
                   fill = alpha(c("white"),0),
                   show.legend = F) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_station

#### partial dbrda correcting for sampling and MEM ####

dbrda_part <- capscale(mntd ~ mean_sss_1year+mean_SST_1year+NoViolence_mean+pH_mean+neartt+HDI2019+bathy+conflicts+dist_to_CT+mean_npp_1year+MarineEcosystemDependency+depth_sampling+Corruption_mean+NGO+distCoast+mean_DHW_1year+Gravity +Condition(volume+MEM4+MEM3+MEM5+MEM2+MEM1), df_mem) 


RsquareAdj(dbrda_part)
anova(dbrda_part)
anova(dbrda_part, by = "term", permutations = 99)
anova(dbrda_part, by = "margin", permutations = 99)


# variation partitioning
#
env_var <- df_mem[,c("mean_sss_1year", "mean_npp_1year", "mean_SST_1year", "mean_DHW_1year", "pH_mean")]
geo_var <- df_mem[, c("distCoast", "bathy", "dist_to_CT", "depth_sampling")]
socio_var <- df_mem[,c("HDI2019", "neartt", "Gravity", "MarineEcosystemDependency", "Corruption_mean", "conflicts", "NGO", "NoViolence_mean")]
mntd <- as.dist(mntd)


varpart_part <- varpart(mntd, env_var, geo_var, socio_var)
varpart_part

plot(varpart_part, digits = 2, Xnames = c('environment', 'geography', 'socio-economy'), bg = c('navy', 'tomato', 'yellow'))

# boxplot partition per variable type

partition <- data.frame(environment=0.112+0.147+0.007+0.011, 
                        geography=0.054+0.007+0.011+0.036, 
                        socioeconomy=0.093+0.147+0.011+0.036) 
                        

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("socioeconomy", "environment", "geography"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())

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
identical(as.character(rownames(data)), rownames(station_scores)) # verify that data in same order
station_scores_met <- cbind(station_scores, data)

grda_station <- ggplot(station_scores_met, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey", show.legend = F) +
  geom_vline(xintercept = 0, lty = 2, col = "grey", show.legend = F) +
  geom_encircle(aes(group = province, fill= province), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = TRUE) + # hull area 
  geom_point(col = "black", cex = 1, show.legend = F) +
  scale_fill_brewer(palette="Paired", direction = 1, aesthetics = "fill") +
  geom_segment(data= var_scores, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "grey",
               arrow=arrow(length=unit(0.01,"npc")), show.legend = F) + # all variables
  geom_segment(data= var_scores_diff75, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc")), show.legend = F) + # most differentiated variables
  geom_label_repel(data= var_scores_diff75, 
                   aes(x= CAP1, y=CAP2, 
                       fontface=3),
                   label = rownames(var_scores_diff75),
                   label.size = NA, 
                   size = 4,
                   fill = alpha(c("white"),0),
                   show.legend = F) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_station


# plot Axis1~latitude

ggplot(station_scores_met, aes(x=latitude_start, y=CAP1))+
  geom_point(aes(col=province))+
  scale_color_brewer(palette="Paired", direction = 1, aesthetics = "col")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
