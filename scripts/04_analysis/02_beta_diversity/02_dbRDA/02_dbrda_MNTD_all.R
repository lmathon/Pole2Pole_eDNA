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
library(ade4)

load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
load("Rdata/MNTD_pairwise_station.rdata")
load("Rdata/db_mem.rdata")

# transform data

data <- exp_var
df <- data %>%
  dplyr::select(-c("station")) # remove station
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

# check for colinearity and select variables
mctest::imcdiag(dbrda_full, method="VIF")


#### partial dbrda correcting for sampling and MEM ####
dbrda_part <- capscale(mntd ~ mean_DHW_1year+mean_SST_1year+mean_sss_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast +Condition(volume+MEM1), df_mem) 



RsquareAdj(dbrda_part)
anova(dbrda_part)
anova(dbrda_part, by = "term", permutations = 99)
anova(dbrda_part, by = "margin", permutations = 99)


# variation partitioning
#
env_var <- df_mem[,c("mean_DHW_1year","mean_SST_1year", "mean_sss_1year", "mean_npp_1year")]
geo_var <- df_mem[, c("bathy", "dist_to_CT", "depth_sampling", "distCoast")]
socio_var <- df_mem[,c("Voice_mean", "HDI2019", "Gravity", "MarineEcosystemDependency", "conflicts")]
mntd <- as.dist(mntd)


varpart_part <- varpart(mntd, env_var, geo_var, socio_var)

plot(varpart_part, digits = 2, Xnames = c('environment', 'geography', 'socio-economy'), bg = c('navy', 'tomato', 'yellow'))

# boxplot partition per variable type

partition <- data.frame(environment=0.116+0.017+0.054+0.202, 
                        geography=0.017+0.054+0.091, 
                        socioeconomy=0.124+0.054+0.202)
                        

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "socioeconomy", "geography"))

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
var_scores_diff75 <- var_scores_diff75[-c(6,9),]

dbrda_MNTD_prov <- ggplot(station_scores_met, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey", show.legend = F) +
  geom_vline(xintercept = 0, lty = 2, col = "grey", show.legend = F) +
  geom_point(cex = 2, show.legend = T, aes(col=province)) +
  #scale_color_gradient(low="blue", high="red")+
  scale_fill_brewer(palette="Paired", direction = 1, aesthetics = "col") +
  geom_segment(data= var_scores_diff75, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc")), show.legend = F) + # most differentiated variables
  geom_label_repel(data= var_scores_diff75, 
                   aes(x= CAP1, y=CAP2, 
                       fontface=3),
                   label = rownames(var_scores_diff75),
                   label.size = NA, 
                   size = 4,
                   fill = NA,
                   show.legend = F) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",             # position in top left corner
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_text(),
        legend.text = element_text(size=9),
        legend.key.height = unit(0.5, "cm"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
dbrda_MNTD_prov

ggsave("outputs/dbRDA/MNTD_all/dbrda_prov.png")
save(dbrda_MNTD_prov, file="Rdata/dbrda_MNTD_province.rdata")
