library(tidyverse)
library(rcompanion)
library(ggplot2)

source("scripts/04_analysis/00_functions.R")

load("Rdata/all_explanatory_variables.rdata")
load("Rdata/MNTD_pairwise_station.rdata")

# Formate data
mntd <- as.matrix(mntd)
df <- exp_var %>% filter(station %in% rownames(mntd))
rownames(df) <- df$station
df <- df[rownames(mntd), ]
df <- df[,-62]
df <- df[, colSums(df != 0) > 0]


# calculate correlation between variables
cor_var <- mixed_assoc(df)

cor_sign <- cor_var %>%
  filter(assoc >= 0.5 | assoc <= -0.5)

plot <- ggplot(cor_sign, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient(low="blue", high="red")

