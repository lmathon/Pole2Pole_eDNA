library(ggplot2)
library(tidyverse)
library(funk)

load("Rdata/motu_effectsize.rdata")
load("Rdata/FDq0_effectsize.rdata")
load("Rdata/crypto_effectsize.rdata")
load("Rdata/large_effectsize.rdata")


effectsize <- rbind(FDq0_effectsize, motus_effectsize, crypto_effectsize, large_effectsize)

effectsize$taxa <- gsub("Richness - Cryptobenthics", "Cryptobenthic a-diversity", effectsize$taxa)
effectsize$taxa <- gsub("Richness - Large fish", "Large fish a-diversity", effectsize$taxa)
effectsize$taxa <- gsub("Richness - all MOTUs", "All fish a-diversity", effectsize$taxa)
effectsize$taxa <- gsub("Functional a-diversity", "Sequence a-diversity", effectsize$taxa)
effectsize$vargroup <- gsub("socio", "Socio-economy", effectsize$vargroup)
effectsize$vargroup <- gsub("environment", "Environment", effectsize$vargroup)
effectsize$vargroup <- gsub("geography", "Geography", effectsize$vargroup)
effectsize$vargroup <- gsub("sampling", "Sampling", effectsize$vargroup)
effectsize$Parameter <- gsub("bathy", "bathymetry", effectsize$Parameter)
effectsize$Parameter <- gsub("dist_to_CT", "distance to CT", effectsize$Parameter)
effectsize$Parameter <- gsub("distCoast", "distance to shore", effectsize$Parameter)
effectsize$Parameter <- gsub("depth_sampling", "depth of sampling", effectsize$Parameter)
effectsize$Parameter <- gsub("sample_method2transect", "method_transect", effectsize$Parameter)

for (i in 1:nrow(effectsize)) {
  if(effectsize[i, "vargroup"]=="Sampling"){
    effectsize[i,"Std_Coefficient"] <- effectsize[i,"Std_Coefficient"]/max(effectsize$CI_high)
    effectsize[i,"CI_low"] <- effectsize[i,"CI_low"]/max(effectsize$CI_high)
    effectsize[i,"CI_high"] <- effectsize[i,"CI_high"]/max(effectsize$CI_high)
  }
}


effectsize$color <- NA
for (i in 1:nrow(effectsize)) {
  if((effectsize[i, "CI_low"]<0) | (effectsize[i,"CI_high"]<0)){
    effectsize[i,"color"] <- "black"
  }
  if ((effectsize[i, "CI_low"]>0) & (effectsize[i,"CI_high"]>0)){
    effectsize[i,"color"] <- "green4"
  }
  if((effectsize[i, "CI_low"]<0) & (effectsize[i,"CI_high"]<0))
    effectsize[i,"color"] <- "red"
}


effectsize <- effectsize %>%
  mutate(across(vargroup, factor, levels=c("Environment","Socio-economy","Geography", "Sampling")))

ggplot(data = effectsize, 
       aes(x = Parameter, y = Std_Coefficient)) +
  geom_hline(aes(yintercept = 0), colour = "black") + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                colour = "black", size = 0.5, width = 0) +
  geom_point(size = 2, col=effectsize$color) +
  coord_flip() + 
  theme_sleek(base_size = 24) + 
  facet_grid(vargroup ~ taxa, scales = "free_y", space = "free_y", switch = "y") + 
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1)) + 
  ylab("Standardized effect size") +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12),
        strip.text.y = element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12,face="bold"),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        strip.placement = "outside")

ggsave("outputs/GLS/FDq0_effect_size.png", width = 10, height = 6.5)
ggsave("outputs/Figures_papier/Fig2.pdf", width = 11, height = 6.5)
