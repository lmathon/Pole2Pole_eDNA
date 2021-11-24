library(ggplot2)
library(tidyverse)
library(funk)

load("Rdata/motu_effectsize.rdata")
load("Rdata/mntd_effectsize.rdata")
load("Rdata/crypto_effectsize.rdata")
load("Rdata/large_effectsize.rdata")


effectsize <- rbind(MNTD_effectsize, motus_effectsize, crypto_effectsize, large_effectsize)

effectsize$taxa <- gsub("Richness - Cryptobenthics", "Richness - Crypto", effectsize$taxa)
effectsize$vargroup <- gsub("socio", "Socio-economy", effectsize$vargroup)
effectsize$vargroup <- gsub("environment", "Environment", effectsize$vargroup)
effectsize$vargroup <- gsub("geography", "Geography", effectsize$vargroup)
effectsize$vargroup <- gsub("sampling", "Samp.", effectsize$vargroup)
effectsize$Parameter <- gsub("bathy", "bathymetry", effectsize$Parameter)
effectsize$Parameter <- gsub("dist_to_CT", "distance to CT", effectsize$Parameter)
effectsize$Parameter <- gsub("distCoast", "distance to shore", effectsize$Parameter)
effectsize$Parameter <- gsub("depth_sampling", "depth of sampling", effectsize$Parameter)

effectsize$color <- NA
for (i in 1:nrow(effectsize)) {
  if((effectsize[i, "CI_low"]<0) | (effectsize[i,"CI_high"]<0)){
    effectsize[i,"color"] <- "black"
  }
  if ((effectsize[i, "CI_low"]>0) & (effectsize[i,"CI_high"]>0)){
    effectsize[i,"color"] <- "forestgreen"
  }
  if((effectsize[i, "CI_low"]<0) & (effectsize[i,"CI_high"]<0))
    effectsize[i,"color"] <- "red"
}


effectsize <- effectsize %>%
  mutate(across(vargroup, factor, levels=c("Environment","Socio-economy","Geography", "Samp.")))

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
        strip.text.x = element_text(size=10,face="bold"),
        strip.text.y = element_text(size=10,face="bold"),
        strip.placement = "outside")

ggsave("outputs/GLS/effect_size.png", width = 10, height = 6.5)
