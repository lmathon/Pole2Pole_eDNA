library(ggplot2)
library(tidyverse)
install.packages("remotes")
remotes::install_github("jpwrobinson/funk")
library(funk)

load("Rdata/motu_effectsize.rdata")
load("Rdata/crypto_effectsize.rdata")
load("Rdata/large_effectsize.rdata")

effectsize <- rbind(motus_effectsize, crypto_effectsize, large_effectsize)


effectsize$vargroup <- gsub("socio", "Socio-economy", effectsize$vargroup)
effectsize$vargroup <- gsub("environment", "Environment", effectsize$vargroup)
effectsize$vargroup <- gsub("geography", "Geography", effectsize$vargroup)
effectsize$vargroup <- gsub("sampling", "Samp.", effectsize$vargroup)
effectsize$Parameter <- gsub("bathy", "bathymetry", effectsize$Parameter)
effectsize$Parameter <- gsub("dist_to_CT", "distance to CT", effectsize$Parameter)
effectsize$Parameter <- gsub("distCoast", "distance to shore", effectsize$Parameter)
effectsize$Parameter <- gsub("depth_sampling", "depth of sampling", effectsize$Parameter)

color <- c("black", "black","red", "red", "black", "black", "red", "red", "red", "red", "red", "black", "black", "red", "red",
           "black", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "black", "dodgerblue4", "black", "dodgerblue4", "black", "dodgerblue4", "black", "dodgerblue4", "black", "dodgerblue4",
           "black", "black", "forestgreen", "forestgreen", "forestgreen", "black", "black", "forestgreen", "forestgreen", "black", "forestgreen", "forestgreen", "black", "forestgreen","forestgreen")

effectsize <- effectsize %>%
  mutate(across(vargroup, factor, levels=c("Environment","Socio-economy","Geography", "Samp.")))

ggplot(data = effectsize, 
       aes(x = Parameter, y = Std_Coefficient)) +
  geom_hline(aes(yintercept = 0), colour = "black") + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                colour = "black", size = 0.5, width = 0) +
  geom_point(size = 2, col=color) +
  coord_flip() + 
  theme_sleek(base_size = 24) + 
  facet_grid(vargroup ~ taxa, scales = "free_y", space = "free_y", switch = "y") + 
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1)) + 
  ylab("Standardized effect size") +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12),
        strip.text.x = element_text(size=12,face="bold"),
        strip.text.y = element_text(size=10,face="bold"),
        strip.placement = "outside")

ggsave("outputs/GLS/effect_size.png", width = 6, height = 6.5)
