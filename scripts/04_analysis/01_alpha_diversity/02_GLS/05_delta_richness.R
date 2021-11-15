library(ggplot2)
library(grid)
library(cowplot)
library(tidyverse)

delta_rich <- data.frame(taxa=character(), vargroup=character(), variable=character(), delta=numeric(), CI_low=numeric(), CI_high=numeric())

#### gravity ####

load("Rdata/fit.grav.motus.rdata")
load("Rdata/fit.grav.crypto.rdata")
load("Rdata/fit.grav.large.rdata")

delta_rich[c(1:3),"taxa"] <- c("All MOTUs", "Cryptobenthics", "Large fish")
delta_rich[c(1:3),"vargroup"] <- "Socio-economy"
delta_rich[c(1:3),"variable"] <- "Gravity"

fit_motu <- fit.grav.motus$fit

delta_rich[1,"delta"] <- (10^(fit_motu[101,"visregFit"]) - 10^(fit_motu[1, "visregFit"]))*100 / 10^(fit_motu[1, "visregFit"])
delta_rich[1,"CI_high"] <- (10^(fit_motu[101,"visregUpr"]) - 10^(fit_motu[1, "visregLwr"]))*100 / 10^(fit_motu[1, "visregLwr"])
delta_rich[1,"CI_low"] <- (10^(fit_motu[101,"visregLwr"]) - 10^(fit_motu[1, "visregUpr"]))*100 / 10^(fit_motu[1, "visregUpr"])

fit_crypto <- fit.grav.crypto$fit

delta_rich[2,"delta"] <- (10^(fit_crypto[101,"visregFit"]) - 10^(fit_crypto[1, "visregFit"]))*100 / 10^(fit_crypto[1, "visregFit"])
delta_rich[2,"CI_high"] <- (10^(fit_crypto[101,"visregUpr"]) - 10^(fit_crypto[1, "visregLwr"]))*100 / 10^(fit_crypto[1, "visregLwr"])
delta_rich[2,"CI_low"] <- (10^(fit_crypto[101,"visregLwr"]) - 10^(fit_crypto[1, "visregUpr"]))*100 / 10^(fit_crypto[1, "visregUpr"])

fit_large <- fit.grav.large$fit

delta_rich[3,"delta"] <- (10^(fit_large[101,"visregFit"]) - 10^(fit_large[1, "visregFit"]))*100 / 10^(fit_large[1, "visregFit"])
delta_rich[3,"CI_high"] <- (10^(fit_large[101,"visregUpr"]) - 10^(fit_large[1, "visregLwr"]))*100 / 10^(fit_large[1, "visregLwr"])
delta_rich[3,"CI_low"] <- (10^(fit_large[101,"visregLwr"]) - 10^(fit_large[1, "visregUpr"]))*100 / 10^(fit_large[1, "visregUpr"])

#### Marine Ecosystem Dependency ####

load("Rdata/fit.med.motus.rdata")
load("Rdata/fit.med.crypto.rdata")
load("Rdata/fit.med.large.rdata")

delta_rich[c(4:6),"taxa"] <- c("All MOTUs", "Cryptobenthics", "Large fish")
delta_rich[c(4:6),"vargroup"] <- "Socio-economy"
delta_rich[c(4:6),"variable"] <- "Marine Ecosystem Dependency"

fit_motu <- fit.MED.motus$fit

delta_rich[4,"delta"] <- (10^(fit_motu[101,"visregFit"]) - 10^(fit_motu[1, "visregFit"]))*100 / 10^(fit_motu[1, "visregFit"])
delta_rich[4,"CI_high"] <- (10^(fit_motu[101,"visregUpr"]) - 10^(fit_motu[1, "visregLwr"]))*100 / 10^(fit_motu[1, "visregLwr"])
delta_rich[4,"CI_low"] <- (10^(fit_motu[101,"visregLwr"]) - 10^(fit_motu[1, "visregUpr"]))*100 / 10^(fit_motu[1, "visregUpr"])

fit_crypto <- fit.MED.crypto$fit

delta_rich[5,"delta"] <- (10^(fit_crypto[101,"visregFit"]) - 10^(fit_crypto[1, "visregFit"]))*100 / 10^(fit_crypto[1, "visregFit"])
delta_rich[5,"CI_high"] <- (10^(fit_crypto[101,"visregUpr"]) - 10^(fit_crypto[1, "visregLwr"]))*100 / 10^(fit_crypto[1, "visregLwr"])
delta_rich[5,"CI_low"] <- (10^(fit_crypto[101,"visregLwr"]) - 10^(fit_crypto[1, "visregUpr"]))*100 / 10^(fit_crypto[1, "visregUpr"])

fit_large <- fit.MED.large$fit

delta_rich[6,"delta"] <- (10^(fit_large[101,"visregFit"]) - 10^(fit_large[1, "visregFit"]))*100 /10^(fit_large[1, "visregFit"]) 
delta_rich[6,"CI_high"] <- (10^(fit_large[101,"visregUpr"]) - 10^(fit_large[1, "visregLwr"]))*100 / 10^(fit_large[1, "visregLwr"])
delta_rich[6,"CI_low"] <- (10^(fit_large[101,"visregLwr"]) - 10^(fit_large[1, "visregUpr"]))*100 / 10^(fit_large[1, "visregUpr"])


#### Gravity + Marine Ecosystem Dependency ####

load("Rdata/fit.grav_med.motus.rdata")
load("Rdata/fit.grav_med.crypto.rdata")
load("Rdata/fit.grav_med.large.rdata")

delta_rich[c(7:9),"taxa"] <- c("All MOTUs", "Cryptobenthics", "Large fish")
delta_rich[c(7:9),"vargroup"] <- "Socio-economy"
delta_rich[c(7:9),"variable"] <- "Gravity + MED"

fit_motu <- as.data.frame(fit.grav_med.motus$z)

delta_rich[7,"delta"] <- (10^(fit_motu[99,99]) - 10^(fit_motu[1,1]))*100 / 10^(fit_motu[1,1])
delta_rich[7,"CI_high"] <- NA
delta_rich[7,"CI_low"] <- NA

fit_crypto <- as.data.frame(fit.grav_med.crypto$z)

delta_rich[8,"delta"] <- (10^(fit_crypto[99,99]) - 10^(fit_crypto[1,1]))*100 / 10^(fit_crypto[1,1])
delta_rich[8,"CI_high"] <- NA
delta_rich[8,"CI_low"] <- NA

fit_large <- as.data.frame(fit.grav_med.large$z)

delta_rich[9,"delta"] <- (10^(fit_large[99,99]) - 10^(fit_large[1,1]))*100 / 10^(fit_large[1,1])
delta_rich[9,"CI_high"] <- NA
delta_rich[9,"CI_low"] <- NA

save(delta_rich, file="Rdata/delta_richness.rdata")

#### plot ####
delta_rich <- delta_rich %>%
  mutate(across(variable, factor, levels=c("Gravity + MED","Marine Ecosystem Dependency","Gravity")))

color <- c("red", "black", "black", "red", "red", "red", "red", "red", "red")

ggplot(data = delta_rich, 
       aes(x = variable, y = delta)) +
  geom_hline(aes(yintercept = 0), colour = "black") + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                colour = "black", size = 0.5, width = 0) +
  geom_point(size = 2, col=color) +
  coord_flip() + 
  theme_sleek(base_size = 24) + 
  facet_grid(. ~ taxa, scales = "free_y", space = "free_y", switch = "y") + 
  scale_y_continuous(breaks = c(-100,-50,0,50,100)) + 
  ylab(expression(paste(Delta,"% MOTU richness")))+
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12),
        strip.text.x = element_text(size=10,face="bold"),
        strip.text.y = element_text(size=10,face="bold"),
        strip.placement = "outside")

ggsave("outputs/GLS/delta_socio.png", width = 6.5, height = 2)

