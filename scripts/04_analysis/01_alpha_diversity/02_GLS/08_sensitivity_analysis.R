library(funk)
library(vegan)
library(ade4)
library(DHARMa)
library(car)
library(visreg)
library(ecospat)
library(modEvA)
library(psych)
library(MASS)
library(lmeInfo)
library(pscl)
library(spdep)
library(tidyverse)
library(nlme)
library(MuMIn)
library(rcompanion)
library(ggpubr)
library(ggplot2)
library(effectsize)

load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data_init <- left_join(exp_var_num, rich_station, by="station")
data_init <- data_init %>%
  dplyr::select(-c(station))

data_init$MOTUs <- log10(data_init$MOTUs +1)
data_init$crypto_MOTUs <- log10(data_init$crypto_MOTUs +1)
data_init$largefish_MOTUs <- log10(data_init$largefish_MOTUs +1)

# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- meta[, c("longitude_start", "latitude_start")]
data_init <- cbind(data_init, coor)

delta_rich_fin <- vector("list", 10)
effectsize_fin <- vector("list", 10)

for (i in 1:10) {
  # select 80%
  
  data <- data_init %>%
    sample_frac(0.8)
  
  # GLS MOTUs
  
  gls.motus <- gls(MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  fit.grav.motus <- visreg(gls.motus,"Gravity",scale="response")
  fit.MED.motus <- visreg(gls.motus,"MarineEcosystemDependency",scale="response")
  fit.grav_med.motus <- visreg2d(gls.motus, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(MOTUs richness +1)", xlab="log10(Gravity +1)")
  
  motus_effectsize <- effectsize(gls.motus)
  motus_effectsize <- motus_effectsize[-1,]
  motus_effectsize$taxa <- "all MOTUs"
  motus_effectsize$vargroup <- c("environment","environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")
  
  # GLS crypto
  
  gls.crypto <- gls(crypto_MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  fit.grav.crypto <- visreg(gls.crypto,"Gravity",scale="response")
  fit.MED.crypto <- visreg(gls.crypto,"MarineEcosystemDependency",scale="response")
  fit.grav_med.crypto <- visreg2d(gls.crypto, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(Crypto richness +1)", xlab="lg10(Gravity +1)")
  
  crypto_effectsize <- effectsize(gls.crypto)
  crypto_effectsize <- crypto_effectsize[-1,]
  crypto_effectsize$taxa <- "Cryptobenthics"
  crypto_effectsize$vargroup <- c("environment","environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")
  
  # GLS large fish
  
  gls.largefish <- gls(largefish_MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  fit.grav.large <- visreg(gls.largefish,"Gravity",scale="response")
  fit.MED.large <- visreg(gls.largefish,"MarineEcosystemDependency",scale="response")
  fit.grav_med.large <- visreg2d(gls.largefish, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(Large fish richness +1)", xlab="lg10(Gravity +1)")
  
  large_effectsize <- effectsize(gls.largefish)
  large_effectsize <- large_effectsize[-1,]
  large_effectsize$taxa <- "Large fish"
  large_effectsize$vargroup <- c("environment","environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")
  
  effectsize_fin[[i]] <- as.data.frame(rbind(motus_effectsize, crypto_effectsize, large_effectsize))
  
  # delta richness 
  delta_rich <- data.frame(taxa=character(), vargroup=character(), variable=character(), delta=numeric(), CI_low=numeric(), CI_high=numeric())
  
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
  
  delta_rich_fin[[i]] <- as.data.frame(delta_rich)
}

#### delta richness ####

delta_rich_fin <- bind_rows(delta_rich_fin)

delta_rich_fin <- delta_rich_fin %>%
  mutate(across(variable, factor, levels=c("Gravity + MED","Marine Ecosystem Dependency","Gravity")))

color <- c("red", "black", "forestgreen", "red", "dodgerblue4", "forestgreen", "red", "dodgerblue4", "forestgreen",
           "red", "black", "forestgreen", "red", "dodgerblue4", "forestgreen", "red", "dodgerblue4", "forestgreen",
           "red", "black", "black", "red", "dodgerblue4", "black", "red", "dodgerblue4", "forestgreen",
           "red", "black", "black", "red", "dodgerblue4", "forestgreen", "red", "dodgerblue4", "forestgreen",
           "red", "black", "black", "red", "dodgerblue4", "forestgreen", "red", "dodgerblue4", "forestgreen",
           "red", "black", "black", "red", "dodgerblue4", "forestgreen", "red", "dodgerblue4", "forestgreen",
           "red", "black", "black", "red", "dodgerblue4", "forestgreen", "red", "dodgerblue4", "forestgreen",
           "red", "black", "forestgreen", "red", "dodgerblue4", "forestgreen", "red", "dodgerblue4", "forestgreen",
           "red", "black", "forestgreen", "red", "dodgerblue4", "forestgreen", "red", "dodgerblue4", "forestgreen",
           "red", "black", "forestgreen", "red", "dodgerblue4", "forestgreen", "red", "dodgerblue4", "forestgreen")

ggplot(data = delta_rich_fin, 
       aes(x = variable, y = delta)) +
  geom_hline(aes(yintercept = 0), colour = "black") + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                colour = "black", size = 0.5, width = 0) +
  geom_point(size = 2, col=color) +
  coord_flip() + 
  ylim(-110, 150)+
  theme_sleek(base_size = 24) + 
  facet_grid(. ~ taxa, scales = "free_y", space = "free_y", switch = "y") + 
  #scale_y_continuous(breaks = c(-100,-50,0,50,100)) + 
  ylab(expression(paste(Delta,"% MOTU richness")))+
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12),
        strip.text.x = element_text(size=10,face="bold"),
        strip.text.y = element_text(size=10,face="bold"),
        strip.placement = "outside")

ggsave("outputs/GLS/delta_socio_sensitivity.png", width=6.5, height = 2)

save(delta_rich_fin, file = "Rdata/delta_richness_sensitivity.rdata")

#### effect size x5 ####


effectsize_fin <- bind_rows(effectsize_fin)

effectsize_fin$vargroup <- gsub("socio", "Socio-economy", effectsize_fin$vargroup)
effectsize_fin$vargroup <- gsub("environment", "Environment", effectsize_fin$vargroup)
effectsize_fin$vargroup <- gsub("geography", "Geography", effectsize_fin$vargroup)
effectsize_fin$vargroup <- gsub("sampling", "Samp.", effectsize_fin$vargroup)
effectsize_fin$Parameter <- gsub("bathy", "bathymetry", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("dist_to_CT", "distance to CT", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("distCoast", "distance to shore", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("depth_sampling", "depth of sampling", effectsize_fin$Parameter)

effectsize_fin <- effectsize_fin %>%
  mutate(across(vargroup, factor, levels=c("Environment","Socio-economy","Geography", "Samp.")))

ggplot(data = effectsize_fin, 
       aes(x = Parameter, y = Std_Coefficient)) +
  geom_hline(aes(yintercept = 0), colour = "black") + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                colour = "black", size = 0.5, width = 0) +
  geom_point(size = 2) +
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

ggsave("outputs/GLS/effect_size_sensitivity.png", width = 6, height = 6.5)
