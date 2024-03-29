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

'%ni%' <- Negate("%in%")

load("Rdata/richness_station.rdata")
load("Rdata/FD_Hill_alpha.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data <- left_join(exp_var, rich_station, by="station")
data <- left_join(data, FD_Hill[,c("FD_q0", "station")], by="station")

data <- data %>%
  dplyr::select(-c(station))

data$MOTUs <- log10(data$MOTUs +1)
data$crypto_MOTUs <- log10(data$crypto_MOTUs +1)
data$largefish_MOTUs <- log10(data$largefish_MOTUs +1)

# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- meta[, c("longitude_start", "latitude_start")]
data <- cbind(data, coor)
data$sample_method2 <- as.factor(data$sample_method2)

data <- data  %>%
  filter(province %ni% c("Scotia_Sea", "Arctic"))%>%
  dplyr::select(-province)
  
  
  # GLS MOTUs
  
  gls.motus <- gls(MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  motus_effectsize <- effectsize(gls.motus)
  motus_effectsize <- motus_effectsize[-1,]
  motus_effectsize$taxa <- "Richness - all MOTUs"
  motus_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","geography","geography","geography","geography","sampling","sampling")
  
  # GLS crypto
  
  gls.crypto <- gls(crypto_MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
   
  crypto_effectsize <- effectsize(gls.crypto)
  crypto_effectsize <- crypto_effectsize[-1,]
  crypto_effectsize$taxa <- "Richness - Crypto"
  crypto_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","geography","geography","geography","geography","sampling","sampling")
  
  # GLS large fish
  
  gls.largefish <- gls(largefish_MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  
  large_effectsize <- effectsize(gls.largefish)
  large_effectsize <- large_effectsize[-1,]
  large_effectsize$taxa <- "Richness - Large fish"
  large_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","geography","geography","geography","geography","sampling","sampling")
  
  
  #### GLS FD ####
  gls.FDq0 <- gls(FD_q0 ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+Gravity+MarineEcosystemDependency+dist_to_CT+bathy+depth_sampling+distCoast+volume+sample_method2, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  FDq0_effectsize <- effectsize(gls.FDq0)
  FDq0_effectsize <- FDq0_effectsize[-1,]
  FDq0_effectsize$taxa <- "sequence a-diversity"
  FDq0_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","geography","geography","geography","geography","sampling","sampling")
  
  
  effectsize_fin <- as.data.frame(rbind(FDq0_effectsize, motus_effectsize, crypto_effectsize, large_effectsize))
  

#### effect size ####

effectsize_fin$taxa <- gsub("Richness - Cryptobenthics", "Cryptobenthic a-diversity", effectsize_fin$taxa)
effectsize_fin$taxa <- gsub("Richness - Large fish", "Large fish a-diversity", effectsize_fin$taxa)
effectsize_fin$taxa <- gsub("Richness - all MOTUs", "All fish a-diversity", effectsize_fin$taxa)
effectsize_fin$taxa <- gsub("Functional a-diversity", "Sequence a-diversity", effectsize_fin$taxa)
effectsize_fin$vargroup <- gsub("socio", "Socio-economy", effectsize_fin$vargroup)
effectsize_fin$vargroup <- gsub("environment", "Environment", effectsize_fin$vargroup)
effectsize_fin$vargroup <- gsub("geography", "Geography", effectsize_fin$vargroup)
effectsize_fin$vargroup <- gsub("sampling", "Sampling", effectsize_fin$vargroup)
effectsize_fin$Parameter <- gsub("bathy", "bathymetry", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("dist_to_CT", "distance to CT", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("distCoast", "distance to shore", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("depth_sampling", "depth of sampling", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("MarineEcosystemDependency", "MarineEcosystemDependence", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("sample_method2transect", "method_transect", effectsize_fin$Parameter)

for (i in 1:nrow(effectsize_fin)) {
  if(effectsize_fin[i, "vargroup"]=="Sampling"){
    effectsize_fin[i,"Std_Coefficient"] <- effectsize_fin[i,"Std_Coefficient"]/max(effectsize_fin$CI_high)
    effectsize_fin[i,"CI_low"] <- effectsize_fin[i,"CI_low"]/max(effectsize_fin$CI_high)
    effectsize_fin[i,"CI_high"] <- effectsize_fin[i,"CI_high"]/max(effectsize_fin$CI_high)
  }
}

effectsize_fin$color <- NA
for (i in 1:nrow(effectsize_fin)) {
  if((effectsize_fin[i, "CI_low"]<0) | (effectsize_fin[i,"CI_high"]<0)){
    effectsize_fin[i,"color"] <- "black"
  }
  if ((effectsize_fin[i, "CI_low"]>0) & (effectsize_fin[i,"CI_high"]>0)){
    effectsize_fin[i,"color"] <- "forestgreen"
  }
  if((effectsize_fin[i, "CI_low"]<0) & (effectsize_fin[i,"CI_high"]<0))
    effectsize_fin[i,"color"] <- "red"
}

effectsize_fin <- effectsize_fin %>%
  mutate(across(vargroup, factor, levels=c("Environment","Socio-economy","Geography", "Sampling")))

ggplot(data = effectsize_fin, 
       aes(x = Parameter, y = Std_Coefficient)) +
  geom_hline(aes(yintercept = 0), colour = "black") + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                colour = "black", size = 0.5, width = 0) +
  geom_point(size = 2, col=effectsize_fin$color) +
  coord_flip() + 
  theme_sleek(base_size = 24) + 
  facet_grid(vargroup ~ taxa, scales = "free_y", space = "free_y", switch = "y") + 
  #scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1)) + 
  ylab("Standardized effect size") +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12),
        strip.text.x = element_text(size=10,face="bold"),
        strip.text.y = element_text(size=10,face="bold"),
        strip.placement = "outside")

ggsave("outputs/GLS/effect_size_without_poles.png", width = 10, height = 6.5)
