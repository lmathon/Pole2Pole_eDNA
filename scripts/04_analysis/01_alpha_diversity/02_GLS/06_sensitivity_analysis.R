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

load("Rdata/MNTD_station.rdata")
load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data_init <- left_join(exp_var_num, rich_station, by="station")
data_init <- left_join(data_init, mntd_stations[,c("MNTD", "station")], by="station")

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

effectsize_fin <- vector("list", 10)

for (i in 1:10) {
  # select 80%
  
  data <- data_init %>%
    sample_frac(0.8)
  
  # GLS MOTUs
  
  gls.motus <- gls(MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Voice_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  
  motus_effectsize <- effectsize(gls.motus)
  motus_effectsize <- motus_effectsize[-1,]
  motus_effectsize$taxa <- "all MOTUs"
  motus_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")
  
  # GLS crypto
  
  gls.crypto <- gls(crypto_MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Voice_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  
  crypto_effectsize <- effectsize(gls.crypto)
  crypto_effectsize <- crypto_effectsize[-1,]
  crypto_effectsize$taxa <- "Cryptobenthics"
  crypto_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")
  
  # GLS large fish
  
  gls.largefish <- gls(largefish_MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Voice_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  large_effectsize <- effectsize(gls.largefish)
  large_effectsize <- large_effectsize[-1,]
  large_effectsize$taxa <- "Large fish"
  large_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")
  
  # GLS MNTD 
  gls.MNTD <- gls(MNTD ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Voice_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
  
  MNTD_effectsize <- effectsize(gls.MNTD)
  MNTD_effectsize <- MNTD_effectsize[-1,]
  MNTD_effectsize$taxa <- "ses.MNTD"
  MNTD_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")
  
   
  effectsize_fin[[i]] <- as.data.frame(rbind(MNTD_effectsize, motus_effectsize, crypto_effectsize, large_effectsize))
  
}


#### effect size x10 ####


effectsize_fin <- bind_rows(effectsize_fin)

effectsize_fin$vargroup <- gsub("socio", "Socio-economy", effectsize_fin$vargroup)
effectsize_fin$vargroup <- gsub("environment", "Environment", effectsize_fin$vargroup)
effectsize_fin$vargroup <- gsub("geography", "Geography", effectsize_fin$vargroup)
effectsize_fin$vargroup <- gsub("sampling", "Samp.", effectsize_fin$vargroup)
effectsize_fin$Parameter <- gsub("bathy", "bathymetry", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("dist_to_CT", "distance to CT", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("distCoast", "distance to shore", effectsize_fin$Parameter)
effectsize_fin$Parameter <- gsub("depth_sampling", "depth of sampling", effectsize_fin$Parameter)

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
  mutate(across(vargroup, factor, levels=c("Environment","Socio-economy","Geography", "Samp.")))

save(effectsize_fin, file = "Rdata/effect_size_sensitivity.rdata")
