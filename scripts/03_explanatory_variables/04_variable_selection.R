library(tidyverse)
library(rcompanion)
library(ggplot2)
library(missMethods)
library(DMwR2)
library(fuzzySim)

source("scripts/04_analysis/00_functions.R")

load("Rdata/environmental_variables.rdata")
load("Rdata/sampling_variables.rdata")
load("Rdata/geographic_variables.rdata")
load("Rdata/socioeconomic_variables.rdata")
load("Rdata/FD_Hill_alpha.rdata")

# Formate data
rownames(env_var) <- env_var$station
env_var <- env_var[rownames(FD_Hill),]
env_var <- env_var[, colSums(env_var != 0) > 0]

rownames(geo_var) <- geo_var$station
geo_var <- geo_var[rownames(FD_Hill),]

rownames(samp_var) <- samp_var$station
samp_var <- samp_var[rownames(FD_Hill),]

rownames(socio_var) <- socio_var$station
socio_var <- socio_var[rownames(FD_Hill),]


# calculate correlation between variables and select variables

  # environment

env_var$mean_DHW_1year <- log10(env_var$mean_DHW_1year +1)
env_var$mean_npp_1year <- log10(env_var$mean_npp_1year +1)
env_var$mean_DHW_5year <- log10(env_var$mean_DHW_5year +1)

multicol(vars=env_var[,-1])

env_var <- env_var %>%
  dplyr::select(c("station", "mean_DHW_1year", "mean_DHW_5year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year", "pH_mean"))

multicol(vars=env_var[,-1])



# geographic

colnames(geo_var) <- c("station", "dist_to_CT", "province", "bathy", "depth_sampling", "distCoast")

# transform data to log(x+1)
geo_var$bathy <- gsub("-", "", geo_var$bathy)
geo_var$bathy <- as.numeric(geo_var$bathy)

geo_var$depth_sampling <- gsub("-", "", geo_var$depth_sampling)
geo_var$depth_sampling <- as.numeric(geo_var$depth_sampling)

geo_var$dist_to_CT <- log10(geo_var$dist_to_CT +1)
geo_var$distCoast <- log10(geo_var$distCoast +1)

multicol(vars=geo_var[,-1])


# sampling

samp_var <- samp_var %>%
  dplyr::select(station, sample_method, sequencer, volume)

samp_var$volume <- log10(samp_var$volume +1)


for(i in 1:nrow(samp_var)){
  if (samp_var[i, "sample_method"]=="transect_aller"){
    samp_var[i, "sample_method2"] <- "transect"
  }
  if (samp_var[i, "sample_method"]=="transect_rond"){
    samp_var[i, "sample_method2"] <- "transect"
  }
  if (samp_var[i, "sample_method"]=="transect_aller_retour"){
    samp_var[i, "sample_method2"] <- "transect"
  }
  if (samp_var[i, "sample_method"]=="transect_rectangle"){
    samp_var[i, "sample_method2"] <- "transect"
  }
  if (samp_var[i, "sample_method"]=="transect_benthique"){
    samp_var[i, "sample_method2"] <- "transect"
  }
  if (samp_var[i, "sample_method"]=="bag_underwater"){
    samp_var[i, "sample_method2"] <- "point"
  }
  if (samp_var[i, "sample_method"]=="bottle"){
    samp_var[i, "sample_method2"] <- "point"
  }
}

# socioeco

multicol(vars=socio_var[-c(1,10:12)])

socio_var <- socio_var %>%
  dplyr::select(station, HDI2019, neartt, Gravity, MarineEcosystemDependency)

socio_var$neartt <- log10(socio_var$neartt +1)
socio_var$Gravity <- log10(socio_var$Gravity +1)

multicol(vars=socio_var[,-1])



# save selected variables
save(env_var, file="Rdata/selected_environmental_variables.rdata")
save(geo_var, file="Rdata/selected_geographic_variables.rdata")
save(samp_var, file="Rdata/selected_sampling_variables.rdata")
save(socio_var, file="Rdata/selected_socioeconomic_variables.rdata")

#---------------------------------------------------------------------------------------------------
# Assemble all
#-----------------------------------------------------------------------------------------------------

var1 <- left_join(env_var, geo_var, by="station")
var2 <- left_join(var1, samp_var, by="station")
exp_var <- left_join(var2, socio_var, by="station")

multicol(vars=exp_var[,-c(1,9,13,14,16)])


# separate numeric and categorical
exp_var_num <- exp_var %>%
  dplyr::select( -c(sample_method, sample_method2, sequencer, province))

exp_var_cat <- exp_var %>%
  dplyr::select(sample_method, sample_method2, sequencer, province)

# save all variables
save(exp_var, file="Rdata/all_explanatory_variables.rdata")
# save numeric variables
save(exp_var_num, file="Rdata/all_explanatory_variables_numeric.rdata")
# save categorical variables
save(exp_var_cat, file="Rdata/all_explanatory_variables_categorical.rdata")
