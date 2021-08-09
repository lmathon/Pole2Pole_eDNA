library(tidyverse)
library(rcompanion)
library(ggplot2)
library(missMethods)
library(DMwR2)

source("scripts/04_analysis/00_functions.R")

load("Rdata/environmental_variables.rdata")
load("Rdata/sampling_variables.rdata")
load("Rdata/geographic_variables.rdata")
load("Rdata/socioeconomic_variables.rdata")
load("Rdata/MNTD_pairwise_station.rdata")

# Formate data
rownames(env_var) <- env_var$station
env_var <- env_var[rownames(mntd),]
env_var <- env_var[, colSums(env_var != 0) > 0]

rownames(geo_var) <- geo_var$station
geo_var <- geo_var[rownames(mntd),]

rownames(samp_var) <- samp_var$station
samp_var <- samp_var[rownames(mntd),]

rownames(socio_var) <- socio_var$station
socio_var <- socio_var[rownames(mntd),]

# calculate correlation between variables and select variables

  # environment
cor_env_var <- mixed_assoc(env_var[,-1])

cor_env_sign <- cor_env_var %>%
  filter(assoc >= 0.7 | assoc <= -0.7)

ggplot(cor_env_var, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint=0)


env_var2 <- env_var %>% # selection of non colinear variables
  dplyr::select(station, mean_DHW_1year, mean_DHW_5year, mean_sss_1year, mean_SST_1year, mean_npp_1year, pH_mean)

cor_env_var2 <- mixed_assoc(env_var2[,-1])

ggplot(cor_env_var2, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint=0, breaks=c(-1, -0.5, 0, 0.5, 1))+
  theme_bw()

    # transform data to log(x+1)
hist(env_var2$mean_sss_1year, col = "grey")
hist(env_var2$pH_mean, col = "grey")
hist(env_var2$mean_SST_1year, col = "grey")
hist(log1p(env_var2$mean_npp_1year), col = "grey")
hist(log1p(env_var2$mean_DHW_1year), col = "grey")
hist(log1p(env_var2$mean_DHW_5year), col = "grey")


env_var2$mean_DHW_1year <- log1p(env_var2$mean_DHW_1year)
env_var2$mean_npp_1year <- log1p(env_var2$mean_npp_1year)
env_var2$mean_DHW_5year <- log1p(env_var2$mean_DHW_5year)


# geographic
cor_geo_var <- mixed_assoc(geo_var[,-1]) 

cor_geo_sign <- cor_geo_var %>%
  filter(assoc >= 0.7 | assoc <= -0.7)

ggplot(cor_geo_var, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint=0)

geo_var2 <- geo_var %>%
  dplyr::select(station, province, dist_to_CT, depth_fin, depth_sampling, distCoast)
colnames(geo_var2) <- c("station", "province", "dist_to_CT", "bathy", "depth_sampling", "distCoast")

    # transform data to log(x+1)
geo_var2$bathy <- gsub("-", "", geo_var2$bathy)
geo_var2$bathy <- as.numeric(geo_var2$bathy)

geo_var2$depth_sampling <- gsub("-", "", geo_var2$depth_sampling)
geo_var2$depth_sampling <- as.numeric(geo_var2$depth_sampling)



hist(log1p(geo_var2$dist_to_CT), col="grey")
hist(geo_var2$bathy, col="grey")
hist(geo_var2$depth_sampling, col="grey")
hist(log1p(geo_var2$distCoast), col="grey")

geo_var2$dist_to_CT <- log1p(geo_var2$dist_to_CT)
geo_var2$distCoast <- log1p(geo_var2$distCoast)

# sampling
cor_samp_var <- mixed_assoc(samp_var[,-1]) 

cor_samp_sign <- cor_samp_var %>%
  filter(assoc >= 0.7 | assoc <= -0.7)

ggplot(cor_samp_sign, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint=0)

samp_var2 <- samp_var %>%
  dplyr::select(station, sample_method, sequencer, volume)

    # transform data to log(x+1)

hist(log1p(samp_var2$volume), col="grey")

samp_var2$volume <- log1p(samp_var2$volume)

# socioeco

cor_socio_var <- mixed_assoc(socio_var[,-1])


cor_socio_sign <- cor_socio_var %>%
  filter(assoc >= 0.7 | assoc <= -0.7)

ggplot(cor_socio_sign, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint=0)

socio_var2 <- socio_var %>%
  dplyr::select(station, NoViolence_mean, Corruption_mean, HDI2019, neartt, Gravity, NGO, MarineEcosystemDependency, conflicts)


hist(socio_var2$HDI2019, col="grey")
hist(socio_var2$NoViolence_mean, col="grey")
hist(socio_var2$Corruption_mean, col="grey")
hist(socio_var2$conflicts, col="grey")
hist(log1p(socio_var2$neartt), col="grey")
hist(log1p(socio_var2$Gravity), col="grey")
hist(log1p(socio_var2$NGO), col="grey")
hist(socio_var2$MarineEcosystemDependency, col="grey")


socio_var2$neartt <- log1p(socio_var2$neartt)
socio_var2$Gravity <- log1p(socio_var2$Gravity)
socio_var2$NGO <- log1p(socio_var2$NGO)


# save selected variables
save(env_var2, file="Rdata/selected_environmental_variables.rdata")
save(geo_var2, file="Rdata/selected_geographic_variables.rdata")
save(samp_var2, file="Rdata/selected_sampling_variables.rdata")
save(socio_var2, file="Rdata/selected_socioeconomic_variables.rdata")

#---------------------------------------------------------------------------------------------------
# Assemble all
#-----------------------------------------------------------------------------------------------------

exp_var <- cbind(env_var2, socio_var2[,-1], geo_var2[,-1], samp_var2[,-1])

cor_exp_var <- mixed_assoc(exp_var[,-1])

cor_var_sign <- cor_exp_var %>%
  filter(assoc >= 0.7 | assoc <= -0.7)

ggplot(cor_var_sign, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint=0)



# separate numeric and categorical
exp_var_num <- exp_var %>%
  dplyr::select( -c(sample_method, sequencer, province))

exp_var_cat <- exp_var %>%
  dplyr::select(sample_method, sequencer, province)

# save all variables
save(exp_var, file="Rdata/all_explanatory_variables.rdata")
# save numeric variables
save(exp_var_num, file="Rdata/all_explanatory_variables_numeric.rdata")
# save categorical variables
save(exp_var_cat, file="Rdata/all_explanatory_variables_categorical.rdata")
