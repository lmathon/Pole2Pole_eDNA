library(tidyverse)
library(rcompanion)
library(ggplot2)
library(missMethods)

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
  select(station, mean_DHW_1year, mean_DHW_5year, mean_sss_1year, mean_SST_1year, mean_npp_1year, pH_mean)

cor_env_var2 <- mixed_assoc(env_var2[,-1])

ggplot(cor_env_var2, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint=0, breaks=c(-1, -0.5, 0, 0.5, 1))+
  theme_bw()


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
  select(station, dist_to_CT, depth_fin, depth_sampling, latitude_start, distCoast)
colnames(geo_var2) <- c("station", "dist_to_CT", "bathy", "depth_sampling", "latitude", "distCoast")

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
  select(station, sample_method, sequencer, volume)

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
  select(station, HDI2019, neartt, Gravity, NGO, MarineEcosystemDependency, Naturalresourcesrents)

# replace NA in Antarctica with median imputation
socio_var2 <- impute_median(socio_var2, type = "columnwise", ordered_low = FALSE)

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



# transform factor to numeric
exp_var_num <- exp_var %>%
  mutate(sample_method = case_when(
    sample_method == "control" ~ 0,
    sample_method == "bag_underwater" ~ 1,
    sample_method == "boat" ~ 2,
    sample_method == "niskin" ~ 3,
    sample_method == "transect_aller" ~ 4,
    sample_method == "transect_aller_retour" ~ 5,
    sample_method == "transect_benthique" ~ 6,
    sample_method == "transect_deep" ~ 7,
    sample_method == "transect_rectangle" ~ 8,
    sample_method == "transect_rond" ~ 9,
    sample_method == "bottle" ~ 10)) %>%
   mutate(sequencer = case_when(
    sequencer == "Miseq" ~ 1,
    sequencer == "Hiseq" ~ 2,
    sequencer == "NextSeq" ~ 3,
    sequencer == "IonTorrent" ~ 4))



save(exp_var, file="Rdata/all_explanatory_variables.rdata")
save(exp_var_num, file="Rdata/all_explanatory_variables_numeric.rdata")
