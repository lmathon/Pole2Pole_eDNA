library(tidyverse)
library(rcompanion)
library(ggplot2)

source("scripts/04_analysis/00_functions.R")

load("Rdata/environmental_variables.rdata")
load("Rdata/sampling_variables.rdata")
load("Rdata/geographic_variables.rdata")
#load("Rdata/socioeconomic_variables.rdata")
load("Rdata/MNTD_pairwise_station.rdata")

# Formate data
rownames(env_var) <- env_var$station
env_var <- env_var[rownames(mntd),]
env_var <- env_var[, colSums(env_var != 0) > 0]

rownames(geo_var) <- geo_var$station
geo_var <- geo_var[rownames(mntd),]

rownames(samp_var) <- samp_var$station
samp_var <- samp_var[rownames(mntd),]

#rownames(socio_var) <- socio_var$station
#socio_var <- socio_var[rownames(mntd),]

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
cor_geo_var <- mixed_assoc(geo_var[,-1]) ## everything ok, maybe remove province

cor_geo_sign <- cor_geo_var %>%
  filter(assoc >= 0.7 | assoc <= -0.7)

ggplot(cor_geo_var, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint=0)


# sampling
cor_samp_var <- mixed_assoc(samp_var[,-c(1,5)]) ## everythin ok, maybe remove province

cor_samp_sign <- cor_samp_var %>%
  filter(assoc >= 0.7 | assoc <= -0.7)

ggplot(cor_samp_var, aes(x,y,fill=assoc))+
  geom_tile()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(face="plain", size=10, angle=90, vjust = 0, hjust = 1))+
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint=0)

samp_var2 <- samp_var %>%
  select(station, sample_method, filter, sequencer, sequencing_depth, volume)

# socioeco








save(env_var2, file="Rdata/selected_environmental_variables.rdata")
save(samp_var2, file="Rdata/selected_sampling_variables.rdata")

#---------------------------------------------------------------------------------------------------
# Assemble all
#-----------------------------------------------------------------------------------------------------

exp_var <- cbind(env_var2, geo_var[,-1], samp_var2[,-1])



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
    sample_method == "transect_rond" ~ 9)) %>%
  mutate(province = case_when(
    province == "Arctic" ~ 1,
    province == "Cold_Temperate_Northwest_Pacific" ~ 2,
    province == "Lusitanian" ~ 3,
    province == "Mediterranean_Sea" ~ 4,
    province == "Northern_European_Seas" ~ 5,
    province == "Scotia_Sea" ~ 6,
    province == "Southeast_Polynesia" ~ 7,
    province == "Tropical_East_Pacific" ~ 8,
    province == "Tropical_Northwestern_Atlantic" ~ 9,
    province == "Tropical_Southwestern_Pacific" ~ 10,
    province == "Western_Coral_Triangle" ~ 11,
    province == "Western_Indian_Ocean" ~ 12
  ))



save(exp_var, file="Rdata/all_explanatory_variables.rdata")
save(exp_var_num, file="Rdata/all_explanatory_variables_numeric.rdata")
