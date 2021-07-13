library(tidyverse)
'%ni%' <- Negate("%in%")

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=",")



#---------------------------------------------------------------------------------------------------------------------
# Environmental variables
#---------------------------------------------------------------------------------------------------------------------

env_var <- read.csv("metadata/eDNA-all-env-join.csv", sep=";")


env_var <- left_join(meta[,c("code_spygen", "station")], env_var, by="code_spygen")

env_var <- env_var %>%
  distinct(station, .keep_all=T) %>%
  filter(station != "")

env_var <- env_var[, -c(1,3)]

save(env_var, file="Rdata/environmental_variables.rdata")


#---------------------------------------------------------------------------------------------------------------------
# Geographic variables
#---------------------------------------------------------------------------------------------------------------------

geo_var <- meta[,c("station", "dist_to_coast", "dist_to_CT", "province", "depth_fin", "depth_sampling", "latitude_start")] # complete depth !

geo_var <- geo_var %>%
  distinct(station, .keep_all=T) %>%
  filter(station != "")

save(geo_var, file="Rdata/geographic_variables.rdata")


#---------------------------------------------------------------------------------------------------------------------
# Sampling variables
#---------------------------------------------------------------------------------------------------------------------

samp_var <- meta[,c("station", "sample_type", "sample_method", "filter", "sequencer", "sequencing_depth")] 

meta$volume <- as.numeric(meta$volume)
vol <- meta %>%
  dplyr::select(station, volume) %>%
  group_by(station) %>% 
  summarise_all(funs(sum))

samp_var <- samp_var %>%
  distinct(station, .keep_all=T) %>%
  filter(station != "")

samp_var <- left_join(samp_var, vol)

save(samp_var, file="Rdata/sampling_variables.rdata")


#---------------------------------------------------------------------------------------------------------------------
# Socioeconomic variables variables
#---------------------------------------------------------------------------------------------------------------------




