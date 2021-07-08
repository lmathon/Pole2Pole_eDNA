library(tidyverse)


#---------------------------------------------------------------------------------------------------------------------
# Environmental variables
#---------------------------------------------------------------------------------------------------------------------

env_var <- read.csv("metadata/eDNA-all-env-join.csv", sep=";")
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")

env_var <- left_join(env_var, meta[,c("code_spygen", "station")], by="code_spygen")

env_var <- env_var %>%
  distinct(station, .keep_all=T)

save(env_var, file="Rdata/environmental_variables.rdata")


#---------------------------------------------------------------------------------------------------------------------
# Geographic variables
#---------------------------------------------------------------------------------------------------------------------

geo_var <- meta[,c("station", "dist_to_coast", "dist_to_CT", "province", "depth_seafloor", "depth_sampling", "latitude_start")] # complete depth !

geo_var <- geo_var %>%
  distinct(station, .keep_all=T)

save(geo_var, file="Rdata/geographic_variables.rdata")


#---------------------------------------------------------------------------------------------------------------------
# Sampling variables
#---------------------------------------------------------------------------------------------------------------------

samp_var <- meta[,c("station", "sample_type", "sample_method", "filter", "sequencer", "sequencing_depth")] 

meta$volume <- as.numeric(meta$volume)
vol <- meta %>%
  select(station, volume) %>%
  group_by(station) %>% 
  summarise_all(funs(sum))

samp_var <- samp_var %>%
  distinct(station, .keep_all=T)

samp_var <- left_join(samp_var, vol)

save(samp_var, file="Rdata/sampling_variables.rdata")


#---------------------------------------------------------------------------------------------------------------------
# Human variables
#---------------------------------------------------------------------------------------------------------------------

