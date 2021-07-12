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
# Human variables
#---------------------------------------------------------------------------------------------------------------------




#---------------------------------------------------------------------------------------------------
# Assemble all
#-----------------------------------------------------------------------------------------------------

exp_var <- left_join(env_var, geo_var, by="station")
exp_var <- left_join(exp_var, samp_var, by="station")


# transform factor to numeric
exp_var_num <- exp_var %>%
  mutate(sample_type = case_when(
    sample_type == "control" ~ 0,
    sample_type == "transect" ~ 1,
    sample_type == "point" ~ 2)) %>%
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
