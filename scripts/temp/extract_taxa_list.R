library(tidyverse)
load("Rdata/03-filter-data.Rdata")

site <- c("Transect_Drake_6", "Transect_Drake_7", "Transect_Drake_10", "Transect_Drake_4", "Transect_Drake_3", "Transect_Drake_5", "Transect_Drake_9",
          "SVAL_21", "SVAL_10", "SVAL_14", "SVAL_24", "SVAL_22", "NOR_25")

species_list <- vector("list", 13)

for (i in 1:length(site)) {
  species_list[[i]] <- df_filtered %>%
    filter(station==site[i]) %>%
    distinct(sequence, .keep_all=T) %>%
    select(station, scientific_name_ncbi_corrected)
}


species_list <- bind_rows(species_list)

load("Rdata/MPD_station.rdata")
load("Rdata/MNTD_station.rdata")

species_list <- left_join(species_list, mpd_stations)
species_list <- left_join(species_list, mntd_stations)

write.csv(species_list, "c://Users/mathon/Desktop/species_list.csv", row.names = F)
