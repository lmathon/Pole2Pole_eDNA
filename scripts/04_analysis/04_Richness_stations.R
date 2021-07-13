library(tidyverse)

load("Rdata/03-filter-data.Rdata")

df <- df_filtered %>%
  subset(!is.na(station))

stations <- unique(df$station)

rich_station <- data.frame(station=stations, MOTUs=numeric(length(stations)), Family=numeric(length(stations)), Species=numeric(length(stations)))

for (i in 1:length(stations)) {
  st <- as.data.frame(df %>%
    subset(station==stations[[i]]))
  rich_station[i,2] <- n_distinct(st$sequence)
  fam <- st %>% subset(!is.na(family_name_corrected))
  rich_station[i,3] <-  n_distinct(fam$family_name_corrected)
  sp <- st %>% subset(!is.na(species_name_corrected))
  rich_station[i,4] <- n_distinct(sp$species_name_corrected)
}

save(rich_station, file="Rdata/richness_station.rdata")
