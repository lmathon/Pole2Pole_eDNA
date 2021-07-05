library(tidyverse)

load("Rdata/03-filter-data.Rdata")

df <- df_filtered %>%
  subset(!is.na(station))

stations <- unique(df$station)

rich_station <- data.frame(station=stations, MOTUs=numeric(length(stations)), Family=numeric(length(stations)))

for (i in 1:length(stations)) {
  st <- df %>%
    subset(station==stations[[i]])
  rich_station[i,2] <- n_distinct(st$sequence)
  rich_station[i,3] <- n_distinct(st$family_name_corrected)
}

save(rich_station, file="Rdata/richness_station.rdata")
