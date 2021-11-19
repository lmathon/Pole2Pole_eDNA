library(tidyverse)

load("Rdata/03-filter-data.Rdata")

df <- df_filtered %>%
  subset(!is.na(station))

stations <- unique(df$station)

rich_station <- data.frame(station=stations, MOTUs=numeric(length(stations)), Family=numeric(length(stations)), Species=numeric(length(stations)), crypto_MOTUs=numeric(length(stations)), chondri_MOTUs=numeric(length(stations)), largefish_MOTUs=numeric(length(stations)), predator_MOTUs=numeric(length(stations)))

for (i in 1:length(stations)) {
  st <- as.data.frame(df %>%
    subset(station==stations[[i]]))
  rich_station[i,2] <- n_distinct(st$sequence)
  fam <- st %>% subset(!is.na(family_name_corrected))
  rich_station[i,3] <-  n_distinct(fam$family_name_corrected)
  sp <- st %>% subset(!is.na(species_name_corrected))
  rich_station[i,4] <- n_distinct(sp$species_name_corrected)
}




#### crypto richness ####
cryptic_family <- c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae", "Kurtidae")
cryptic_order <- c("Kurtiformes", "Gobiiformes", "Blenniiformes", "Syngnathiformes")
df_crypto <- filter(df_filtered, order_name %in% cryptic_order | family_name_corrected %in% cryptic_family)


for (i in 1:length(stations)) {
  st <- as.data.frame(df_crypto %>%
      subset(station==stations[[i]]))
  rich_station[i,5] <- n_distinct(st$sequence)
}

#### Chondri richness  ####
chondri_family <- c("Dasyatidae", "Potamotrygonidae", "Urotrygonidae", "Myliobatidae", "Gymnuridae", "Hexatrygonidae", 
                    "Plesiobatidae", "Urolophidae", "Anacanthobatidae", "Arhynchobatidae", "Rajidae", "Glaucostegidae",
                    "Pristidae", "Rhinidae", "Rhinobatidae", "Rhynchobatidae", "Zanobatidae", "Hypnidae", "Narcinidae",
                    "Narkidae", "Torpedinidae", "Platyrhinidae", "Carcharhinidae", "Hemigaleidae", "Leptochariidae",
                    "Proscylliidae", "Pseudotriakidae", "Scyliorhinidae", "Sphyrnidae", "Triakidae", "Alopiidae",
                    "Megachasmidae", "Mitsukurinidae", "Odontaspididae", "Pseudocarchariidae", "Brachaeluridae",
                    "Ginglymostomatidae", "Hemiscylliidae", "Orectolobidae", "Parascylliidae", "Rhincodontidae",
                    "Stegostomatidae", "Heterodontidae", "Echinorhinidae", "Chlamydoselachidae", "Hexanchidae", 
                    "Pristiophoridae", "Centrophoridae", "Dalatiidae", "Etmopteridae", "Oxynotidae", "Somniosidae", 
                    "Squalidae", "Squatinidae", "Callorhinchidae", "Chimaeridae", "Rhinochimaeridae")

chondri_order <- c("Myliobatiformes", "Rajiformes", "Rhinopristiformes", "Torpediniformes", "Chimaeriformes",
                   "Carcharhiniformes", "Lamniformes", "Orectolobiformes", "Heterodontiformes", "Echinorhiniformes",
                   "Hexanchiformes", "Pristiophoriformes", "Squaliformes", "Squatiniformes")

df_chondri <- filter(df_filtered, order_name %in% chondri_order | family_name_corrected %in% chondri_family)

for (i in 1:length(stations)) {
  st <- as.data.frame(df_chondri %>%
                        subset(station==stations[[i]]))
  rich_station[i,6] <- n_distinct(st$sequence)
}

#### Large fish richness ####

load("Rdata/large_families.rdata")
load("Rdata/large_orders.rdata")


df_large <- filter(df_filtered, order_name %in% large_orders | family_name_corrected %in% large_families)

for (i in 1:length(stations)) {
  st <- as.data.frame(df_large %>%
                        subset(station==stations[[i]]))
  rich_station[i,7] <- n_distinct(st$sequence)
}

#### Predators richness ####

pred_families <- c("Carangidae", "Carcharhinidae", "Ginglymostomatidae", "Heterodontidae", "Lutjanidae", "Serranidae", "Sphyraenidae", "Sphyrnidae")

df_predator <- filter(df_filtered, family_name_corrected %in% pred_families | class_name == "Chondrichthyes")

for (i in 1:length(stations)) {
  st <- as.data.frame(df_predator %>%
                        subset(station==stations[[i]]))
  rich_station[i,8] <- n_distinct(st$sequence)
}

save(rich_station, file="Rdata/richness_station.rdata")
