library(tidyverse)
library(ape)
library(vegan)

load("Rdata/03-filter-data.Rdata")

# prepare community matrix

df <- df_filtered %>%
  subset(!is.na(station))

amplicon <- unique(df$sequence)

stations <- unique(df$station)
st <- vector("list")
com <- data.frame(amplicon=amplicon)

for (i in 1:length(stations)) {
  st[[i]] <- df %>%
    subset(station == stations[i])%>%
    dplyr::select(sequence, count_reads) %>%
    group_by(sequence) %>% 
    summarise_all(funs(sum))
  colnames(st[[i]]) <- c("amplicon", stations[[i]])
  com <- left_join(com, st[[i]], by="amplicon")
  
}

rownames(com) <- com$amplicon
com <- com[,-1]
com[is.na(com)] <- 0
save(com, file="Rdata/station_presence_absence.rdata")
com <- as.data.frame(t(com))


#### calculate dissimilarity matrix with jaccard ####
jaccard_motu <- as.matrix(vegdist(com, method = "jaccard"))



save(jaccard_motu, file="Rdata/Jaccard_MOTU_dissimilarity.rdata")

# plot dissimilarity matrix

dissplot(jaccard_motu, method=NA, 
         upper_tri = TRUE, 
         lower_tri = FALSE, 
         reverse_columns=TRUE,
         main="Jaccard dissimilarity in MOTU composition",
         col=bluered(100))


#### Jaccard on chondri MOTUs ####
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
chondri_motu <- unique(df_chondri$definition)

com_chondri <- com[chondri_motu]

jaccard_chondri <- as.matrix(vegdist(com_chondri, method = "jaccard"))

save(jaccard_chondri, file="Rdata/Jaccard_chondri_dissimilarity.rdata")

#### MNTD on crypto MOTUs ####
cryptic_family <- c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae", "Kurtidae")
cryptic_order <- c("Kurtiformes", "Gobiiformes", "Blenniiformes", "Syngnathiformes")
df_crypto <- filter(df_filtered, order_name %in% cryptic_order | family_name_corrected %in% cryptic_family)

crypto_motu <- unique(df_crypto$definition)

com_crypto <- com[crypto_motu]

jaccard_crypto <- as.matrix(vegdist(com_crypto, method = "jaccard"))

save(jaccard_crypto, file="Rdata/Jaccard_crypto_dissimilarity.rdata")
