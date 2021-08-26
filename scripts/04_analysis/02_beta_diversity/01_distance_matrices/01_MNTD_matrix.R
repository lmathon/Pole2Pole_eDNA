library(tidyverse)
library(ape)
library(picante)
library(seriation)

load("Rdata/03-filter-data.Rdata")

# prepare pairwise sequence dissimilarity matrix
seq <- df_filtered %>%
  distinct(sequence, definition)
seq2 <- seq[,1]
rownames(seq2) <- seq$definition
seq2 <- strsplit(seq2$sequence, split = character(0))
seq3 <- data.frame()

for (i in 1:length(seq2)) {
  seq3 <- qpcR:::rbind.na(seq3, seq2[[i]])
  
}

dist_gen <- as.matrix(dist.gene(seq3, method = "percentage"))
rownames(dist_gen) <- seq$definition
colnames(dist_gen) <- seq$definition


# prepare community matrix

df <- df_filtered %>%
  subset(!is.na(station))

stations <- unique(df$station)
st <- vector("list")
com <- data.frame(amplicon=seq$definition)

for (i in 1:length(stations)) {
  st[[i]] <- df %>%
    subset(station == stations[i])%>%
    dplyr::select(definition, count_reads) %>%
    group_by(definition) %>% 
    summarise_all(funs(sum))
  colnames(st[[i]]) <- c("amplicon", stations[[i]])
  com <- left_join(com, st[[i]], by="amplicon")
  
}

rownames(com) <- com$amplicon
com <- com[,-1]
com <- as.data.frame(t(com))
com[is.na(com)] <- 0

#### calculate MNTD between pairs of stations ####

mntd <- comdistnt(com, dist_gen, abundance.weighted = FALSE, exclude.conspecifics = FALSE)
mntd <- as.matrix(mntd)


# save Rdata

save(mntd, file="Rdata/MNTD_pairwise_station.rdata")

# plot distance matrix

dissplot(mntd, method=NA, 
         upper_tri = TRUE, 
         lower_tri = FALSE, 
         reverse_columns=TRUE,
         main="Mean Nearest Taxon Distance between stations",
         col=bluered(100))


#### MNTD for chondri MOTUs ####
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

# compute MNTD
mntd_chondri <- comdistnt(com_chondri, dist_gen, abundance.weighted = FALSE, exclude.conspecifics = FALSE)
mntd_chondri <- as.matrix(mntd_chondri)


# save Rdata
save(mntd_chondri, file="Rdata/MNTD_chondri_pairwise_station.rdata")

#### MNTD on crypto MOTUs ####
cryptic_family <- c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae", "Kurtidae")
cryptic_order <- c("Kurtiformes", "Gobiiformes", "Blenniiformes", "Syngnathiformes")
df_crypto <- filter(df_filtered, order_name %in% cryptic_order | family_name_corrected %in% cryptic_family)

crypto_motu <- unique(df_crypto$definition)

com_crypto <- com[crypto_motu]

# compute MNTD
mntd_crypto <- comdistnt(com_crypto, dist_gen, abundance.weighted = FALSE, exclude.conspecifics = FALSE)
mntd_crypto <- as.matrix(mntd_crypto)

# save Rdata
save(mntd_crypto, file="Rdata/MNTD_crypto_pairwise_station.rdata")


