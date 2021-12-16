library(tidyverse)
library(funrar)

#import species lists

actino <- read.csv("data/liste_actino.csv")
chondri <- read.csv("data/liste_chondri.csv")
species <- rbind(actino, chondri) 

species <- as.factor(species$species_name_corrected)


# Import functional traits and compute functional distance matrix

traits=read.table("data/traits.txt")
head(traits)
traits$Species <- rownames(traits)
traits$Species <- gsub("_", " ", traits$Species)

traits_species <- traits %>%
  filter(Species %in% species)
traits_species <- traits_species[order(traits_species$Species), ]
rownames(traits_species) <- traits_species$Species
traits_species <- traits_species[,-8]

traits_species=na.omit(traits_species)

dist_trait = compute_dist_matrix(traits_species,metric = "euclidean", center = T, scale = T)


# Compute genetic distance matrix


library(ape)
library(picante)
library(seriation)
library(mFD)

load("Rdata/03-filter-data.Rdata")

seq_species <- df_filtered %>%
  filter(species_name_corrected %in% species)

seq <- seq_species %>%
  distinct(species_name_corrected,sequence)
seq <- seq %>%
  distinct(species_name_corrected, .keep_all=T)
seq <- seq[order(seq$species_name_corrected), ]
seq2 <- seq[,1]
rownames(seq2) <- seq$species_name_corrected

seq2 <- strsplit(seq2$sequence, split = character(0))
seq3 <- data.frame()

for (i in 1:length(seq2)) {
  seq3 <- qpcR:::rbind.na(seq3, seq2[[i]])
  
}

dist_gen <- as.matrix(dist.gene(seq3, method = "percentage"))
rownames(dist_gen) <- seq$species_name_corrected
colnames(dist_gen) <- seq$species_name_corrected
