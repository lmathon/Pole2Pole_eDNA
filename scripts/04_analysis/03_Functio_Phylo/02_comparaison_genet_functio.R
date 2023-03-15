library(tidyverse)
library(funrar)
library(fishtree)
library(picante)
library(geiger)
library(rfishbase)
library(ape)
library(seriation)
library(ade4)
library(cowplot)
library(mFD)
library(ecodist)
library(coRanking)
library(funk)


#import species lists

actino <- read.csv("data/liste_actino.csv")
chondri <- read.csv("data/liste_chondri.csv")
species <- rbind(actino, chondri)
species <- as.character(species$species_name_corrected)
species <- species[order(species)]
species_val <- validate_names(species)

species_val[8] <- "Acanthostracion polygonius"
species_val[88] <- "Haemulon chrysargyreum"
species_val[255] <- "Encrasicholina heteroloba"
species_val[510] <- "Nectamia fusca"

species <- species[-c(53,76,242,283)]
species_val <- species_val[-c(53,76,242,283)]


# prepare genetic data

load("Rdata/03-filter-data.Rdata")

seq_species <- df_filtered %>%
  filter(species_name_corrected %in% species)

seq <- seq_species %>%
  distinct(species_name_corrected,sequence)
seq <- seq %>%
  distinct(species_name_corrected, .keep_all=T)
seq <- seq[order(seq$species_name_corrected), ]
identical(seq$species_name_corrected, species)
seq$species_name_corrected <- species_val
seq <- seq %>%
  distinct(species_name_corrected, .keep_all=T)
rownames(seq) <- seq$species_name_corrected

# Import functional traits 

load("Rdata/Data_trait_fb_3OCC_bathyend_TL_CML.rdata")
traits <- Data_trait_fb_3OCC_bathyend_TL_CML

traits$Species <- rownames(traits)
traits$Species <- gsub("_", " ", traits$Species)

traits_species <- traits %>%
  filter(Species %in% species_val)

traits_species <- traits_species[order(traits_species$Species), ]
rownames(traits_species) <- traits_species$Species
traits_species <- traits_species[,-c(1:6,15:19,21:23,25:27,29:50)]
traits_species <- traits_species[,-c(12:14,31:36)]
traits_species$DemersPelag <- as.factor(traits_species$DemersPelag)



# compute functional distance

dist_trait = compute_dist_matrix(traits_species, metric = "gower")

maxs <- apply(dist_trait, 2, max)
mins <- apply(dist_trait, 2, min)
dist_trait <- scale(dist_trait, center = mins, scale = maxs - mins)


# Compute genetic distance matrix

seq <- seq %>%
  filter(species_name_corrected %in% rownames(traits_species))
seq <- seq[order(seq$species_name_corrected), ]
rownames(seq) <- seq$species_name_corrected

identical(rownames(seq), rownames(traits_species))

seq2 <- strsplit(seq$sequence, split = character(0))
seq3 <- data.frame()

for (i in 1:length(seq2)) {
  seq3 <- qpcR:::rbind.na(seq3, seq2[[i]])
  
}

dist_gen <- as.matrix(dist.gene(seq3, method = "percentage"))
rownames(dist_gen) <- seq$species_name_corrected
colnames(dist_gen) <- seq$species_name_corrected

identical(rownames(dist_gen), rownames(dist_trait))



# Correlation between matrix

dist_gen <- as.matrix(dist_gen)
dist_trait <- as.matrix(dist_trait)

pairwise <- data.frame(gen=as.vector(dist_gen), trait=as.vector(dist_trait))
cor.test(pairwise$trait, pairwise$gen, method = "pearson")

ecodist::mantel(dist_gen ~ dist_trait)


co_rank <- coranking(dist_trait, dist_gen, input_Xi = "dist")
NX <- coRanking::R_NX(co_rank)
AUC <- coRanking::AUC_ln_K(NX)

# biplot

plot_functio_gen <- ggplot(pairwise, aes(gen, trait))+
  geom_point()+
  labs(y= "Functional distance", x="Genetic distance")+
  annotate(geom="text", x=0, y=1, label="Mantel=0.04", hjust=0, size=6, fontface = "bold")+
  theme_bw(base_size = 24)+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14))

save(plot_functio_gen, file="Rdata/plot_functio_gen.rdata")



# prepare community matrix

df <- df_filtered %>%
  subset(!is.na(station))

stations <- unique(df$station)
st <- vector("list")
com <- data.frame(species=species)

for (i in 1:length(stations)) {
  st[[i]] <- df %>%
    subset(station == stations[i])%>%
    dplyr::select(species_name_corrected, count_reads) %>%
    group_by(species_name_corrected) %>% 
    summarise_all(funs(sum))
  colnames(st[[i]]) <- c("species", stations[[i]])
  com <- left_join(com, st[[i]], by="species")
  
}

com$species <- species_val
com <- com %>%
  distinct(species, .keep_all=T)
rownames(com) <- com$species
com <- com[,-1]
com[is.na(com)] <- 0

com <- as.data.frame(t(com))
com[com > 1] <- 1

com <- com[, (colnames(com) %in% colnames(dist_gen))]
com <- com[rowSums(com)>0, ]


# Calculate alpha HILL for genet and functio per sample


com <- as.matrix(com)
Hill_gen <- alpha.fd.hill(com, dist_gen, q=0, tau = "mean")

dist_trait <- as.matrix(dist_trait)
Hill_trait <- alpha.fd.hill(com, dist_trait, q=0, tau = "mean")

Hill <- data.frame(station=rownames(com), genet=Hill_gen$asb_FD_Hill, trait=Hill_trait$asb_FD_Hill)
colnames(Hill) <- c("station", "genet", "trait")

cor.test(Hill$trait, Hill$genet, method = "pearson")


plot_alpha_trait_gen <- ggplot(Hill, aes(trait, genet))+
  geom_point()+
  labs(x= expression(paste("Functional ", alpha,"-diversity")), y=expression(paste("Sequence ", alpha,"-diversity")))+
  annotate(geom="text", x=1, y=28, label="Pearson cor=0.91 \n p<0.001", hjust=0, size=6, fontface = "bold")+
  theme_bw(base_size = 24)+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1, fill=NA))

save(plot_alpha_trait_gen, file="Rdata/plot_alpha_trait_gen.rdata")

# calculate beta HILL for genet and functio between samples

beta_hill_gen <- beta.fd.hill(com, dist_gen, q=0, tau = "mean", beta_type="Jaccard")
beta_hill_gen <- beta_hill_gen$beta_fd_q$q0

beta_hill_trait <- beta.fd.hill(com, dist_trait, q=0, tau = "mean", beta_type="Jaccard")
beta_hill_trait <- beta_hill_trait$beta_fd_q$q0


co_rank <- coranking(beta_hill_gen, beta_hill_trait, input_Xi = "dist")
NX <- coRanking::R_NX(co_rank)
AUC <- coRanking::AUC_ln_K(NX)


ecodist::mantel(as.dist(beta_hill_gen) ~ as.dist(beta_hill_trait))

beta_hill <- data.frame(gen=as.vector(beta_hill_gen), trait=as.vector(beta_hill_trait))
cor.test(beta_hill$trait, beta_hill$gen, method = "pearson")

plot_beta_trait_gen <- ggplot(beta_hill, aes(trait, gen))+
  geom_point()+
  labs(x= expression(paste("Functional ", beta,"-diversity")), y=expression(paste("Sequence ", beta,"-diversity")))+
  annotate(geom="text", x=0, y=1, label="Mantel=0.87", hjust=0, size=6, fontface = "bold")+
  theme_bw(base_size = 24)+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1, fill=NA))

save(plot_beta_trait_gen, file="Rdata/plot_beta_trait_gen.rdata")
