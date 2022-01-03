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


# Import phylogenetic data
phylo <- read.csv("data/rabosky_with_codes.csv")

phylo_species <- phylo %>%
  filter(sci_name %in% species_val)



# compute phylogenetic distance

rownames(phylo_species) <- phylo_species$sci_name
phylo_species <- phylo_species[order(phylo_species$sci_name),]

phy <- fishtree_phylogeny(species = phylo_species$sci_name)

dist_phylo=cophenetic.phylo(phy)
dist_phylo <- dist_phylo[order(rownames(dist_phylo)),]
dist_phylo <- dist_phylo[,order(colnames(dist_phylo))]

rownames(dist_phylo) <- gsub("_", " ", rownames(dist_phylo))
colnames(dist_phylo) <- gsub("_", " ", colnames(dist_phylo))

maxs <- apply(dist_phylo, 2, max)
mins <- apply(dist_phylo, 2, min)
dist_phylo <- scale(dist_phylo, center = mins, scale = maxs - mins)

# Compute genetic distance matrix

seq <- seq %>%
  filter(species_name_corrected %in% rownames(phylo_species))
seq <- seq[order(seq$species_name_corrected), ]
rownames(seq) <- seq$species_name_corrected

identical(rownames(seq), rownames(dist_phylo))

seq2 <- strsplit(seq$sequence, split = character(0))
seq3 <- data.frame()

for (i in 1:length(seq2)) {
  seq3 <- qpcR:::rbind.na(seq3, seq2[[i]])
  
}

dist_gen <- as.matrix(dist.gene(seq3, method = "percentage"))
rownames(dist_gen) <- seq$species_name_corrected
colnames(dist_gen) <- seq$species_name_corrected

identical(rownames(dist_gen), rownames(dist_phylo))

# plot 
dist_gen <- as.dist(dist_gen)
dist_phylo <- as.dist(dist_phylo)

plot_phylo <- dissplot(dist_phylo, method=NA, 
         upper_tri = TRUE, 
         lower_tri = FALSE, 
         reverse_columns=TRUE,
         main="Distance phylogenetic",
         col=bluered(100))

plot_seq <- dissplot(dist_gen, method=NA, 
                       upper_tri = TRUE, 
                       lower_tri = FALSE, 
                       reverse_columns=TRUE,
                       main="Distance sequence",
                       col=bluered(100))

# correlation between matrix
dist_gen <- as.matrix(dist_gen)
dist_phylo <- as.matrix(dist_phylo)

ecodist::mantel(dist_phylo ~ dist_gen)


co_rank <- coranking(dist_phylo, dist_gen, input_Xi = "dist")
NX <- coRanking::R_NX(co_rank)
AUC <- coRanking::AUC_ln_K(NX)

# biplot

plot(dist_phylo ~ dist_gen, xlab="Genetic distance", ylab="Phylogenetic distance", col="black")


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


# Calculate alpha HILL for genet and phylo per sample

dist_gen <- as.matrix(dist_gen)
com <- as.matrix(com)
Hill_gen <- alpha.fd.hill(com, dist_gen, q=2, tau = "mean")

dist_phylo <- as.matrix(dist_phylo)
Hill_phylo <- alpha.fd.hill(com, dist_phylo, q=2, tau = "mean")

Hill <- data.frame(station=rownames(com), genet=Hill_gen$asb_FD_Hill, trait=Hill_phylo$asb_FD_Hill)
colnames(Hill) <- c("station", "genet", "phylo")

cor.test(Hill$phylo, Hill$genet, method = "spearman")


plot_alpha_phylo_gen <- ggplot(Hill, aes(genet, phylo))+
  geom_point()+
  labs(y= expression(paste("Phylogenetic ", alpha,"-diversity")), x=expression(paste("Sequence ", alpha,"-diversity")))+
  annotate(geom="text", x=1, y=14, label="Spearman rho=0.92 \n p<0.001", hjust=0, size=6, fontface = "bold")+
  theme_sleek(base_size = 24)+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14))

save(plot_alpha_phylo_gen, file="Rdata/plot_alpha_phylo_gen.rdata")


# calculate beta HILL for genet and phylo between samples

beta_hill_gen <- beta.fd.hill(com, dist_gen, q=2, tau = "mean", beta_type="Jaccard")
beta_hill_gen <- beta_hill_gen$beta_fd_q$q2

beta_hill_phylo <- beta.fd.hill(com, dist_phylo, q=2, tau = "mean", beta_type="Jaccard")
beta_hill_phylo <- beta_hill_phylo$beta_fd_q$q2


ecodist::mantel(beta_hill_gen ~ beta_hill_phylo)


co_rank <- coranking(beta_hill_gen, beta_hill_phylo, input_Xi = "dist")
NX <- coRanking::R_NX(co_rank)
AUC <- coRanking::AUC_ln_K(NX)


plot(beta_hill_phylo ~ beta_hill_gen, 
     xlab=expression(paste("Sequence ", beta,"-diversity")), 
     ylab=expression(paste("Phylogenetic  ", beta,"-diversity")),
     cex.lab = 1.2, cex.axis = 1.2)
