library(tidyverse)
library(ape)
library(vegan)

load("Rdata/03-filter-data.Rdata")

# prepare community matrix

df <- df_filtered %>%
  subset(!is.na(station))

species <- unique(df$species_name_corrected)

stations <- unique(df$station)
st <- vector("list")
com <- data.frame(species=species)

for (i in 1:length(stations)) {
  st[[i]] <- df %>%
    subset(station == stations[i])%>%
    select(species_name_corrected, count_reads) %>%
    group_by(species_name_corrected) %>% 
    summarise_all(funs(sum))
  colnames(st[[i]]) <- c("species", stations[[i]])
  com <- left_join(com, st[[i]], by="species")
  
}
com <- com %>%
  filter(!is.na(species))
rownames(com) <- com$species
com <- com[,-1]
com[is.na(com)] <- 0
save(com, file="Rdata/station_species_presence_absence.rdata")
com <- as.data.frame(t(com))


# calculate dissimilarity matrix with jaccard (PA) or Bray-curtis (abund)
dist_jac_sp <- as.matrix(vegdist(com, method = "jaccard"))



save(dist_jac_sp, file="Rdata/Jaccard_species_dissimilarity.rdata")

# plot dissimilarity matrix

dissplot(dist_jac_sp, method=NA, 
         upper_tri = TRUE, 
         lower_tri = FALSE, 
         reverse_columns=TRUE,
         main="Jaccard dissimilarity in species composition",
         col=bluered(100))
