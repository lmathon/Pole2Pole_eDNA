library(tidyverse)
library(ape)
library(vegan)

load("Rdata/03-filter-data.Rdata")

# prepare community matrix

df <- df_filtered %>%
  subset(!is.na(station))

amplicon <- unique(df$definition)

stations <- unique(df$station)
st <- vector("list")
com <- data.frame(amplicon=amplicon)

for (i in 1:length(stations)) {
  st[[i]] <- df %>%
    subset(station == stations[i])%>%
    select(definition, count_reads) %>%
    group_by(definition) %>% 
    summarise_all(funs(sum))
  colnames(st[[i]]) <- c("amplicon", stations[[i]])
  com <- left_join(com, st[[i]], by="amplicon")
  
}

rownames(com) <- com$amplicon
com <- com[,-1]
com <- as.data.frame(t(com))
com[is.na(com)] <- 0

# calculate dissimilarity matrix with jaccard (PA) or Bray-curtis (abund)
dist_jac <- as.matrix(vegdist(com, method = "jaccard"))

dist_bc <- as.matrix(vegdist(com, method = "bray"))


save(dist_bc, file="Rdata/Bray_dissimilarity_station.rdata")
save(dist_jac, file="Rdata/Jaccard_dissimilarity_station.rdata")
