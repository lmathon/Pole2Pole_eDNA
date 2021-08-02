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
    select(sequence, count_reads) %>%
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


# calculate dissimilarity matrix with jaccard (PA) or Bray-curtis (abund)
dist_jac_mo <- as.matrix(vegdist(com, method = "jaccard"))



save(dist_jac_mo, file="Rdata/Jaccard_MOTU_dissimilarity.rdata")

# plot dissimilarity matrix

dissplot(dist_jac_mo, method=NA, 
         upper_tri = TRUE, 
         lower_tri = FALSE, 
         reverse_columns=TRUE,
         main="Jaccard dissimilarity in MOTU composition",
         col=bluered(100))
