library(tidyverse)
library(ape)
library(vegan)

load("Rdata/03-filter-data.Rdata")

# prepare community matrix

df <- df_filtered %>%
  subset(!is.na(station))

family <- unique(df$family_name_corrected)

stations <- unique(df$station)
st <- vector("list")
com <- data.frame(family=family)

for (i in 1:length(stations)) {
  st[[i]] <- df %>%
    subset(station == stations[i])%>%
    dplyr::select(family_name_corrected, count_reads) %>%
    group_by(family_name_corrected) %>% 
    summarise_all(funs(sum))
  colnames(st[[i]]) <- c("family", stations[[i]])
  com <- left_join(com, st[[i]], by="family")
  
}

com <- com %>%
  filter(!is.na(family))
rownames(com) <- com$family
com <- com[,-1]
com[is.na(com)] <- 0
save(com, file="Rdata/station_family_presence_absence.rdata")
com <- as.data.frame(t(com))


# calculate dissimilarity matrix with jaccard (PA) or Bray-curtis (abund)
dist_jac_fam <- as.matrix(vegdist(com, method = "jaccard"))



save(dist_jac_fam, file="Rdata/Jaccard_family_dissimilarity.rdata")

# plot dissimilarity matrix

dissplot(dist_jac_fam, method=NA, 
         upper_tri = TRUE, 
         lower_tri = FALSE, 
         reverse_columns=TRUE,
         main="Jaccard dissimilarity in Family composition",
         col=bluered(100))
