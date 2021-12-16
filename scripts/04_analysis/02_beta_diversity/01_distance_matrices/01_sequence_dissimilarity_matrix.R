library(tidyverse)
library(ape)
library(picante)
library(seriation)
library(mFD)

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


#### calculate beta-FD between pairs of stations
com[com>0] <- 1
com <- as.matrix(com)
beta_FD <- beta.fd.hill(com, dist_gen, q=2, tau = "mean", beta_type = "Jaccard") 
beta_FD <- beta_FD$beta_fd_q$q2
beta_FD <- as.matrix(beta_FD)

save(beta_FD, file = "Rdata/beta_FD.rdata")










#######################################################################
# UNUSED
#######################################################################

#### calculate MNTD between pairs of stations ####

mntd <- comdistnt(com, dist_gen, abundance.weighted = FALSE, exclude.conspecifics = FALSE)
mntd <- as.matrix(mntd)


# save Rdata

save(mntd, file="Rdata/MNTD_pairwise_station.rdata")



##### plot distance matrix ####

dissplot(mntd, method=NA, 
         upper_tri = TRUE, 
         lower_tri = FALSE, 
         reverse_columns=TRUE,
         main="Mean Nearest Taxon Distance between stations",
         col=bluered(100))


