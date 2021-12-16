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
com[com>0] <- 1



#### calculate alpha FD Hill ####
com <- as.matrix(com)
test <- alpha.fd.hill(com, dist_gen, q=c(0,1,2), tau = "mean")
FD_Hill <- as.data.frame(test$asb_FD_Hill)
FD_Hill$station <- rownames(FD_Hill)

save(FD_Hill, file = "Rdata/FD_Hill_alpha.rdata")








######################################################################################
# UNUSED
######################################################################################



#### calculate MNTD for each stations ####

mntd_stations <- ses.mntd(com, dist_gen, abundance.weighted = FALSE, null.model = "independentswap", runs = 999, iterations = 100)

mntd_stations <- data.frame(station=stations, MNTD=mntd_stations$mntd.obs.z)
save(mntd_stations, file="Rdata/MNTD_station.rdata")

mntd_stations <- mntd(com, dist_gen, abundance.weighted = FALSE)
mntd_stations_nc <- data.frame(station=stations, MNTD=mntd_stations)
save(mntd_stations_nc, file="Rdata/MNTD_station_NC.rdata")

