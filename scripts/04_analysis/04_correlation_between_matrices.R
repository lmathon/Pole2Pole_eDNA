library(ade4)
library(tidyverse)

load("Rdata/Bray_dissimilarity_station.rdata")
load("Rdata/Jaccard_dissimilarity_station.rdata")
load("Rdata/MNTD_pairwise_station.rdata")

load("Rdata/geographic_distance_stations.rdata")

# select same stations in big geographic distance matrix
dist_km <- dist_km[, colnames(dist_bc)]
dist_km <- dist_km[rownames(dist_bc),]


dist_jac <- as.dist(dist_jac)
dist_km <- as.dist(dist_km)
mntd <- as.dist(mntd)

mantel.rtest(dist_jac, mntd) # Motu composition dissimilarity & mntd correlated
mantel.rtest(dist_km, mntd) # geographic distance & mntd correlated
mantel.rtest(dist_km, dist_jac) # geographic distance & Motu composition dissimilarity correlated
