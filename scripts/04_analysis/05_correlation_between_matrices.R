library(ade4)
library(tidyverse)


load("Rdata/Jaccard_MOTU_dissimilarity.rdata")
load("Rdata/Jaccard_species_dissimilarity.rdata")
load("Rdata/MNTD_pairwise_station.rdata")

load("Rdata/geographic_distance_stations.rdata")

# select same stations in big geographic distance matrix
dist_km <- dist_km[rownames(dist_bc), colnames(dist_bc)]


dist_jac_mo <- as.dist(dist_jac_mo)
dist_jac_sp <- as.dist(dist_jac_sp)
dist_km <- as.dist(dist_km)
mntd <- as.dist(mntd)

mantel.rtest(dist_jac_mo, mntd) # Motu composition dissimilarity & mntd correlated
mantel.rtest(dist_jac_sp, mntd)
mantel.rtest(dist_km, mntd) # geographic distance & mntd correlated
mantel.rtest(dist_km, dist_jac_mo) # geographic distance & Motu composition dissimilarity correlated
mantel.rtest(dist_km, dist_jac_sp)