library(ade4)
library(tidyverse)
library(vegan)
library(seriation)
library(ape)


load("Rdata/Jaccard_MOTU_dissimilarity.rdata")
load("Rdata/Jaccard_species_dissimilarity.rdata")
load("Rdata/MNTD_pairwise_station.rdata")

load("Rdata/geographic_distance_stations.rdata")
load("Rdata/selected_environmental_variables.rdata")
load("Rdata/selected_socioeconomic_variables.rdata")

# select same stations in big geographic distance matrix
dist_km <- dist_km[rownames(mntd), colnames(mntd)]

dist_jac_mo <- as.dist(dist_jac_mo)
dist_jac_sp <- as.dist(dist_jac_sp)
dist_km <- as.dist(dist_km)
mntd <- as.dist(mntd)


# environmental distance matrix
dist_env <- vegdist(env_var2[,-1], "mahalanobis", na.rm = TRUE)
dist_socio <- vegdist(socio_var2[,-1], "mahalanobis", na.rm = TRUE)

save(dist_env, file="Rdata/environmental_distance_matrix.rdata")
save(dist_socio, file="Rdata/socioeconomic_distance_matrix.rdata")

# plot environmental distance matrix
dissplot(dist_env, method=NA, 
         upper_tri = TRUE, 
         lower_tri = FALSE, 
         reverse_columns=TRUE,
         main="Environmental distance between stations",
         col=bluered(100))

# plot socioeconomic distance matrix
dissplot(dist_socio, method=NA, 
         upper_tri = TRUE, 
         lower_tri = FALSE, 
         reverse_columns=TRUE,
         main="Socioeconomic distance between stations",
         col=bluered(100))


# calculate correlations
mantel.rtest(dist_jac_mo, mntd) # Motu composition dissimilarity & mntd correlated
mantel.rtest(dist_jac_mo, dist_km) # geographic distance & Motu composition dissimilarity correlated
mantel.rtest(dist_km, mntd) # geographic distance & mntd correlated
mantel.rtest(dist_env, mntd)
mantel.rtest(dist_socio, mntd)
mantel.rtest(dist_env, dist_jac_mo)
mantel.rtest(dist_socio, dist_jac_mo)


