library(vegan)
library(ade4)
library(DHARMa)
library(car)
library(visreg)
library(ecospat)
library(modEvA)
library(psych)
library(MASS)
library(AER)
library(pscl)
library(spdep)
library(tidyverse)

load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data <- left_join(exp_var, rich_station, by="station")
data <- data %>%
  dplyr::select(c(MOTUs, mean_sss_1year, mean_SST_1year, Gravity, NGO, Naturalresourcesrents, dist_to_CT, bathy, latitude, distCoast, volume, sample_method, sequencer))

MOTUs <- rich_station$MOTUs


#### full model ####
glm_full <- glm(MOTUs ~ ., data=data, family="poisson")
summary(glm_full)
anova(glm_full)


# check for colinearity and select variables
mctest::imcdiag(glm_full, method="VIF")

# select best model
glm_select=step(glm_full)
summary(glm_select)
dispersiontest(glm_select) # overdispersed
pchisq(glm_select$deviance,glm_select$df.residual,lower.tail = F) # not well adjusted to data

# GLM negative binomial
glm_nb <- glm.nb(data=data, MOTUs ~ . -Gravity)
summary(glm_nb)
pchisq(glm_nb$deviance,glm_nb$df.residual,lower.tail = F) # well adjusted to data


# check residuals :and overdispersion
simulateResiduals(glm_nb, plot=TRUE, refit = T)
testDispersion(glm_nb, alternative = "greater") 


#### check spatial autocorrelation ####

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- as.matrix(cbind(meta[,"longitude_start"], meta[,"latitude_start"]))

nb <- knn2nb(knearneigh(coor, 1, longlat = TRUE)) 
lstw <- nb2listw(nb)

# Compute Moran's I using residuals of model and also raw data
moran.test(glm_nb$residuals, lstw) # spatial autocorrelation in residuals

# create variable of autocorrelation
bw <- max(unlist(nbdists(nb, coor)))

autocor <- autocov_dist(MOTUs, nbs=bw, coor, longlat = T, zero.policy = T)
autocor[which(!is.finite(autocor))] <- 0

# regress autocorrelation variable with all other variables
lm_autocor <- lm(autocor~ mean_sss_1year+mean_SST_1year+NGO+Naturalresourcesrents+dist_to_CT+bathy+latitude+distCoast+volume+sample_method+sequencer, data=data)
autocorVar <- residuals(lm_autocor)


## GLM with autocorVar
glm_nb2 <- glm.nb(data=data, MOTUs ~ scale(autocorVar)+mean_sss_1year+mean_SST_1year+NGO+Naturalresourcesrents+dist_to_CT+bathy+latitude+distCoast+volume+sample_method+sequencer)
summary <- summary(glm_nb2)

pchisq(glm_nb2$deviance,glm_nb2$df.residual,lower.tail = F) # well adjusted to data
simulateResiduals(glm_nb2, plot=TRUE, refit = T)
testDispersion(glm_nb2, alternative = "greater") 


#### Variation partitioning ####
env_var <- data[,c("mean_sss_1year", "mean_SST_1year")]
geo_var <- data[, c("distCoast", "latitude", "bathy", "dist_to_CT")]
socio_var <- data[,c("NGO", "Naturalresourcesrents")]
samp_var <- data[, c("volume")]

varpart <- varpart(glm_nb2$fitted.values, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))
