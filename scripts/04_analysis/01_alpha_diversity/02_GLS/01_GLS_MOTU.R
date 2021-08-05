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
library(nlme)
library(MuMIn)
library(rcompanion)
library(ggpubr)

load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data <- left_join(exp_var_num, rich_station[,c("MOTUs", "station")], by="station")
data <- data %>%
  select(-c(station))
#data <- data %>%
  #dplyr::select(c(MOTUs, mean_sss_1year, mean_SST_1year, Gravity, NGO, Naturalresourcesrents, dist_to_CT, bathy, latitude, distCoast, volume, sample_method, sequencer))

MOTUs <- rich_station$MOTUs

# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- meta[, c("longitude_start", "latitude_start")]
data <- cbind(data, coor)

#### full model -- GLM negative binomial ####
glm_nb <- glm.nb(MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+pH_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+NGO+Naturalresourcesrents+dist_to_CT+bathy+depth_sampling+latitude+distCoast+volume+sample_method+sequencer, data=data, family="poisson")
summary(glm_nb)
pchisq(glm_nb$deviance,glm_nb$df.residual,lower.tail = F) # well adjusted to data

# check for colinearity and select variables
mctest::imcdiag(glm_full, method="VIF")

# check residuals :and overdispersion
simulateResiduals(glm_nb, plot=TRUE, refit = T)
testDispersion(glm_nb, alternative = "greater") 


# check spatial autocorrelation 
# Compute Moran's I using residuals of model and also raw data
moran.test(glm_nb$residuals, lstw) # spatial autocorrelation in residuals


#### GLS to account for spatial autocorrelation ####

mexp <- gls(MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+pH_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+NGO+Naturalresourcesrents+dist_to_CT+bathy+depth_sampling+latitude+distCoast+volume, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+pH_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+NGO+Naturalresourcesrents+dist_to_CT+bathy+depth_sampling+latitude+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+pH_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+NGO+Naturalresourcesrents+dist_to_CT+bathy+depth_sampling+latitude+distCoast+volume, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+pH_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+NGO+Naturalresourcesrents+dist_to_CT+bathy+depth_sampling+latitude+distCoast+volume, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+pH_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+NGO+Naturalresourcesrents+dist_to_CT+bathy+depth_sampling+latitude+distCoast+volume, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
            

# select best model
AIC(mexp, mgau, msph, mlin, mrat)

gls.full <- mgau
summary(gls.full)
Anova(gls.full)

stepAIC(gls.full)

gls.final <- gls(MOTUs ~ mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+HDI2019+neartt+Gravity+MarineEcosystemDependency+NGO+bathy+depth_sampling, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
AIC(gls.final)

summary(gls.final)

# calculate Rsquare
nagelkerke(gls.final)



Anova(gls.final)

plot(gls.final)

visreg(gls.final,"mean_DHW_5year",scale="response")
visreg(gls.final,"mean_SST_1year",scale="response")
visreg(gls.final,"mean_sss_1year",scale="response")
visreg(gls.final,"mean_npp_1year",scale="response")
visreg(gls.final,"HDI2019",scale="response")
visreg(gls.final,"neartt",scale="response")
visreg(gls.final,"Gravity",scale="response")
visreg(gls.final,"MarineEcosystemDependency",scale="response")
visreg(gls.final,"NGO",scale="response")
visreg(gls.final,"bathy",scale="response")
visreg(gls.final,"depth_sampling",scale="response")



#### Variation partitioning ####
env_var <- data[,c("mean_DHW_5year","mean_SST_1year","mean_npp_1year","mean_sss_1year")]
geo_var <- data[, c("bathy","depth_sampling")]
socio_var <- data[,c("HDI2019","neartt","MarineEcosystemDependency","Gravity","NGO")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.final$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))
