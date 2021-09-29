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
library(ggplot2)


load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data <- left_join(exp_var_num, rich_station[,c("crypto_MOTUs", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

hist(data$crypto_MOTUs, main = "MOTUs_crypto", xlab ="MOTUs_crypto")

data$crypto_MOTUs <- log1p(data$crypto_MOTUs)
hist(data$crypto_MOTUs, main = "log(MOTUs_crypto)", xlab ="log(MOTUs_crypto)")

# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
identical(as.character(rownames(meta)), rownames(rich_station))
coor <- meta[, c("longitude_start", "latitude_start")]
data <- cbind(data, coor)


#### GLS to account for spatial autocorrelation ####

mexp <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(crypto_MOTUs ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, mlin, msph, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.crypto <- gls(crypto_MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

summary(gls.crypto)
anova(gls.crypto, type = "marginal")
AIC(gls.crypto)


# R² for GLS
MOTU_pred <- predict(gls.crypto)
fit <- lm(MOTU_pred ~ data$crypto_MOTUs)
RsquareAdj(fit)

shapiro.test(gls.crypto$residuals)
hist(gls.crypto$residuals)



fit.grav.crypto <- visreg(gls.crypto,"Gravity",scale="response")
save(fit.grav.crypto, file="Rdata/fit.grav.crypto.rdata")
fit.SST.crypto <- visreg(gls.crypto,"mean_SST_1year",scale="response")
save(fit.SST.crypto, file="Rdata/fit.SST.crypto.rdata")
fit.SSS.crypto <- visreg(gls.crypto,"mean_sss_1year",scale="response")
save(fit.SSS.crypto, file="Rdata/fit.SSS.crypto.rdata")
fit.MED.crypto <- visreg(gls.crypto,"MarineEcosystemDependency",scale="response")
save(fit.MED.crypto, file="Rdata/fit.MED.crypto.rdata")
fit.CT.crypto <- visreg(gls.crypto,"dist_to_CT",scale="response")
save(fit.CT.crypto, file="Rdata/fit.CT.crypto.rdata")
fit.coast.crypto <- visreg(gls.crypto,"distCoast",scale="response")
save(fit.coast.crypto, file="Rdata/fit.coast.crypto.rdata")
fit.vol.crypto <- visreg(gls.crypto,"volume",scale="response")
save(fit.vol.crypto, file="Rdata/fit.vol.crypto.rdata")
fit.samp.crypto <- visreg(gls.crypto,"depth_sampling",scale="response")
save(fit.samp.crypto, file="Rdata/fit.samp.crypto.rdata")
fit.DHW.crypto <- visreg(gls.crypto,"mean_DHW_5year",scale="response")
save(fit.DHW.crypto, file="Rdata/fit.DHW.crypto.rdata")



#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_DHW_5year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year")]
geo_var <- data[, c("bathy", "dist_to_CT", "distCoast","depth_sampling")]
socio_var <- data[,c("HDI2019", "conflicts", "Corruption_mean", "Gravity", "MarineEcosystemDependency")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.crypto$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.164+0.214+0.125+0.140+0.072+0.025+0.099, 
                        geography=0.130+0.214+0.125+0.025+0.004, 
                        socioeconomy=0.071+0.018+0.214+0.140+0.125+0.072, 
                        sampling=0.009+0.018+0.004+0.125+0.072+0.025+0.099)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "socioeconomy", "geography", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


