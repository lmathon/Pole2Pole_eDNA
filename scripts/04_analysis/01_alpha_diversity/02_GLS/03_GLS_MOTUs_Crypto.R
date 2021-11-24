library(funk)
library(vegan)
library(ade4)
library(DHARMa)
library(car)
library(visreg)
library(ecospat)
library(modEvA)
library(psych)
library(MASS)
library(lmeInfo)
library(pscl)
library(spdep)
library(tidyverse)
library(nlme)
library(MuMIn)
library(rcompanion)
library(ggpubr)
library(ggplot2)
library(effectsize)


load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

data <- left_join(exp_var_num, rich_station[,c("crypto_MOTUs", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

hist(data$crypto_MOTUs, main = "MOTUs_crypto", xlab ="MOTUs_crypto")

data$crypto_MOTUs <- log10(data$crypto_MOTUs+1)
hist(data$crypto_MOTUs, main = "log10(MOTUs_crypto+1)", xlab ="log10(MOTUs_crypto+1)")

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
gls.crypto <- gls(crypto_MOTUs ~ mean_DHW_1year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Voice_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

save(gls.crypto, file="Rdata/gls_crypto.rdata")

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
fit.DHW.crypto <- visreg(gls.crypto,"mean_DHW_1year",scale="response")
save(fit.DHW.crypto, file="Rdata/fit.DHW.crypto.rdata")

fit.grav_med.crypto <- visreg2d(gls.crypto, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(Crypto richness +1)", xlab="lg10(Gravity +1)")
save(fit.grav_med.crypto, file="Rdata/fit.grav_med.crypto.rdata")


#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year")]
geo_var <- data[, c("bathy", "dist_to_CT", "distCoast","depth_sampling")]
socio_var <- data[,c("HDI2019", "conflicts", "Voice_mean", "Gravity", "MarineEcosystemDependency")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.crypto$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.142+0.233+0.148+0.118+0.155+0.047, 
                        geography=0.159+0.233+0.155+0.012, 
                        socioeconomy=0.090+0.024+0.233+0.155+0.148+0.118, 
                        sampling=0.006+0.024+0.012+0.155+0.118+0.047)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "socioeconomy", "geography", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())

#### effect size ####

crypto_effectsize <- effectsize(gls.crypto)
crypto_effectsize <- crypto_effectsize[-1,]
crypto_effectsize$taxa <- "Richness - Cryptobenthics"
crypto_effectsize$vargroup <- c("environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")

save(crypto_effectsize, file = "Rdata/crypto_effectsize.rdata")

