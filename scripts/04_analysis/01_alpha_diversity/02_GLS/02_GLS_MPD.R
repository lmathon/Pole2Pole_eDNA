library(funk)
library(vegan)
library(ade4)
library(lme4)
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

load("Rdata/MPD_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(mpd_stations) <- mpd_stations$station

data <- left_join(exp_var_num, mpd_stations[,c("MPD", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

hist(data$MPD, main = "MPD", xlab ="MPD")

# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(mpd_stations),]
identical(as.character(rownames(meta)), rownames(mpd_stations))
coor <- meta[, c("longitude_start", "latitude_start", "province")]
data <- cbind(data, coor)


#### GLS to account for spatial autocorrelation ####

mexp <- gls(MPD ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(MPD ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(MPD ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(MPD ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(MPD ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, msph, mlin, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.MPD <- gls(MPD ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

AIC(gls.MPD)
summary(gls.MPD)
anova(gls.MPD, type = "marginal")

# R² for GLS
MPD_pred <- predict(gls.MPD)
fit <- lm(MPD_pred ~ data$MPD)
RsquareAdj(fit)

hist(gls.MPD$residuals)


fit.grav.MPD <- visreg(gls.MPD,"Gravity",scale="response")
save(fit.grav.MPD, file="Rdata/fit.grav.MPD.rdata")
fit.SST.MPD <- visreg(gls.MPD,"mean_SST_1year",scale="response")
save(fit.SST.MPD, file="Rdata/fit.SST.MPD.rdata")
fit.SSS.MPD <- visreg(gls.MPD,"mean_sss_1year",scale="response")
save(fit.SSS.MPD, file="Rdata/fit.SSS.MPD.rdata")
fit.MED.MPD <- visreg(gls.MPD,"MarineEcosystemDependency",scale="response")
save(fit.MED.MPD, file="Rdata/fit.MED.MPD.rdata")
fit.CT.MPD <- visreg(gls.MPD,"dist_to_CT",scale="response")
save(fit.CT.MPD, file="Rdata/fit.CT.MPD.rdata")
fit.coast.MPD <- visreg(gls.MPD,"distCoast",scale="response")
save(fit.coast.MPD, file="Rdata/fit.coast.MPD.rdata")
fit.vol.MPD <- visreg(gls.MPD,"volume",scale="response")
save(fit.vol.MPD, file="Rdata/fit.vol.MPD.rdata")
fit.samp.MPD <- visreg(gls.MPD,"depth_sampling",scale="response")
save(fit.samp.MPD, file="Rdata/fit.samp.MPD.rdata")
fit.DHW.MPD <- visreg(gls.MPD,"mean_DHW_5year",scale="response")
save(fit.DHW.MPD, file="Rdata/fit.DHW.MPD.rdata")


fit.grav_med.MPD <- visreg2d(gls.MPD, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(MPD richness +1)", xlab="log10(Gravity +1)")
save(fit.grav_med.MPD, file="Rdata/fit.grav_med.MPD.rdata")

#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_DHW_5year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year")]
geo_var <- data[, c("bathy", "dist_to_CT", "distCoast","depth_sampling")]
socio_var <- data[,c("HDI2019", "conflicts", "Corruption_mean", "Gravity", "MarineEcosystemDependency")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.MPD$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.301+0.075+0.027+0.283+0.045+0.042, 
                        geography=0.103+0.027+0.024+0.075+0.042+0.206, 
                        socioeconomy=0.027+0.027+0.024+0.283, 
                        sampling=0.026+0.206+0.042+0.045)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment","geography", "socioeconomy",  "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


#### effect size ####

MPD_effectsize <- effectsize(gls.MPD)
MPD_effectsize <- MPD_effectsize[-1,]
MPD_effectsize$taxa <- "ses.MPD"
MPD_effectsize$vargroup <- c("environment","environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")

save(MPD_effectsize, file = "Rdata/mpd_effectsize.rdata")
