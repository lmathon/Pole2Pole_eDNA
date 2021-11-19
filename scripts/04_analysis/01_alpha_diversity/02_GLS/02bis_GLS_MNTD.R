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

load("Rdata/MNTD_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(mntd_stations) <- mntd_stations$station

data <- left_join(exp_var_num, mntd_stations[,c("MNTD", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

hist(data$MNTD, main = "MNTD", xlab ="MNTD")

# join longitude & latitude
meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(mntd_stations),]
identical(as.character(rownames(meta)), rownames(mntd_stations))
coor <- meta[, c("longitude_start", "latitude_start", "province")]
data <- cbind(data, coor)


#### GLS to account for spatial autocorrelation ####

mexp <- gls(MNTD ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(MNTD ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(MNTD ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(MNTD ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(MNTD ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, msph, mlin, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.MNTD <- gls(MNTD ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

AIC(gls.MNTD)
summary(gls.MNTD)
anova(gls.MNTD, type = "marginal")

# R² for GLS
MNTD_pred <- predict(gls.MNTD)
fit <- lm(MNTD_pred ~ data$MNTD)
RsquareAdj(fit)

hist(gls.MNTD$residuals)


fit.grav.MNTD <- visreg(gls.MNTD,"Gravity",scale="response")
save(fit.grav.MNTD, file="Rdata/fit.grav.MNTD.rdata")
fit.SST.MNTD <- visreg(gls.MNTD,"mean_SST_1year",scale="response")
save(fit.SST.MNTD, file="Rdata/fit.SST.MNTD.rdata")
fit.SSS.MNTD <- visreg(gls.MNTD,"mean_sss_1year",scale="response")
save(fit.SSS.MNTD, file="Rdata/fit.SSS.MNTD.rdata")
fit.MED.MNTD <- visreg(gls.MNTD,"MarineEcosystemDependency",scale="response")
save(fit.MED.MNTD, file="Rdata/fit.MED.MNTD.rdata")
fit.CT.MNTD <- visreg(gls.MNTD,"dist_to_CT",scale="response")
save(fit.CT.MNTD, file="Rdata/fit.CT.MNTD.rdata")
fit.coast.MNTD <- visreg(gls.MNTD,"distCoast",scale="response")
save(fit.coast.MNTD, file="Rdata/fit.coast.MNTD.rdata")
fit.vol.MNTD <- visreg(gls.MNTD,"volume",scale="response")
save(fit.vol.MNTD, file="Rdata/fit.vol.MNTD.rdata")
fit.samp.MNTD <- visreg(gls.MNTD,"depth_sampling",scale="response")
save(fit.samp.MNTD, file="Rdata/fit.samp.MNTD.rdata")
fit.DHW.MNTD <- visreg(gls.MNTD,"mean_DHW_5year",scale="response")
save(fit.DHW.MNTD, file="Rdata/fit.DHW.MNTD.rdata")


fit.grav_med.MNTD <- visreg2d(gls.MNTD, "Gravity", "MarineEcosystemDependency", scale = "response", type = "conditional", main="log10(MNTD richness +1)", xlab="log10(Gravity +1)")
save(fit.grav_med.MNTD, file="Rdata/fit.grav_med.MNTD.rdata")

#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_DHW_5year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year")]
geo_var <- data[, c("bathy", "dist_to_CT", "distCoast","depth_sampling")]
socio_var <- data[,c("HDI2019", "conflicts", "Corruption_mean", "Gravity", "MarineEcosystemDependency")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.MNTD$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.135+0.0155+0.315+0.006+0.104, 
                        geography=0.089+0.155+0.019+0.006+0.336, 
                        socioeconomy=0.025+0.055+0.019+0.006+0.315, 
                        sampling=0.055+0.104+0.104+0.006+0.336)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("geography",  "sampling", "environment", "socioeconomy"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


#### effect size ####

MNTD_effectsize <- effectsize(gls.MNTD)
MNTD_effectsize <- MNTD_effectsize[-1,]
MNTD_effectsize$taxa <- "MNTD"
MNTD_effectsize$vargroup <- c("environment","environment","environment","environment","environment","socio","socio","socio","socio","socio","geography","geography","geography","geography","sampling")

save(MNTD_effectsize, file = "Rdata/MNTD_effectsize.rdata")
