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

data <- left_join(exp_var_num, rich_station[,c("chondri_MOTUs", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

hist(data$chondri_MOTUs, main = "MOTUs_chondri", xlab ="MOTUs_chondri")

#data$chondri_MOTUs <- log10(data$chondri_MOTUs+1)
#hist(data$chondri_MOTUs, main = "log10(MOTUs_chondri+1)", xlab ="log10(MOTUs_chondri+1)")


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

mexp <- gls(chondri_MOTUs ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(chondri_MOTUs ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(chondri_MOTUs ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(chondri_MOTUs ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(chondri_MOTUs ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")


# Info full model
AIC(mexp, mgau, mlin, msph, mrat)

gls.full <- mgau

# remove colinear variables from VIF
gls.chondri <- gls(chondri_MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

AIC(gls.chondri)
summary(gls.chondri)
anova(gls.chondri, type ="marginal")

shapiro.test(gls.chondri$residuals)
hist(gls.chondri$residuals)

# R² for GLS
MOTU_pred <- predict(gls.chondri)
fit <- lm(MOTU_pred ~ data$chondri_MOTUs)
RsquareAdj(fit)



fit.grav.chondri <- visreg(gls.chondri,"Gravity",scale="response")
save(fit.grav.chondri, file="Rdata/fit.grav.chondri.rdata")
fit.SST.chondri <- visreg(gls.chondri,"mean_SST_1year",scale="response")
save(fit.SST.chondri, file="Rdata/fit.SST.chondri.rdata")
fit.SSS.chondri <- visreg(gls.chondri,"mean_sss_1year",scale="response")
save(fit.SSS.chondri, file="Rdata/fit.SSS.chondri.rdata")
fit.MED.chondri <- visreg(gls.chondri,"MarineEcosystemDependency",scale="response")
save(fit.MED.chondri, file="Rdata/fit.MED.chondri.rdata")
fit.CT.chondri <- visreg(gls.chondri,"dist_to_CT",scale="response")
save(fit.CT.chondri, file="Rdata/fit.CT.chondri.rdata")
fit.coast.chondri <- visreg(gls.chondri,"distCoast",scale="response")
save(fit.coast.chondri, file="Rdata/fit.coast.chondri.rdata")
fit.vol.chondri <- visreg(gls.chondri,"volume",scale="response")
save(fit.vol.chondri, file="Rdata/fit.vol.chondri.rdata")
fit.samp.chondri <- visreg(gls.chondri,"depth_sampling",scale="response")
save(fit.samp.chondri, file="Rdata/fit.samp.chondri.rdata")
fit.DHW.chondri <- visreg(gls.chondri,"mean_DHW_5year",scale="response")
save(fit.DHW.chondri, file="Rdata/fit.DHW.chondri.rdata")




#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_DHW_5year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year")]
geo_var <- data[, c("bathy", "dist_to_CT", "distCoast","depth_sampling")]
socio_var <- data[,c("HDI2019", "conflicts", "Corruption_mean", "Gravity", "MarineEcosystemDependency")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.chondri$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.106+0.259+0.055+0.055+0.023+0.007+0.046, 
                        geography=0.297+0.259+0.031+0.023+0.007+0.005, 
                        socioeconomy=0.116+0.031+0.05+0.259+0.055+0.055+0.023,
                        sampling=0.007+0.05+0.023+0.055+0.046+0.007+0.005)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("geography", "socioeconomy", "environment", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())




#### Prediction reference ####

# 0 human impact

noimpact <- data %>%
  dplyr::select("mean_DHW_1year", "mean_DHW_5year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year", "Corruption_mean", "HDI2019", "Gravity", "MarineEcosystemDependency", "conflicts", "dist_to_CT", "bathy", "depth_sampling", "distCoast", "volume")

noimpact$Corruption_mean <- 2.5
noimpact$HDI2019 <- 1
noimpact$Gravity <- 0
noimpact$MarineEcosystemDependency <- 0
noimpact$conflicts <- 0

chondri_noimpact <- expm1(predict(gls.chondri, noimpact))

rich_station$chondri_noimpact <- ceiling(chondri_noimpact)


# least regional human impact

leastimpact <- exp_var %>%
  dplyr::select("province", "mean_DHW_1year", "mean_DHW_5year", "mean_sss_1year", "mean_SST_1year", "mean_npp_1year", "Corruption_mean", "HDI2019", "Gravity", "MarineEcosystemDependency", "conflicts", "dist_to_CT", "bathy", "depth_sampling", "distCoast", "volume")

province <- as.character(unique(exp_var$province))

for (i in 1:length(province)) {
  df <- exp_var[exp_var$province == province[i],]
  leastimpact[leastimpact$province == province[i], "Gravity"] <- min(df$Gravity)
  leastimpact[leastimpact$province == province[i], "conflicts"] <- min(df$conflicts)
  leastimpact[leastimpact$province == province[i], "MarineEcosystemDependency"] <- min(df$MarineEcosystemDependency)
  leastimpact[leastimpact$province == province[i], "Corruption_mean"] <- max(df$Corruption_mean)
  leastimpact[leastimpact$province == province[i], "HDI2019"] <- max(df$HDI2019)
}


chondri_leastimpact <- expm1(predict(gls.chondri, leastimpact))

rich_station$chondri_leastimpact <- ceiling(chondri_leastimpact)

# MOTUs richness prediction by gls

rich_station$chondri_predicted <- ceiling(expm1(MOTU_pred))


save(rich_station, file="Rdata/richness_station.rdata")
