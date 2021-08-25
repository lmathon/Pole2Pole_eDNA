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

hist(data$chondri_MOTUs)

data$chondri_MOTUs <- log1p(data$chondri_MOTUs)


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

gls.full <- msph
summary(gls.full)
anova(gls.full)
AIC(gls.full)

# R² for GLS
MOTU_pred <- predict(gls.full)
fit <- lm(MOTU_pred ~ data$chondri_MOTUs)
RsquareAdj(fit)

plot(gls.full$residuals)


# remove colinear variables from VIF
gls.final <- gls(chondri_MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
stepAIC(gls.final)

gls.final <- gls(chondri_MOTUs ~ mean_sss_1year+mean_npp_1year+HDI2019+neartt+Gravity+MarineEcosystemDependency+conflicts+bathy+depth_sampling, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

AIC(gls.final)
summary(gls.final)
anova(gls.final)

shapiro.test(gls.final$residuals)
plot(gls.final$residuals ~ gls.final$fitted)
hist(gls.final$residuals)

# R² for GLS
MOTU_pred <- predict(gls.final)
fit <- lm(MOTU_pred ~ data$chondri_MOTUs)
RsquareAdj(fit)



visreg(gls.final,"mean_sss_1year",scale="response")
visreg(gls.final,"mean_npp_1year",scale="response")
visreg(gls.final,"HDI2019",scale="response")
visreg(gls.final,"neartt",scale="response")
visreg(gls.final,"Gravity",scale="response")
visreg(gls.final,"MarineEcosystemDependency",scale="response")
visreg(gls.final,"conflicts",scale="response")
visreg(gls.final,"bathy",scale="response")
visreg(gls.final,"depth_sampling",scale="response")



#### Variation partitioning ####
env_var <- data[,c("mean_npp_1year", "mean_sss_1year")]
geo_var <- data[, c("depth_sampling", "bathy")]
socio_var <- data[,c("HDI2019","neartt","Gravity","MarineEcosystemDependency","conflicts")]


varpart <- varpart(gls.final$fitted, env_var, geo_var, socio_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy'), bg = c('navy', 'tomato', 'yellow'))


# boxplot partition per variable type

partition <- data.frame(environment=0.145+0.003, 
                        geography=0.056+0.003+0.285, 
                        socioeconomy=0.547+0.285)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("socioeconomy", "geography", "environment"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


