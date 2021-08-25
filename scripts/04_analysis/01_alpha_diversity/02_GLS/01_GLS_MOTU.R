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

data <- left_join(exp_var_num, rich_station[,c("MOTUs", "station")], by="station")
data <- data %>%
  dplyr::select(-c(station))

data$MOTUs <- log1p(data$MOTUs)

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

mexp <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corExp(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mgau <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

msph <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corSpher(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mlin <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corLin(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

mrat <- gls(MOTUs ~ . -latitude_start - longitude_start, correlation = corRatio(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
            

# Info full model
AIC(mexp, mgau, msph, mlin, mrat)

gls.full <- mgau
summary(gls.full)
anova(gls.full)
AIC(gls.full)

# R² for GLS
MOTU_pred <- predict(gls.full)
fit <- lm(MOTU_pred ~ data$MOTUs)
RsquareAdj(fit)

plot(gls.full$residuals)


env_var <- data[,c("mean_DHW_5year","mean_SST_1year","mean_npp_1year","mean_sss_1year", "mean_DHW_5year", "pH_mean")]
geo_var <- data[, c("bathy","depth_sampling", "dist_to_CT", "distCoast")]
socio_var <- data[,c("NoViolence_mean", "Corruption_mean", "HDI2019","neartt","MarineEcosystemDependency","Gravity","NGO", "conflicts")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.full$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# remove colinear variables from VIF
gls.final <- gls(MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
stepAIC(gls.final)

gls.final <- gls(MOTUs ~ mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

AIC(gls.final)
summary(gls.final)
anova(gls.final)

# R² for GLS
MOTU_pred <- predict(gls.final)
fit <- lm(MOTU_pred ~ data$MOTUs)
RsquareAdj(fit)

hist(gls.final$residuals)

visreg(gls.final,"mean_DHW_5year",scale="response")
visreg(gls.final,"mean_SST_1year",scale="response")
visreg(gls.final,"mean_sss_1year",scale="response")
visreg(gls.final,"mean_npp_1year",scale="response")
visreg(gls.final,"Corruption_mean",scale="response")
visreg(gls.final,"HDI2019",scale="response")
visreg(gls.final,"neartt",scale="response")
visreg(gls.final,"Gravity",scale="response")
visreg(gls.final,"MarineEcosystemDependency",scale="response")
visreg(gls.final,"conflicts",scale="response")
visreg(gls.final,"dist_to_CT",scale="response")
visreg(gls.final,"depth_sampling",scale="response")
visreg(gls.final,"distCoast",scale="response")
visreg(gls.final,"volume",scale="response")



#### Variation partitioning ####
env_var <- data[,c("mean_DHW_5year","mean_SST_1year","mean_sss_1year", "mean_npp_1year")]
geo_var <- data[, c("dist_to_CT", "depth_sampling","distCoast" )]
socio_var <- data[,c("Corruption_mean","HDI2019","neartt","MarineEcosystemDependency","Gravity","conflicts")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.final$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.124+0.083+0.109+0.075+0.034+0.149+0.214, 
                        geography=0.075+0.214+0.149+0.034+0.045+0.027, 
                        socioeconomy=0.046+0.027+0.214+0.083+0.109+0.149+0.032, 
                        sampling=0.032+0.149+0.109+0.075+0.034+0.045+0.031)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "socioeconomy", "geography", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())

