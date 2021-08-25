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

hist(data$crypto_MOTUs)

data$crypto_MOTUs <- log1p(data$crypto_MOTUs)

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
summary(gls.full)
anova(gls.full)
AIC(gls.full)

# R² for GLS
MOTU_pred <- predict(gls.full)
fit <- lm(MOTU_pred ~ data$crypto_MOTUs)
RsquareAdj(fit)




# remove colinear variables from VIF
gls.final <- gls(crypto_MOTUs ~ mean_DHW_1year+mean_DHW_5year+mean_sss_1year+mean_SST_1year+mean_npp_1year+Corruption_mean+HDI2019+neartt+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+bathy+depth_sampling+distCoast+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")
stepAIC(gls.final)

gls.final <- gls(crypto_MOTUs ~ mean_DHW_1year+mean_DHW_5year+HDI2019+neartt+conflicts+depth_sampling+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

AIC(gls.final)
summary(gls.final)
anova(gls.final)

shapiro.test(gls.final$residuals)
hist(gls.final$residuals)

# R² for GLS
MOTU_pred <- predict(gls.final)
fit <- lm(MOTU_pred ~ data$crypto_MOTUs)
RsquareAdj(fit)



visreg(gls.final,"mean_DHW_1year",scale="response")
visreg(gls.final,"mean_DHW_5year",scale="response")
visreg(gls.final,"HDI2019",scale="response")
visreg(gls.final,"neartt",scale="response")
visreg(gls.final,"volume",scale="response")
visreg(gls.final,"conflicts",scale="response")
visreg(gls.final,"depth_sampling",scale="response")



#### Variation partitioning ####
env_var <- data[,c("mean_DHW_1year", "mean_DHW_5year")]
geo_var <- data[, c("depth_sampling")]
socio_var <- data[,c("HDI2019","neartt", "conflicts")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.final$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.043+0.243+0.005+0.001, 
                        geography=0.005+0.005+0.001+0.021+0.07, 
                        socioeconomy=0.426+0.219+0.07+0.021+0.005+0.243+0.001, 
                        sampling=0.219+0.021+0.001+0.046)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("socioeconomy", "environment", "sampling", "geography"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


