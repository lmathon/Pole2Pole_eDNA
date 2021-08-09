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
glm_nb <- glm.nb(MOTUs ~ . -latitude_start - longitude_start -NGO -pH_mean -NoViolence_mean, data=data)
summary(glm_nb)
pchisq(glm_nb$deviance,glm_nb$df.residual,lower.tail = F) # well adjusted to data

# check for colinearity and select variables
mctest::imcdiag(glm_nb, method="VIF")

# check residuals :and overdispersion
simulateResiduals(glm_nb, plot=TRUE, refit = T)
testDispersion(glm_nb, alternative = "greater") 


# check spatial autocorrelation 
# Compute Moran's I using residuals of model and also raw data
moran.test(glm_nb$residuals, lstw) # spatial autocorrelation in residuals


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
fit <- lm(MOTU_pred ~ MOTUs)
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

gls.final <- gls(MOTUs ~ mean_DHW_5year+mean_sss_1year+mean_SST_1year+HDI2019+neartt+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

AIC(gls.final)
summary(gls.final)
anova(gls.final)

# R² for GLS
MOTU_pred <- predict(gls.final)
fit <- lm(MOTU_pred ~ MOTUs)
RsquareAdj(fit)



visreg(gls.final,"mean_DHW_5year",scale="response")
visreg(gls.final,"mean_SST_1year",scale="response")
visreg(gls.final,"mean_sss_1year",scale="response")
visreg(gls.final,"HDI2019",scale="response")
visreg(gls.final,"neartt",scale="response")
visreg(gls.final,"Gravity",scale="response")
visreg(gls.final,"MarineEcosystemDependency",scale="response")
visreg(gls.final,"conflicts",scale="response")
visreg(gls.final,"dist_to_CT",scale="response")
visreg(gls.final,"volume",scale="response")



#### Variation partitioning ####
env_var <- data[,c("mean_DHW_5year","mean_SST_1year","mean_sss_1year")]
geo_var <- data[, c("dist_to_CT")]
socio_var <- data[,c("HDI2019","neartt","MarineEcosystemDependency","Gravity","conflicts")]
samp_var <- data[, c("volume")]

varpart <- varpart(gls.final$fitted, env_var, geo_var, socio_var, samp_var)
plot(varpart, digits = 2, Xnames = c('environment', 'geography', 'socio-economy', 'sampling'), bg = c('navy', 'tomato', 'yellow', 'lightgreen'))


# boxplot partition per variable type

partition <- data.frame(environment=0.122+0.3+0.021+0.2+0.103, 
                        geography=0.25+0.3+0.002, 
                        socioeconomy=0.165+0.049+0.3+0.021+0.2, 
                        sampling=0.028+0.002+0.2+0.103)

partition <- as.data.frame(t(partition))
partition$variables <- rownames(partition)
partition$variables2 <- factor(partition$variables, levels = c("environment", "socioeconomy", "geography", "sampling"))

ggplot(partition, aes(x=variables2,y = V1))+
  geom_col(width = 0.2)+
  xlab("Variable type")+
  ylab("cumulated variance explained")+
  theme(legend.position="none", panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"), panel.grid.major = element_blank())


#### GLS on log(MOTUs) ####
data$MOTUs <- log1p(data$MOTUs)
gls.final <- gls(MOTUs ~ mean_DHW_5year+mean_sss_1year+mean_SST_1year+HDI2019+neartt+Gravity+MarineEcosystemDependency+conflicts+dist_to_CT+volume, correlation = corGaus(form = ~longitude_start + latitude_start, nugget = TRUE), data = data,method="ML")

AIC(gls.final)
summary(gls.final)
anova(gls.final)

# R² for GLS
MOTU_pred <- predict(gls.final)
fit <- lm(MOTU_pred ~ MOTUs)
RsquareAdj(fit)
