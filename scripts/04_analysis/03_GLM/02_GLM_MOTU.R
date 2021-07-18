library(vegan)
library(ade4)
library(DHARMa)
library(car)
library(visreg)
library(ecospat)
library(modEvA)

load("Rdata/richness_station.rdata")
load("Rdata/all_explanatory_variables.rdata")
load("Rdata/all_explanatory_variables_numeric.rdata")
rownames(rich_station) <- rich_station$station

df <- exp_var%>%
  select(-c("station", "province"))
df_num <- exp_var_num%>%
  select(-c("station", "province"))
MOTUs <- rich_station$MOTUs


#### full model ####
glm_full <- glm(MOTUs ~ ., data=df_num, family="poisson")

summary(glm_full)
anova(glm_full)
AIC(glm_full)

# check for colinearity and select variables
vif(glm_full)
mctest::imcdiag(glm_full, method="VIF")

df_sel <- df_num %>%
  select(-c("pH_mean", "NGO","latitude"))

#### GLM with selected variables ####
glm_sel <- glm(MOTUs ~ ., data=df_sel, family="poisson")

summary(glm_sel)
anova(glm_sel)
AIC(glm_sel)

vif(glm_sel)
mctest::imcdiag(glm_sel, method="VIF") # should be ok, but if not remove correlated variables

# check residuals : if weird test overdispersion
simulateResiduals(glm_sel, plot=TRUE)

testOverdispersion(glm_sel) # si significatif --> quasi-poisson

# visualise regression
visreg(glm_sel)


#### check spatial autocorrelation ####

meta <- read.csv("metadata/Metadata_eDNA_Pole2Pole_v4.csv", sep=";")
meta <- meta %>%
  distinct(station, .keep_all=T)
rownames(meta) <- meta$station
meta <- meta[rownames(rich_station),]
coor <- as.matrix(cbind(meta[,"longitude_start"], meta[,"latitude_start"]))

nb <- knn2nb(knearneigh(coor, 1)) 
lstw <- nb2listw((knn2nb(knearneigh(coor, k=1))))

# Compute Moran's I using residuals of model and also raw data
moran.test(glm_sel$residuals, lstw) 


#### Variation partitioning ####
env_var <- df_num[,c("mean_sss_1year", "mean_npp_1year", "mean_SST_1year", "pH_mean", "mean_DHW_1year")]
geo_var <- df_num[, c("dist_to_coast", "latitude_start", "depth_fin", "dist_to_CT", "province", "depth_sampling")]
samp_var <- df_num[, c("sample_method", "volume")]


varPart()
