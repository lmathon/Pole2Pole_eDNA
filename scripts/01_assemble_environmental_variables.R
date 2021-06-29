library(tidyverse)


meta <- read.csv("c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/Metadata_eDNA_Pole2Pole.csv", sep=";")

met <- meta[,c("code_spygen", "station", "latitude_start", "longitude_start")]

sst <- read.csv("c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA-sst-spatio-temporal.csv", sep=",")
met_sst <- left_join(met, sst, by=c("latitude_start", "longitude_start"))                
met_sst <- met_sst %>%
  distinct(code_spygen, .keep_all = T)
write.csv(met_sst, "c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA_sst.csv", row.names = F)


ph <- read.csv("c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA-pH.csv", sep=",")
met_ph <- left_join(met, ph, by=c("latitude_start", "longitude_start"))
write.csv(met_ph, "c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA_ph.csv", row.names = F)

sss <- read.csv("c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA-sss-spatio-temporal.csv", sep = ",")
met_sss <- left_join(met, sss, by=c("latitude_start", "longitude_start"))
met_sss <- met_sss %>%
  distinct(code_spygen, .keep_all = T)
write.csv(met_sss, "c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA_sss.csv", row.names = F)


npp <- read.csv("c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA-npp-spatio-temporal.csv", sep = ",")
met_npp <- left_join(met, npp, by=c("latitude_start", "longitude_start"))
met_npp <- met_npp %>%
  distinct(code_spygen, .keep_all = T)
write.csv(met_npp, "c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA_npp.csv", row.names = F)


dhw <- read.csv("c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA-dhw-spatio-temporal.csv", sep = ",")
met_dhw <- left_join(met, dhw, by=c("latitude_start", "longitude_start"))
met_dhw <- met_dhw %>%
  distinct(code_spygen, .keep_all = T)
write.csv(met_dhw, "c:/Users/mathon/Desktop/linux/Pole2Pole_eDNA/metadata/eDNA_dhw.csv", row.names = F)

