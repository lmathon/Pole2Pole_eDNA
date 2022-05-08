library(tidyverse)
'%ni%' <- Negate("%in%")

load("Rdata/02-clean-data.Rdata")
df_all_filters$depth_sampling <- as.numeric(df_all_filters$depth_sampling)

df_filtered <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(habitat=="marine")%>%
  filter(sample_method !="control") %>%
  filter(depth_sampling < 40) %>%
  filter(site_name != "Kirkenes")
  
  
list_read_step4 <- lapply(list_read_step4, function(x){
    x %>%
      filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
      filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
      filter(habitat=="marine")%>%
      filter(sample_method !="control") %>%
      filter(depth_sampling < 40) %>%
      filter(site_name != "Kirkenes") #%>%
      
})

save(list_read_step3, list_read_step3_0PCR, list_read_step4, df_filtered, file = "Rdata/03-filter-data.Rdata")



list_stat <- split(df_filtered, df_filtered$station)

nb_samp <- lapply(list_stat, function(x){
     length(unique(x$sample_name_all_pcr))
 })

nb_samp <- as.data.frame(t(bind_rows(nb_samp)))
mean(nb_samp$V1)
sd(nb_samp$V1)