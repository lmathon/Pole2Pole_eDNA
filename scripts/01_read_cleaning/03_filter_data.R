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
  #filter(family_name_corrected %ni% "Salmonidae") %>%
  #filter(site_name!="kirkenes") %>%
  #filter(site_name %ni% c("mallorca", "menorca"))
  
list_read_step4 <- lapply(list_read_step4, function(x){
    x %>%
      filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
      filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
      filter(habitat=="marine")%>%
      filter(sample_method !="control") %>%
      filter(depth_sampling < 40) %>%
      filter(site_name != "Kirkenes") #%>%
      #filter(family_name_corrected %ni% "Salmonidae") %>%
      #filter(site_name!="kirkenes") %>%
      #filter(site_name %ni% c("mallorca", "menorca"))
})

save(list_read_step3, list_read_step3_0PCR, list_read_step4, df_filtered, file = "Rdata/03-filter-data.Rdata")
