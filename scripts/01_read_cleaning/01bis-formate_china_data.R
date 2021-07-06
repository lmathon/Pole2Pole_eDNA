library(tidyverse)
library(data.table)
library(reshape)
library(dplyr)

# Source functions
source("scripts/01_read_cleaning/00_functions.R")
load("Rdata/archive_class_ncbi.Rdata")

# load and formate df china
df_china <- read.csv("data/swarm/chine/otutable_formated.csv", sep=";")
vars <- colnames(df_china[,1:30])
df_melt <- melt(df_china, id.vars=vars, variable_name = "sample_name")

colnames(df_melt) <- c("sequence", "id", "definition", "best_match.db_embl_std", "count", "family", "family_name", "genus", "genus_name" , "match_count.db_embl_std", "order", "order_name", "rank", "scientific_name", "species", "species_list.db_embl_std", "species_name", "taxid", "taxid_by_db.db_embl_std", "OTU", "total", "cloud", "length", "abundance", "chimera", "best_identity_database", "spread", "identity", "taxonomy", "references", "sample_name", "count_reads")

columns_to_remove <- c("amplicon", "family", "genus", "order", "species", "taxid", "OTU", "total", 
                       "cloud", "length", "abundance", "spread", "identity", "taxonomy", "references")

df_short <- df_melt %>%
  select(-one_of(columns_to_remove))

df_short <- df_short %>% relocate(best_identity_database, .after = last_col())

df_short$plaque <- NA
df_short$run <- NA
df_short$project <- "chine"
df_short$marker <- "teleo"
df_short$somme_tot <- NA
df_short$seuil <- NA
df_short$project_i <- "chine"

# clean taxonomy
df_taxo <- clean_taxonomy(df_short)

df_taxo$class_name <- NA

##fonctionne pas chez moi !!
  #df_output <- add_class_name_archive(df_taxo, archive_class_ncbi)
  #file_taxo_all <- list_output[[1]]
  #archive_class_ncbi <- list_output[[2]]

df_taxo$length_sequence <- df_china$length
df_taxo$n_PCR <- NA
df_taxo$sample_name_all_pcr <- df_taxo$sample_name



# assemble with metadata

metadata_field <- read.csv("metadata/Metadata_eDNA_Pole2Pole.csv", sep=";", stringsAsFactors = F)

columns_delete_field_metadata <- c("turbidity", "gps_start", "gps_b", "lat_gps_b", "long_gps_b", "gps_c", "long_gps_c", "lat_gps_d", "gps_half_turn", "longitude_turn", "latitude_end", "longitude_end", 
                                   "gps_end", "long_gps_d", "gps_d", "lat_gps_c", "latitude_turn", "data_manager", "gps_owner", "project")
metadata_field <- select(metadata_field, -c(columns_delete_field_metadata))

df_china_fin <- left_join(df_taxo, metadata_field, by=c("sample_name_all_pcr" = "code_spygen"))

save(df_china_fin, file="Rdata/df_china_formated.rdata")
