# Add lot to metadata NC 

meta <- read.table("data/swarm/nouvellecaledonie/metadata/keep_data/all_samples.csv", h=F, sep=";", stringsAsFactors = F)

lot1 <- c("200414_NB501473_A_L1-4_AIMI-219", 
          "200414_NB501473_A_L1-4_AIMI-220",
          "200414_NB501473_A_L1-4_AIMI-221",
          "200414_NB501473_A_L1-4_AIMI-222",
          "200416_NB501473_A_L1-4_AIMI-223",
          "200416_NB501473_A_L1-4_AIMI-224",
          "200416_NB501473_A_L1-4_AIMI-225",
          "200416_NB501473_A_L1-4_AIMI-226")

lot2 <- c("200608_SN6662_A_L001_AIMI-236",
          "200615_SN234_A_L001_AIMI-237",
          "200629_SN6662_A_L001_AIMI-243")

lot3 <- c("201012_SN6662_A_L001_AIMI-278")

lot4 <- c("210205_NB501473_A_L1-4_AIMI-315",
          "210205_NB501473_A_L1-4_AIMI-316",
          "210205_NB501850_A_L1-4_AIMI-317",
          "210205_NB501850_A_L1-4_AIMI-318")

meta <- meta %>%
  mutate(lot = case_when(
    V2 %in% lot1 ~ "LotA", 
    V2 %in% lot2 ~ "LotB", 
    V2 %in% lot3 ~ "LotC", 
    V2 %in% lot4 ~ "LotD"))

write.table(meta,file="data/swarm/nouvellecaledonie/metadata/all_samples.csv",sep=";",col.names = F, row.names = F)
