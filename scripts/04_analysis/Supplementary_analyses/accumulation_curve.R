# Accumulation curves on genus, family, order on eDNA data


# Lib 
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(conflicted)

# data
load("Rdata/03-filter-data.Rdata")

# 
'%ni%' <- Negate("%in%")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("filter", "dplyr")

# Functions
source('scripts/04_analysis/Supplementary_analyses/00_functions.R')




# Split by region
list_read_step4 <- split(df_filtered, df_filtered$province)




# ------------------------------------------------------------------------------- # 
#### On MOTUs ----
# ------------------------------------------------------------------------------- # 

# rank_specify
rank_choice = 'sequence'

# accumlation all plots
liste_accumulation <- lapply(list_read_step4, accumulation_curve_df, species_unit = rank_choice)

# Asymptote of all plots 
liste_asymptote <- lapply(list_read_step4, asymptote_mm, species_unit = rank_choice)

# Unlist
df_accumulation <- bind_rows(liste_accumulation, .id = "project_name")
df_asymptote <- bind_rows(liste_asymptote, .id = "project_name")

# Dataset for plot
df_accumulation_all <- left_join(df_accumulation, df_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))

# Add All samples
# Add a global saturation curve, i.e. all samples together?
all_accumulation_motu <- accumulation_curve_df(df_filtered, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)


# Asymptote of all plots 
all_asymptote_motu <- asymptote_mm(df_filtered, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, asymptote, slope)


# Bind together
df_all_accumulation <- rbind(df_accumulation, all_accumulation_motu)
df_all_asymptote <- rbind(df_asymptote, all_asymptote_motu)

# 
df_join_all <- df_all_accumulation %>%
  left_join(., df_all_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 1.05*asymptote, 
         position_asymptote_x = max(sites),
         position_slope_y = 0.30 * max(asymptote)) 

# Plot with facet
colnames(df_join_all) <- c("Province", "richness", "sd", "sites", "asymptote", "slope", "position_asymptote_y", "position_asymptote_x", "position_slope_y")
plot_acc_motus <- ggplot(df_join_all, aes(fill = Province)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.7, show.legend = F) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  facet_wrap(~Province, scales = "free") +
  scale_fill_manual(values=c("grey", "#A6CEE3","#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#B15928", "#FF7F00", "#CAB2D6","#6A3D9A","#FFD92F"))+ 
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("MOTUs") + 
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste(round(asymptote, 0), "MOTUs")), col = "black", size=3)
  
plot_acc_motus

ggsave("outputs/maps & plot sup/accumulation_curves.png", width=12, height = 8)

