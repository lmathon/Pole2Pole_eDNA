library(tidyverse)
library(rfishbase)
library(ggpubr)

europa <- read.csv("data/checklist_Europa.csv", sep=";", header = F)
colnames(europa) <- c("Family", "nb_species")

load("Rdata/03-filter-data.Rdata")

eparses <- df_filtered %>%
  filter(province=="Western_Indian_Ocean")%>%
  distinct(sequence, .keep_all=T)%>%
  select(sequence, family_name_corrected)%>%
  filter(!is.na(family_name_corrected))

eparses <- as.data.frame(table(eparses$family_name_corrected))
colnames(eparses) <- c("Family", "nb_MOTUs")


MOTU_checklist <- full_join(eparses, europa)
MOTU_checklist[is.na(MOTU_checklist)] <- 0

cor.test(MOTU_checklist$nb_MOTUs, MOTU_checklist$nb_species, method="pearson")
lm <- lm(nb_MOTUs~nb_species, data=MOTU_checklist)
summary(lm)

plot_eparses <- ggplot(MOTU_checklist, aes(x=nb_species, y=nb_MOTUs))+
  geom_point()+
  geom_abline(slope = 0.6, intercept = 0.85, size=0.8)+
  annotate(geom="text", x=45, y=3, label="Pearson\nr = 0.82\np < 0.001", hjust=1, size=3.2) +
  ggtitle("Scattered Islands")+
  xlab("Number of species per family in the checklist")+
  ylab("Number of MOTUs per family")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.1,0.1,0.2), "cm"))

###################################################################################################################################

nc <- read.csv("data/checklist-new-caledonia.csv", sep=";", header = T)
nc <- nc %>% select(Species)

x <- load_taxa()
taxo <- collect(x)

nc <- left_join(nc, taxo[,c("Species", "Family")])
nc <- nc %>% filter(!is.na(Family))
nc <- as.data.frame(table(nc$Family))

colnames(nc) <- c("Family", "nb_species")


nc_MOTUs <- df_filtered %>%
  filter(province=="Tropical_Southwestern_Pacific")%>%
  distinct(sequence, .keep_all=T)%>%
  select(sequence, family_name_corrected)%>%
  filter(!is.na(family_name_corrected))

nc_MOTUs <- as.data.frame(table(nc_MOTUs$family_name_corrected))
colnames(nc_MOTUs) <- c("Family", "nb_MOTUs")


MOTU_checklist2 <- full_join(nc_MOTUs, nc)
MOTU_checklist2[is.na(MOTU_checklist2)] <- 0


cor.test(MOTU_checklist2$nb_MOTUs, MOTU_checklist2$nb_species, method="pearson")
lm <- lm(nb_MOTUs~nb_species, data=MOTU_checklist2)
summary(lm)

plot_nc <- ggplot(MOTU_checklist2, aes(x=nb_species, y=nb_MOTUs))+
  geom_point()+
  geom_abline(slope = 0.45, intercept = -0.92, size=0.8)+
  annotate(geom="text", x=200, y=10, label="Pearson\nr = 0.91\np < 0.001", hjust=1, size=3.2) +
  ggtitle("New-Caledonia")+
  xlab("Number of species per family in the checklist")+
  ylab("Number of MOTUs per family")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.1,0.1,0.2), "cm"))


##################################################################################################################################

lengguru <- read.csv("data/checklist_lengguru.csv", sep=";", header = T)

lengguru <- lengguru %>% select("Species", "Family")
lengguru <- as.data.frame(table(lengguru$Family))

colnames(lengguru) <- c("Family", "nb_species")


leng_MOTUs <- df_filtered %>%
  filter(province=="Western_Coral_Triangle")%>%
  distinct(sequence, .keep_all=T)%>%
  select(sequence, family_name_corrected)%>%
  filter(!is.na(family_name_corrected))

leng_MOTUs <- as.data.frame(table(leng_MOTUs$family_name_corrected))
colnames(leng_MOTUs) <- c("Family", "nb_MOTUs")


MOTU_checklist3 <- full_join(leng_MOTUs, lengguru)
MOTU_checklist3[is.na(MOTU_checklist3)] <- 0


cor.test(MOTU_checklist3$nb_MOTUs, MOTU_checklist3$nb_species, method="pearson")
lm <- lm(nb_MOTUs~nb_species, data=MOTU_checklist3)
summary(lm)

plot_leng <- ggplot(MOTU_checklist3, aes(x=nb_species, y=nb_MOTUs))+
  geom_point()+
  geom_abline(slope = 0.3, intercept = 0.8, size=0.8)+
  annotate(geom="text", x=350, y=10, label="Pearson\nr = 0.84\np < 0.001", hjust=1, size=3.2) +
  ggtitle("Lengguru")+
  xlab("Number of species per family in the checklist")+
  ylab("Number of MOTUs per family")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.1,0.1,0.2), "cm"))


############################################################################################

med <- read.csv("data/checklist_mediterranee.csv", sep=";", header = T)

med <- med %>% select("Species", "Family")
med <- as.data.frame(table(med$Family))

colnames(med) <- c("Family", "nb_species")


med_MOTUs <- df_filtered %>%
  filter(province=="Mediterranean_Sea")%>%
  distinct(sequence, .keep_all=T)%>%
  select(sequence, family_name_corrected)%>%
  filter(!is.na(family_name_corrected))

med_MOTUs <- as.data.frame(table(med_MOTUs$family_name_corrected))
colnames(med_MOTUs) <- c("Family", "nb_MOTUs")


MOTU_checklist4 <- full_join(med_MOTUs, med)
MOTU_checklist4[is.na(MOTU_checklist4)] <- 0


cor.test(MOTU_checklist4$nb_MOTUs, MOTU_checklist4$nb_species, method="pearson")
lm <- lm(nb_MOTUs~nb_species, data=MOTU_checklist4)
summary(lm)

plot_med <- ggplot(MOTU_checklist4, aes(x=nb_species, y=nb_MOTUs))+
  geom_point()+
  ylim(0,30)+
  xlim(0,30)+
  geom_abline(slope = 1.3, intercept = 0.55, size=0.8)+
  annotate(geom="text", x=30, y=3, label="Pearson\nr = 0.81\np < 0.001", hjust=1, size=3.2) +
  ggtitle("Mediterrannean Sea")+
  xlab("Number of species per family in the checklist")+
  ylab("Number of MOTUs per family")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.1,0.1,0.2), "cm"))



all <- ggarrange(plot_eparses, plot_med, plot_nc, plot_leng, nrow=2, ncol=2)

ggsave(all, file="outputs/maps & plot sup/comparison_checklists.png")
