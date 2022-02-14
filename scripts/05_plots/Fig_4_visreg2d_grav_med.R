library(ggpubr)

load("Rdata/fit.grav_med.crypto.rdata")
load("Rdata/fit.grav_med.FDq0.rdata")
load("Rdata/fit.grav_med.motus.rdata")
load("Rdata/fit.grav_med.large.rdata")

ggarrange(fit.grav_med.motus, fit.grav_med.crypto, fit.grav_med.large, fit.grav_med.FDq0, labels = c("A", "B", "C", "D"))


ggsave("outputs/Figures_papier/Fig4.pdf", width = 7, height = 6)
