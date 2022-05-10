library(ggpubr)

load("Rdata/fit.grav_med.crypto.rdata")
load("Rdata/fit.grav_med.FDq0.rdata")
load("Rdata/fit.grav_med.motus.rdata")
load("Rdata/fit.grav_med.large.rdata")

fit.grav_med.motus <- fit.grav_med.motus +
  ylab("Marine Ecosystem Dependence")+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2,'cm'))

fit.grav_med.crypto <- fit.grav_med.crypto +
  ylab("Marine Ecosystem Dependence")+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2,'cm'))

fit.grav_med.large <- fit.grav_med.large +
  ylab("Marine Ecosystem Dependence")+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2,'cm'))

fit.grav_med.FDq0 <- fit.grav_med.FDq0 +
  ylab("Marine Ecosystem Dependence")+
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2,'cm'))



ggarrange(fit.grav_med.motus, fit.grav_med.crypto, fit.grav_med.large, fit.grav_med.FDq0, labels = c("A", "B", "C", "D"))


ggsave("outputs/Figures_papier/Fig4.pdf", width = 7, height = 8)
ggsave("outputs/Figures_papier/Fig4.png", width = 7, height = 8)
