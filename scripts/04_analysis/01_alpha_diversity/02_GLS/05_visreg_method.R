library(ggplot2)
library(ggpubr)

load("Rdata/fit.method.FDq0.rdata")
load("Rdata/fit.method.motus.rdata")
load("Rdata/fit.method.crypto.rdata")
load("Rdata/fit.method.large.rdata")


fit.method.motus <- fit.method.motus+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

fit.method.crypto <- fit.method.crypto+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

fit.method.large <- fit.method.large+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

fit.method.FDq0 <- fit.method.FDq0+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))


ggarrange(fit.method.motus, fit.method.crypto, fit.method.large, fit.method.FDq0,
          nrow=2, ncol=2, labels = c("a", "b", "c", "d"))

ggsave("outputs/GLS/visreg_method.png", width=7, height = 7)
