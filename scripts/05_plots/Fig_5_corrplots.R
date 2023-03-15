library(ggpubr)
library(gridExtra)

load("Rdata/plot_alpha_phylo_gen.rdata")
load("Rdata/plot_alpha_trait_gen.rdata")
load("Rdata/plot_beta_phylo_gen.rdata")
load("Rdata/plot_beta_trait_gen.rdata")


remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}

plot_alpha_phylo_gen <- remove_geom(plot_alpha_phylo_gen, "GeomText")
plot_alpha_phylo_gen <- plot_alpha_phylo_gen +
  coord_flip()+
  annotate(geom="text", x=4, y=20, label="Pearson cor=0.94 \n p<0.001", hjust=0, size=6, fontface = "bold")


plot_alpha_trait_gen <- remove_geom(plot_alpha_trait_gen, "GeomText")
plot_alpha_trait_gen <- plot_alpha_trait_gen +
  ylim(0,22)+
  coord_flip()+
  annotate(geom="text", x=4, y=12, label="Pearson cor=0.91 \n p<0.001", hjust=0, size=6, fontface = "bold")


plot_beta_phylo_gen <- remove_geom(plot_beta_phylo_gen, "GeomText")
plot_beta_phylo_gen <- plot_beta_phylo_gen +
  coord_flip()+
  annotate(geom="text", x=0, y=0.6, label="Mantel=0.91", hjust=0, size=6, fontface = "bold")

plot_beta_trait_gen <- remove_geom(plot_beta_trait_gen, "GeomText")
plot_beta_trait_gen <- plot_beta_trait_gen +
  coord_flip()+
  annotate(geom="text", x=0, y=0.6, label="Mantel=0.85", hjust=0, size=6, fontface = "bold")


ggarrange(plot_alpha_phylo_gen, plot_alpha_trait_gen, plot_beta_phylo_gen, plot_beta_trait_gen,
          nrow=2, ncol=2, labels = c("A", "B", "C", "D"), font.label = list(size = 24))

ggsave("outputs/Figures_papier/Fig5.pdf", width = 13, height = 14, dpi = 600)
ggsave("outputs/Figures_papier/Fig5.png", width = 13, height = 14, dpi = 600)
