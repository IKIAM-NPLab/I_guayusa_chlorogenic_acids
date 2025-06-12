# Library loading
library(gridExtra)
library(ggtext)
library(cowplot)


venn1_grob <- grobTree(venn1)
venn2_grob <- grobTree(venn2)
venn3_grob <- grobTree(venn3)
venn1_neg_grob <- grobTree(venn1_neg)
venn2_neg_grob <- grobTree(venn2_neg)
venn3_neg_grob <- grobTree(venn3_neg)




figure_union_supp_1 <- arrangeGrob(venn1_grob,
                                   venn2_grob,
                                   venn3_grob,
                                   barras1_pos,
                                   venn1_neg_grob,
                                   venn2_neg_grob,
                                   venn3_neg_grob,
                                   barras1_neg,
                                   layout_matrix = rbind(c(1, 2, 3),
                                                         c(4, 4, 4),
                                                         c(4, 4, 4),
                                                         c(5, 6, 7),
                                                         c(8, 8, 8),
                                                         c(8, 8, 8)))
# Adding label to the figures
figure_unionsuppone <- ggpubr::as_ggplot(figure_union_supp_1) +
  cowplot::draw_plot_label(
    label = LETTERS[1:8],
    x = c(0.03, 0.35, 0.67, 0.03, 0.03, 0.35, 0.67, 0.03),  # columnas relativas
    y = c(0.98, 0.98, 0.98, 0.81, 0.48, 0.48, 0.48, 0.31),  # filas relativas desde arriba
    size = 14
  )
# Exporting (*.pdf) file
ggsave(filename = "../Result/venn_results/Figuras/PDF/figure_union.pdf", plot = figure_unionsuppone,
       width = 175, height = 250, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/venn_results/Figuras/PNG/figure_union.png", plot = figure_unionsuppone,
       width = 175, height = 250, units = "mm", dpi = 300, scale = 2.5)




figure_union_supp_2 <- arrangeGrob(barras_filtradas,
                                   barras_filtradas_neg,
                                   layout_matrix = rbind(c(1,1),
                                                         c(2,2)))
# Adding label to the figures
figure_unionsuppone_2 <- ggpubr::as_ggplot(figure_union_supp_2) +
  cowplot::draw_plot_label(
    label = LETTERS[1:2],
    x = c(0.03, 0.03),  # columnas relativas
    y = c(0.98, 0.46),  # filas relativas desde arriba
    size = 14
  )
# Exporting (*.pdf) file
ggsave(filename = "../Result/venn_results/Figuras/PDF/figure_union_2.pdf", plot = figure_unionsuppone_2,
       width = 175, height = 250, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/venn_results/Figuras/PNG/figure_union_2.png", plot = figure_unionsuppone_2,
       width = 175, height = 250, units = "mm", dpi = 300, scale = 2.5)
