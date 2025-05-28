#Negative


# Exporting plot

## General PCA

#Merge of score and loading plots of PCA analysis.



# Library loading
library(gridExtra)
# Figure matrix
figure_1 <- arrangeGrob(figure_1a,
                        figure_1b,
                        figure_1c,
                        figure_1d,
                        layout_matrix = rbind(c(1, 2, 3),
                                              c(4, 4, 4),
                                              c(4, 4, 4)))
# Adding label to the figures
figure_one <- ggpubr::as_ggplot(figure_1) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_1_glog.pdf", plot = figure_one,
      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_1_glog.png", plot = figure_one,
      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)

#PCA score plots of different factors.



# Figure matrix
figure_1_pca <- arrangeGrob(figure_1a,
                            figure_1b,
                            figure_1c,
                            layout_matrix = rbind(c(1, 2, 3)))
# Adding label to the figures
figure_one_pca <- ggpubr::as_ggplot(figure_1_pca) +
  draw_plot_label(label = LETTERS[1:3],
                  x = c(0, 0.332, 0.665),
                  y = c(.99, .99, .99))

# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_1_glog_pca.pdf",
       plot = figure_one_pca, width = 175, height = 45, units = "mm",
       dpi = 300, scale = 2.5)

# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_1_glog_pca.png",
       plot = figure_one_pca, width = 175, height = 45, units = "mm",
       dpi = 300, scale = 2.5)

## PCA by location

### Alto Pano



# Figure matrix
figure_s1 <- arrangeGrob(figure_s1a,
                         figure_s1b,
                         figure_s1c,
                         figure_s1d,
                         layout_matrix = rbind(c(1, 2, 3),
                                               c(4, 4, 4),
                                               c(4, 4, 4)))
# Adding label to the figures
figure_sone <- ggpubr::as_ggplot(figure_s1) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))

# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s1.pdf", plot = figure_sone,
      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)

# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s1.png", plot = figure_sone,
      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)

### Alto Tena



# Figure matrix
figure_s2 <- arrangeGrob(figure_s1a,
                         figure_s1b,
                         figure_s1c,
                         figure_s1e,
                         layout_matrix = rbind(c(1, 2, 3),
                                               c(4, 4, 4),
                                               c(4, 4, 4)))
# Adding label to the figures
figure_stwo <- ggpubr::as_ggplot(figure_s2) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s2.pdf", plot = figure_stwo,
      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s2.png", plot = figure_stwo,
      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)

### Talag



# Figure matrix
figure_s3 <- arrangeGrob(figure_s1a,
                         figure_s1b,
                         figure_s1c,
                         figure_s1f,
                         layout_matrix = rbind(c(1, 2, 3),
                                               c(4, 4, 4),
                                               c(4, 4, 4)))
# Adding label to the figures
figure_sthree <- ggpubr::as_ggplot(figure_s3) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s3.pdf", plot = figure_sthree,
      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s3.png", plot = figure_sthree,
     width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)

## Age plots



# Figure matrix
figure_2 <- arrangeGrob(gcms_hm,
                        vc_plot,
                        layout_matrix = rbind(c(1, 2)))
# Adding label to the figures
figure_two <- ggpubr::as_ggplot(figure_2) +
  draw_plot_label(label = LETTERS[1:2],
                  x = c(0, 0.5),
                  y = c(.99, .99))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_2.pdf", plot = figure_two,
      width = 175, height = 75, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_2.png", plot = figure_two,
      width = 175, height = 75, units = "mm", dpi = 300, scale = 2.5)

#---------------------------------------------------------------------------------

#Positive

# Exporting plot

## General PCA

#Merge of score and loading plots of PCA analysis.



# Library loading
library(gridExtra)
# Figure matrix
figure_1p <- arrangeGrob(figure_1a_p,
                        figure_1b_p,
                        figure_1c_p,
                        figure_1d_p,
                        layout_matrix = rbind(c(1, 2, 3),
                                              c(4, 4, 4),
                                              c(4, 4, 4)))
# Adding label to the figures
figure_onep <- ggpubr::as_ggplot(figure_1p) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))

# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_1_glog_p.pdf", plot = figure_onep,
      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)

# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_1_glog_p.png", plot = figure_onep,
       width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)

#PCA score plots of different factors.



# Figure matrix
figure_1p_pca <- arrangeGrob(figure_1a_p,
                            figure_1b_p,
                            figure_1c_p,
                            layout_matrix = rbind(c(1, 2, 3)))
# Adding label to the figures
figure_onep_pca <- ggpubr::as_ggplot(figure_1p_pca) +
  draw_plot_label(label = LETTERS[1:3],
                  x = c(0, 0.332, 0.665),
                  y = c(.99, .99, .99))

# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_1p_glog_pca.pdf",
       plot = figure_onep_pca, width = 175, height = 45, units = "mm",
       dpi = 300, scale = 2.5)

# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_1p_glog_pca.png",
       plot = figure_onep_pca, width = 175, height = 45, units = "mm",
       dpi = 300, scale = 2.5)


## PCA by location

### Alto Pano



# Figure matrix
figure_s1p <- arrangeGrob(figure_s1a_p,
                         figure_s1b_p,
                         figure_s1c_p,
                         figure_s1d_p,
                         layout_matrix = rbind(c(1, 2, 3),
                                               c(4, 4, 4),
                                               c(4, 4, 4)))
# Adding label to the figures
figure_sonep <- ggpubr::as_ggplot(figure_s1p) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s1p.pdf", plot = figure_sonep,
       width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s1p.png", plot = figure_sonep,
       width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)


### Alto Tena



# Figure matrix
figure_s2p <- arrangeGrob(figure_s1a_p,
                         figure_s1b_p,
                         figure_s1c_p,
                         figure_s1e_p,
                         layout_matrix = rbind(c(1, 2, 3),
                                               c(4, 4, 4),
                                               c(4, 4, 4)))
# Adding label to the figures
figure_stwop <- ggpubr::as_ggplot(figure_s2p) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s2p.pdf", plot = figure_stwop,
       width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s2p.png", plot = figure_stwop,
       width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)


### Talag



# Figure matrix
figure_s3p <- arrangeGrob(figure_s1a_p,
                         figure_s1b_p,
                         figure_s1c_p,
                         figure_s1f_p,
                         layout_matrix = rbind(c(1, 2, 3),
                                               c(4, 4, 4),
                                               c(4, 4, 4)))
# Adding label to the figures
figure_sthreep <- ggpubr::as_ggplot(figure_s3p) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s3p.pdf", plot = figure_sthreep,
       width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_s3p.png", plot = figure_sthreep,
       width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)


## Age plots



# Figure matrix
figure_2p <- arrangeGrob(gcms_hm_p,
                        vc_plot_pos,
                        layout_matrix = rbind(c(1, 2)))
# Adding label to the figures
figure_twop <- ggpubr::as_ggplot(figure_2p) +
  draw_plot_label(label = LETTERS[1:2],
                  x = c(0, 0.5),
                  y = c(.99, .99))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_2p.pdf", plot = figure_twop,
       width = 175, height = 75, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_2p.png", plot = figure_twop,
       width = 175, height = 75, units = "mm", dpi = 300, scale = 2.5)


#---------------

#Combination negative + positive
library(ggtext)

#PCA

figure_union_1 <- arrangeGrob(figure_one,
                        figure_onep,
                        layout_matrix = rbind(c(1, 1),
                                              c(2, 2)))
# Adding label to the figures
figure_uninonone <- ggpubr::as_ggplot(figure_union_1) +
  draw_plot_label(label = c("[M-H]-", "[M+H]+"),
                  x = c(0.02, 0.02),
                  y = c(.83, 0.33))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_union_1.pdf", plot = figure_uninonone,
       width = 175, height = 240, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_union_1.png", plot = figure_uninonone,
       width = 175, height = 240, units = "mm", dpi = 300, scale = 2.5)


#HCA

figure_union_2 <- arrangeGrob(figure_two,
                              figure_twop,
                              layout_matrix = rbind(c(1, 1),
                                                    c(2, 2)))
# Adding label to the figures
figure_uniontwo <- ggpubr::as_ggplot(figure_union_2) +
  draw_plot_label(label = c("[M-H]-", "[M+H]+"),
                  x = c(-0.001, 0),
                  y = c(0.997, 0.497))
# Exporting (*.pdf) file
ggsave(filename = "../Result/notame_results/Figuras/figure_union_2.pdf", plot = figure_uniontwo,
       width = 175, height = 150, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/notame_results/Figuras/figure_union_2.png", plot = figure_uniontwo,
       width = 175, height = 150, units = "mm", dpi = 300, scale = 2.5)

#---------------------------------------------------------------------------------

#Factorial plots
library(cowplot)

figure_union_3 <- arrangeGrob(plot_chakra,
                              plot_age,
                              plot_interaction_ttest,
                              layout_matrix = rbind(c(1, 2),
                                                    c(3, 3)))
# Adding label to the figures
figure_unionthree <- ggpubr::as_ggplot(figure_union_3) +
  draw_plot_label(label = LETTERS[1:3],
                  x = c(0, 0.50, 0.01),
                  y = c(1, 1, .5))
# Exporting (*.pdf) file
ggsave(filename = "../Result/Factorial_results/Figuras/PDF/figure_union_3.pdf", plot = figure_unionthree,
       width = 175, height = 150, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Result/Factorial_results/Figuras/PNG/figure_union_3.png", plot = figure_unionthree,
       width = 175, height = 150, units = "mm", dpi = 300, scale = 2.5)








