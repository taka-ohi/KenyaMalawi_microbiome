####
#### R script for Ohigashi et al (2024)
#### Figure for NMDS plot
#### 2024.07.22 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data (created in 04_NMDS_PERMANOVA.R)
# plot data
nmds_prok <- readRDS("04_NMDS_PERMANOVA_out/NMDS_plot_prok.obj")
nmds_fungi <- readRDS("04_NMDS_PERMANOVA_out/NMDS_plot_fungi.obj")

# extract legend data
nmds_legend <- get_legend(nmds_prok) # just used prokaryotes one (fungi one is also fine)


### combine and save figures
fig1_all <- plot_grid(nmds_prok + theme(legend.position = "none",
                                        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                      nmds_fungi + theme(legend.position = "none",
                                         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                      nmds_legend,
                      ncol = 3,
                      rel_widths = c(1, 1, 0.4),
                      labels = c("a", "b", NA)
                      )
dir.create("FigCode/Fig1_out")
# ggsave(filename = "FigCode/Fig1_out/Fig1_NMDS.pdf", plot = fig1_all, width = 13.2, height = 5.5,
#        device = cairo_pdf()
#        )
# somehow generates garbled texts... so give up pdf here

ggsave(filename = "FigCode/Fig1_out/Fig1_NMDS.png", plot = fig1_all, width = 13.2, height = 5.5)

