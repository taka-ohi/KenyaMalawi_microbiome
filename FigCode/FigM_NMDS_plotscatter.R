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
                      labels = c("a", "b", NA),
                      label_size = 16
                      )

# create a directory
dir.create("FigCode/FigM_NMDS_out")

# save PDF
cairo_pdf("FigCode/FigM_NMDS_out/FigM_NMDS.pdf", width = 14.4, height = 6)
print(fig1_all)
dev.off()

# save png
ggsave(filename = "FigCode/FigM_NMDS_out/FigM_NMDS.png",
       plot = fig1_all, width = 14.4, height = 6, bg = "white")


### save session info
# create a directory
dir.create("FigCode/Fig_SessionInfo")
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigM_NMDS_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

