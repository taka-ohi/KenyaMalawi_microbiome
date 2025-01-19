####
#### R script for Ohigashi et al (2024)
#### Combine the plots and create figures to show in the main article
#### 2024.09.17 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")


### load data and create panels (combine figures)
## 1. NMDS and distance to centroids (within- and across-site)
# NMDS for prokaryotes
nmds_p <- readRDS("04_NMDS_PERMANOVA_out/NMDS_plot_prok.obj") 
# NMDS for fungi
nmds_f <- readRDS("04_NMDS_PERMANOVA_out/NMDS_plot_fungi.obj") 
# get NMDS legend
leg_nmds <- get_legend(nmds_p) 

# distance to centroids for prokaryotes (within-site)
dis.cen.box_p_ws <- readRDS("FigCode/FigM_dist.to.cent_out/proktaxa_within_dist.to.cent.obj") 
# distance to centroids for fungi (within-site)
dis.cen.box_f_ws <- readRDS("FigCode/FigM_dist.to.cent_out/fungitaxa_within_dist.to.cent.obj") 
# get legend for the boxplot
leg_dis.cen.box <- get_legend(dis.cen.box_p_ws)

# distance to centroids for prokaryotes (across-site)
dis.cen.box_p_as <- readRDS("FigCode/FigM_dist.to.cent_out/proktaxa_across_dist.to.cent.obj") 
# distance to centroids for fungi (across-site)
dis.cen.box_f_as <- readRDS("FigCode/FigM_dist.to.cent_out/fungitaxa_across_dist.to.cent.obj") 

# combine figures
figM1 <- plot_grid(nmds_p + theme(legend.position = "none",
                                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   nmds_f + theme(legend.position = "none",
                                      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   leg_nmds,
                   dis.cen.box_p_ws + theme(legend.position = "none",
                                            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   dis.cen.box_f_ws + theme(legend.position = "none",
                                            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   leg_dis.cen.box,
                   dis.cen.box_p_as + theme(legend.position = "none",
                                            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   dis.cen.box_f_as + theme(legend.position = "none",
                                            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   ncol = 3,
                   rel_widths = c(1, 1, 0.4),
                   labels = c("a", "b", NA, "c", "d", NA, "e", "f"),
                   label_size = 20
)


## 2. correlations between taxa heterogeneity and functional heterogeneity, and Mantel correlograms
# correlation btw proktaxahetero x prokaryotic all functions (within-site)
cor_prok_funcall_ws <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/cor_taxaprok_funcall_within.obj")
# correlation btw fungitaxahetero x fungal lifestyles (within-site)
cor_fungi_funcall_ws <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/cor_taxafungi_fungilife_within.obj")
# get legend for the within-site correlations
leg_cor_ws <- get_legend(cor_prok_funcall_ws)

# correlation btw proktaxahetero x prokaryotic all functions (across-site)
cor_prok_funcall_as <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/cor_taxaprok_funcall_across.obj")
# correlation btw fungitaxahetero x fungal lifestyles (across-site)
cor_fungi_funcall_as <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/cor_taxafungi_fungilife_across.obj")
# get legend for the across-site correlations
leg_cor_as <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/FigM_cor_across_legend.obj")

# Mantel correlograms for prokaryotes and fungi (taxa and function)
man_p <- readRDS("FigCode/FigM_Mantel_out/Mantel_prok_taxafunc.obj")
man_f <- readRDS("FigCode/FigM_Mantel_out/Mantel_fungi_taxafunc.obj")
# get legend for the correlograms
leg_mantel <- get_legend(man_p)

# combine figures
figM2 <- plot_grid(cor_prok_funcall_ws + theme(legend.position = "none",
                                               plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   cor_fungi_funcall_ws + theme(legend.position = "none",
                                                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   leg_cor_ws,
                   cor_prok_funcall_as,
                   cor_fungi_funcall_as,
                   leg_cor_as,
                   man_p + theme(legend.position = "none",
                                 plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   man_f + theme(legend.position = "none",
                                 plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   leg_mantel,
                   ncol = 3,
                   rel_widths = c(1, 1, 0.4),
                   labels = c("a", "b", NA, "c", "d", NA, "e", "f", NA),
                   label_size = 20
)


## 3. partitioning the dissimilarity
# skipped (already made in FigM_part_dissimilarity_plottriangle.R)

## 3'. community assembly process and variation partitioning
bNTI_RC <- readRDS("FigCode/FigM_AssemblyProcess_out/assemblyprocess_ratio_plot.obj")
varpart <- readRDS("FigCode/FigM_VariationPartitioning_out/variation_explained_plot.obj")
figM3 <- plot_grid(bNTI_RC + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   varpart + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   nrow = 2,
                   rel_heights =  c(1, 1),
                   labels = c("a", "b"),
                   label_size = 20
)

## 4. multiple regression
# skipped (already made in FigM_multi_reg_plotcaterpillar.R)


## 5. correlation between fungi taxa heterogeneity and relative abundance of pathogen, SIMPER
# correlation btw fungitaxahetero x pathogenic fungi (within-site)
cor_fungi_patho_ws <- readRDS("FigCode/FigM_cor_btw_fungihetero_pathogen_out/fungihetero_pathogen_within.obj")
# correlation btw fungitaxahetero x pathogenic fungi (across-site)
cor_fungi_patho_as <- readRDS("FigCode/FigM_cor_btw_fungihetero_pathogen_out/fungihetero_pathogen_across.obj")

## SIMPER analyzed only considering land uses
# simper_allsite <- readRDS("FigCode/FigM_SIMPER_out/SIMPER_allsites.obj")
# 
# # combine figures
# figM5 <- plot_grid(cor_fungi_patho_ws,
#                    cor_fungi_patho_as,
#                    leg_cor_as,
#                    simper_allsite,
#                    ncol = 3,
#                    rel_widths = c(1.3, 1, 0.3),
#                    labels = c("a", "b", NA, "c"),
#                    label_size = 20
# )

## contribution to heterogeneity (across-site)
contrib_hetero_asac <- readRDS("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_as_ac.obj")

figM5ab <- plot_grid(cor_fungi_patho_ws,
                     cor_fungi_patho_as,
                     leg_cor_as,
                     ncol = 3,
                     rel_widths = c(1.3, 1, 0.3),
                     labels = c("a", "b", NA),
                     label_size = 20
)

figM5_2 <- plot_grid(figM5ab,
                     contrib_hetero_asac +
                       theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                             axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1, vjust = 1)),
                     nrow = 2,
                     rel_heights = c(1, 1),
                     labels = c(NA, "c"),
                     label_size = 20
                     )



### save data
# create a directory for output files
dir.create("FigCode/FigMain_out")

## 1. NMDS and distance to centroids (within- and across-site)
# save PDF
cairo_pdf("FigCode/FigMain_out/FigM_NMDS_dis.cen.pdf", width = 14.4, height = 18)
print(figM1)
dev.off()
# save png
ggsave(filename = "FigCode/FigMain_out/FigM_NMDS_dis.cen.png",
       plot = figM1, width = 14.4, height = 18, bg = "white")


## 2. correlations between taxa heterogeneity and functional heterogeneity, and Mantel correlograms
# save PDF
cairo_pdf("FigCode/FigMain_out/FigM_cor.taxafunchetero_Mantel.pdf", width = 14.4, height = 18)
print(figM2)
dev.off()
# save png
ggsave(filename = "FigCode/FigMain_out/FigM_cor.taxafunchetero_Mantel.png",
       plot = figM2, width = 14.4, height = 18, bg = "white")


## 3. partitioning the dissimilarity
# the figures are already made so just copy and paste
# current_folder <- "FigCode/FigM_part_dissimilarity_out"
# new_folder <- "FigCode/FigMain_out"
# list_of_files <- list.files(current_folder, pattern = "part_dissimilarity")
# file.copy(file.path(current_folder,list_of_files), new_folder, overwrite = TRUE)

## 3'. community assembly process and variation partitioning
# save PDF
cairo_pdf("FigCode/FigMain_out/FigM_assemblyprocess.varpart.pdf", width = 12, height = 13)
print(figM3)
dev.off()
# save png
ggsave(filename = "FigCode/FigMain_out/FigM_assemblyprocess.varpart.png",
       plot = figM3, width = 12, height = 13, bg = "white")



## 4. multiple regression
# the figures are already made so just copy and paste
current_folder2 <- "FigCode/FigM_multi_reg_out"
new_folder <- "FigCode/FigMain_out"
list_of_files2 <- list.files(current_folder2, pattern = "FigM_multireg")
file.copy(file.path(current_folder2,list_of_files2), new_folder, overwrite = TRUE)


## 5. correlation between fungi taxa heterogeneity and relative abundance of pathogen, SIMPER
# save PDF
# cairo_pdf("FigCode/FigMain_out/FigM_cor.fungihetero.patho_SIMPER.pdf", width = 14, height = 12)
cairo_pdf("FigCode/FigMain_out/FigM_cor.fungihetero.patho_contrib.hetero.pdf", width = 14, height = 12)
print(figM5_2)
dev.off()
# save png
# ggsave(filename = "FigCode/FigMain_out/FigM_cor.fungihetero.patho_SIMPER.png",
#        plot = figM5, width = 14, height = 12, bg = "white")
ggsave(filename = "FigCode/FigMain_out/FigM_cor.fungihetero.patho_contrib.hetero.png",
       plot = figM5_2, width = 14, height = 12, bg = "white")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigMain_seisho_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))



