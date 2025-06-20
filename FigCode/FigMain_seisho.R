####
#### R script for Ohigashi et al (2024)
#### Combine the plots and create figures to show in the main article
#### 2025.06.13 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")

# create a directory for output files
dir.create("FigCode/FigMain_out")


### load data, create panels (combine figures), and save
## 1. PCoA of taxa and functions, and Mantel correlogram
# PCoA for prokaryotes
pcoa_p <- readRDS("04_PCoA_PERMANOVA_out/PCoA_prok_taxa.rds") 
# PCoA for fungi
pcoa_f <- readRDS("04_PCoA_PERMANOVA_out/PCoA_fungi_taxa.rds") 

# PCoA for prokaryotic functions
pcoa_pfunc <- readRDS("04_PCoA_PERMANOVA_out/PCoA_prok_functions.rds") 
# PCoA for fungi
pcoa_ffunc <- readRDS("04_PCoA_PERMANOVA_out/PCoA_fungi_functions.rds") 

# get PCoA legend
leg_pcoa <- get_legend(pcoa_p)

# Mantel correlograms for prokaryotes and fungi (taxa and function)
man_p <- readRDS("FigCode/FigM_Mantel_out/Mantel_prok_taxafunc.obj")
man_f <- readRDS("FigCode/FigM_Mantel_out/Mantel_fungi_taxafunc.obj")

# get legend for the correlograms
leg_mantel <- get_legend(man_p)

# combine figures
figM1 <- plot_grid(pcoa_p + theme(legend.position = "none",
                                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   pcoa_f + theme(legend.position = "none",
                                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   leg_pcoa,
                   pcoa_pfunc + theme(legend.position = "none",
                                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   pcoa_ffunc + theme(legend.position = "none",
                                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   NA,
                   man_p + theme(legend.position = "none",
                                 plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   man_f + theme(legend.position = "none",
                                 plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   leg_mantel,
                   ncol = 3,
                   rel_widths = c(1, 1, 0.3),
                   labels = c("a", "b", NA, "c", "d", NA, "e", "f", NA),
                   label_size = 20
)
# save PDF
cairo_pdf("FigCode/FigMain_out/FigM_PCoA_Mantel.pdf", width = 14.4, height = 18)
print(figM1)
dev.off()
# save png
ggsave(filename = "FigCode/FigMain_out/FigM_PCoA_Mantel.png",
       plot = figM1, width = 14.4, height = 18, bg = "white")


## 2. community assembly process and variation partitioning
bNTI_RC <- readRDS("FigCode/FigM_AssemblyProcess_out/assemblyprocess_ratio_plot.obj")
varpart <- readRDS("FigCode/FigM_VariationPartitioning_out/variation_explained_plot.obj")
# figM2 <- plot_grid(bNTI_RC + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
#                    varpart + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
#                    nrow = 2,
#                    rel_heights =  c(1, 1),
#                    labels = c("a", "b"),
#                    label_size = 20
# )
assem_leg <- get_legend(bNTI_RC)
varp_leg <- get_legend(varpart)
figM2 <- plot_grid(bNTI_RC + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                   legend.position = "none"),
                   assem_leg,
                   varpart + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                   legend.position = "none"),
                   varp_leg,
                   nrow = 2,
                   rel_heights = c(1, 1),
                   rel_widths = c(1, 0.25),
                   labels = c("a", NA, "b", NA),
                   label_size = 20
                   )


# save PDF
cairo_pdf("FigCode/FigMain_out/FigM_assemblyprocess.varpart.pdf", width = 15, height = 14)
print(figM2)
dev.off()
# save png
ggsave(filename = "FigCode/FigMain_out/FigM_assemblyprocess.varpart.png",
       plot = figM2, width = 15, height = 14, bg = "white")

## 3. multiple regression
# skipped (already made in FigM_multi_reg_plotcaterpillar.R)
# the figures are already made so just copy and paste
current_folder <- "FigCode/FigM_multi_reg_out"
new_folder <- "FigCode/FigMain_out"
list_of_files <- list.files(current_folder, pattern = "FigM_multireg")
file.copy(file.path(current_folder,list_of_files), new_folder, overwrite = TRUE)


## 4. correlation between fungi taxa heterogeneity and relative abundance of pathogen, SIMPER
# correlation btw fungitaxahetero x pathogenic fungi (within-site)
cor_fungi_patho_ws <- readRDS("FigCode/FigM_cor_btw_fungihetero_pathogen_out/fungihetero_pathogen_within.obj")
# correlation btw fungitaxahetero x pathogenic fungi (across-site)
cor_fungi_patho_as <- readRDS("FigCode/FigM_cor_btw_fungihetero_pathogen_out/fungihetero_pathogen_across.obj")


# get legend for the across-site correlations
leg_cor_as <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/FigM_cor_across_legend.obj")

## contribution to heterogeneity (across-site)
contrib_hetero_asac <- readRDS("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_as_ac.obj")

# combine figures
figM4ab <- plot_grid(cor_fungi_patho_ws,
                     cor_fungi_patho_as,
                     leg_cor_as,
                     ncol = 3,
                     rel_widths = c(1.3, 1, 0.3),
                     labels = c("a", "b", NA),
                     label_size = 20
)

figM4_2 <- plot_grid(figM4ab,
                     contrib_hetero_asac +
                       theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                             axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1, vjust = 1)),
                     nrow = 2,
                     rel_heights = c(1, 1),
                     labels = c(NA, "c"),
                     label_size = 20
)

# save PDF
cairo_pdf("FigCode/FigMain_out/FigM_cor.fungihetero.patho_contrib.hetero.pdf", width = 14, height = 12)
print(figM4_2)
dev.off()
# save png
ggsave(filename = "FigCode/FigMain_out/FigM_cor.fungihetero.patho_contrib.hetero.png",
       plot = figM4_2, width = 14, height = 12, bg = "white")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigMain_seisho_re_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


