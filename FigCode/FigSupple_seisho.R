####
#### R script for Ohigashi et al (2024)
#### Combine the plots and create figures to show in the supplementary materials
#### 2024.12.01 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")


### load data and create panels (combine figures)
## 1. ANCOMBC results
# load data
# Class
b_ancom_class <- readRDS("12_ANCOMBC_out/ANCOM_landuse_prokclass_plot.obj") # prokaryotes
f_ancom_class <- readRDS("12_ANCOMBC_out/ANCOM_landuse_fungiclass_plot.obj") # fungi
# Genus
b_ancom_genus <- readRDS("12_ANCOMBC_out/ANCOM_landuse_prokgenus_plot.obj") # prokaryotes
f_ancom_genus <- readRDS("12_ANCOMBC_out/ANCOM_landuse_fungigenus_plot.obj") # fungi

# get legend
leg_ancom <- get_legend(b_ancom_class) 

# combine figures
figS_ancom <- plot_grid(b_ancom_class +
                          theme(legend.position = "none", plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
                          labs(title = NULL),
                        f_ancom_class + theme(legend.position = "none", plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
                          labs(title = NULL),
                        leg_ancom,
                        b_ancom_genus + 
                          theme(legend.position = "none", plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
                          labs(title = NULL),
                        f_ancom_genus + theme(legend.position = "none", plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
                          labs(title = NULL),
                        ncol = 3,
                        rel_widths = c(1, 1, 0.3),
                        rel_heights = c(1, 1.8),
                        labels = c("a", "b", NA, "c", "d"),
                        label_size = 20
)

## 2. heterogeneity of functions
# load data
prokfunc_ws <- readRDS("FigCode/FigM_dist.to.cent_out/prokallfunc_within_dist.to.cent.obj") # prok function within site
fungifunc_ws <- readRDS("FigCode/FigM_dist.to.cent_out/fungilife_within_dist.to.cent.obj") # fungi function within site

prokfunc_as <- readRDS("FigCode/FigM_dist.to.cent_out/prokallfunc_across_dist.to.cent.obj") # prok function across site
fungifunc_as <- readRDS("FigCode/FigM_dist.to.cent_out/fungilife_across_dist.to.cent.obj") # fungi function across site

# get legend
leg_dist.cent <- get_legend(prokfunc_ws) 

# combine figures
figS_func.dist <- plot_grid(prokfunc_ws + theme(legend.position = "none",
                                                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            fungifunc_ws + theme(legend.position = "none",
                                            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            leg_dist.cent,
                            prokfunc_as + theme(legend.position = "none",
                                            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            fungifunc_as + theme(legend.position = "none",
                                            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            ncol = 3,
                            rel_widths = c(1, 1, 0.4),
                            labels = c("a", "b", NA, "c", "d"),
                            label_size = 20
)


## 3. correlations between environmental and taxonomic heterogeneity
# load data
cor_env.proktaxa_ws <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/cor_env_taxaprok_within.obj")
cor_env.fungitaxa_ws <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/cor_env_taxafungi_within.obj")

cor_env.proktaxa_as <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/cor_env_taxaprok_across.obj")
cor_env.fungitaxa_as <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/cor_env_taxafungi_across.obj")

# get legends
leg_cor_ws <- get_legend(cor_env.proktaxa_ws)
leg_cor_as <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/FigM_cor_across_legend.obj")

# combine figures
figS_cor_env.taxa <- plot_grid(cor_env.proktaxa_ws + theme(legend.position = "none",
                                                           plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                               cor_env.fungitaxa_ws + theme(legend.position = "none",
                                                            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                               leg_cor_ws,
                               cor_env.proktaxa_as,
                               cor_env.fungitaxa_as,
                               leg_cor_as,
                               ncol = 3,
                               rel_widths = c(1, 1, 0.4),
                               labels = c("a", "b", NA, "c", "d"),
                               label_size = 20
)


## 4. contribution of fungal ASVs summarized in lifestyle
# laod data
ct_siteD <- readRDS("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_within_D.obj")
ct_siteE <- readRDS("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_within_E.obj")
ct_siteF <- readRDS("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_within_F.obj")
ct_siteG <- readRDS("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_within_G.obj")
ct_siteH <- readRDS("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_within_H.obj")

ct_Kenya <- readRDS("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_within_Kenya.obj")
ct_Malawi <- readRDS("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_within_Malawi.obj")

# get legend
leg_ct <- get_legend(ct_siteD)

# combine figures
figS_contrib.fungilife <- plot_grid(ct_siteD + theme(legend.position = "none",
                                                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                                    ct_siteE + theme(legend.position = "none",
                                                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                                    ct_siteF + theme(legend.position = "none",
                                                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                                    ct_siteG + theme(legend.position = "none",
                                                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                                    ct_siteH + theme(legend.position = "none",
                                                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                                    leg_ct,
                                    ct_Kenya + theme(legend.position = "none",
                                                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                                    ct_Malawi + theme(legend.position = "none",
                                                      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                                    ncol = 2,
                                    rel_widths = c(1, 1),
                                    labels = c("a", "b", "c", "d", "e", NA, "f", "g"),
                                    label_size = 20
)


## 5. pathotroph abundance
# load data
rel_patho_lifestyle <- readRDS("FigCode/FigS_pathotroph_subanalysis_out/relabun_lifestyle.obj")
rel_patho_genus <- readRDS("FigCode/FigS_pathotroph_subanalysis_out/relabun_genus.obj")

# combine figures
figS_rel_patho <- plot_grid(rel_patho_lifestyle + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            rel_patho_genus + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            nrow = 2,
                            rel_heights =  c(1, 1),
                            labels = c("a", "b"),
                            label_size = 20
)


## 6. bNTI and soil moisture
# load data
bnti.wa.plots <- readRDS("FigCode/FigS_bNTI_cor_out/plot_list.pbj")
# combine figures
figS_bnti.wa <- plot_grid(bnti.wa.plots[["mean_prok"]] + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                          bnti.wa.plots[["mean_fungi"]] + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                          ncol = 2,
                          rel_heights =  c(1, 1),
                          labels = c("a", "b"),
                          label_size = 20
                          )




## 7. NMDS (buffer)
# load data
b_nmds.buffer <- readRDS("FigCode/FigS_NMDS_check.buffer.effect_out/NMDS_plot_prok.buffer.obj")
f_nmds.buffer <- readRDS("FigCode/FigS_NMDS_check.buffer.effect_out/NMDS_plot_fungi.buffer.obj")

# get legend
leg_nmds.buffer <- get_legend(b_nmds.buffer)

# combine figures
figS_nmds <- plot_grid(b_nmds.buffer + theme(legend.position = "none",
                                             plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                       f_nmds.buffer + theme(legend.position = "none",
                                             plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                       leg_nmds.buffer,
                       ncol = 3,
                       rel_widths =  c(1, 1, 0.4),
                       labels = c("a", "b"),
                       label_size = 20
)


## 8. Rarefaction curve
# load data
rare_prok <- readRDS("01_DADA2_out/rarefaction_curve_prok.rds")
rare_fungi <- readRDS("01_DADA2_out/rarefaction_curve_fungi.rds")

# combine data
figS_rare <- plot_grid(rare_prok + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                       rare_fungi + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                       ncol = 2,
                       rel_widths =  c(1, 1),
                       labels = c("a", "b"),
                       label_size = 20
)


########### save data ###########
# create a directory for output files
dir.create("FigCode/FigSupple_out")

## 1. ANCOMBC
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_ANCOMBC.pdf", width = 14.4, height = 15)
print(figS_ancom)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_ANCOMBC.png",
       plot = figS_ancom, width = 14.4, height = 15, bg = "white")


## 2. heterogeneity of functions
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_dis.cen.func.pdf", width = 14.4, height = 12)
print(figS_func.dist)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_dis.cen.func.png",
       plot = figS_func.dist, width = 14.4, height = 12, bg = "white")


## 3. correlations between environmental and taxonomic heterogeneity
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_cor.taxaenvhetero.pdf", width = 14.4, height = 12)
print(figS_cor_env.taxa)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_cor.taxaenvhetero.png",
       plot = figS_cor_env.taxa, width = 14.4, height = 12, bg = "white")


## 4. contribution of fungal ASVs summarized in lifestyle
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_contrib.hetero_fungilife.pdf", width = 14.4, height = 18)
print(figS_contrib.fungilife)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_contrib.hetero_fungilife.png",
       plot = figS_contrib.fungilife, width = 14.4, height = 18, bg = "white")


## 5. pathotroph abundance
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_relabun_pathotrophs.pdf", width = 12, height = 13)
print(figS_rel_patho)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_relabun_pathotrophs.png",
       plot = figS_rel_patho, width = 12, height = 13, bg = "white")


## 6. bNTI and moisture
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_bNTI_moisture.pdf", width = 13, height = 6)
print(figS_bnti.wa)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_bNTI_moisture.png",
       plot = figS_bnti.wa, width = 13, height = 6, bg = "white")


## 7. NMDS (buffer)
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_nmds_buffer.pdf", width = 14.4, height = 6)
print(figS_nmds)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_nmds_buffer.png",
       plot = figS_nmds, width = 14.4, height = 6, bg = "white")


## 8. Rarefaction curve
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_rarecurves.pdf", width = 12, height = 6)
print(figS_rare)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_rarecurves.png",
       plot = figS_rare, width = 12, height = 6, bg = "white")



### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigSupple_seisho_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


