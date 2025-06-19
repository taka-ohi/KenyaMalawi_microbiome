####
#### R script for Ohigashi et al (2024)
#### Combine the plots and create figures to show in the supplementary materials
#### 2025.06.13 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")

# create a directory for output files
dir.create("FigCode/FigSupple_out")


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
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_ANCOMBC.pdf", width = 14.4, height = 15)
print(figS_ancom)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_ANCOMBC.png",
       plot = figS_ancom, width = 14.4, height = 15, bg = "white")


## 2. heterogeneity (distance to centroids) of taxa
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
figS_taxa.dist <- plot_grid(dis.cen.box_p_ws + theme(legend.position = "none",
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
                   labels = c("a", "b", NA, "c", "d", NA),
                   label_size = 20
)
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_dis.cen.taxa.pdf", width = 14.4, height = 12)
print(figS_taxa.dist)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_dis.cen.taxa.png",
       plot = figS_taxa.dist, width = 14.4, height = 12, bg = "white")



## 3. heterogeneity (distance to centroids) of functions
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
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_dis.cen.func.pdf", width = 14.4, height = 12)
print(figS_func.dist)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_dis.cen.func.png",
       plot = figS_func.dist, width = 14.4, height = 12, bg = "white")



## 4. correlations between environmental and taxonomic heterogeneity
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
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_cor.taxaenvhetero.pdf", width = 14.4, height = 12)
print(figS_cor_env.taxa)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_cor.taxaenvhetero.png",
       plot = figS_cor_env.taxa, width = 14.4, height = 12, bg = "white")



## 5. contribution of fungal ASVs summarized in lifestyle
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
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_contrib.hetero_fungilife.pdf", width = 14.4, height = 18)
print(figS_contrib.fungilife)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_contrib.hetero_fungilife.png",
       plot = figS_contrib.fungilife, width = 14.4, height = 18, bg = "white")


## 6. pathotroph abundance
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
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_relabun_pathotrophs.pdf", width = 12, height = 13)
print(figS_rel_patho)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_relabun_pathotrophs.png",
       plot = figS_rel_patho, width = 12, height = 13, bg = "white")



## 7. bNTI and soil moisture
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
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_bNTI_moisture.pdf", width = 13, height = 6)
print(figS_bnti.wa)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_bNTI_moisture.png",
       plot = figS_bnti.wa, width = 13, height = 6, bg = "white")


## 8. PCoA (buffer)
# load data
b_pcoa.buffer <- readRDS("FigCode/FigS_NMDS_check.buffer.effect_out/PCoA_plot_prok.buffer.rds")
f_pcoa.buffer <- readRDS("FigCode/FigS_NMDS_check.buffer.effect_out/PCoA_plot_fungi.buffer.rds")

# get legend
leg_pcoa.buffer <- get_legend(b_pcoa.buffer)

# combine figures
figS_pcoa <- plot_grid(b_pcoa.buffer + theme(legend.position = "none",
                                             plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                       f_pcoa.buffer + theme(legend.position = "none",
                                             plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                       leg_pcoa.buffer,
                       ncol = 3,
                       rel_widths =  c(1, 1, 0.4),
                       labels = c("a", "b"),
                       label_size = 20
)
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_pcoa_buffer.pdf", width = 14.4, height = 6)
print(figS_pcoa)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_pcoa_buffer.png",
       plot = figS_pcoa, width = 14.4, height = 6, bg = "white")



## 9. Rarefaction curve
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
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_rarecurves.pdf", width = 12, height = 6)
print(figS_rare)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_rarecurves.png",
       plot = figS_rare, width = 12, height = 6, bg = "white")



## 10. results of sensitivity analysis
# load data
# PCoA plot
sens_pcoa <- readRDS("FigCode/FigS_FungalTraits_summary_out/pcoa_filtered.rds")

# correlation between distances to centroids of fungal taxa and relative abundance of pathogen
sens_patho_ws <- readRDS("FigCode/FigS_FungalTraits_summary_out/cor_taxahetero_pathogen_ws.rds")
sens_patho_as <- readRDS("FigCode/FigS_FungalTraits_summary_out/cor_taxahetero_pathogen_as.rds")

# Mantel correlogram
sens_mantel_list <- readRDS("FigCode/FigS_FungalTraits_summary_out/sensitivity_Mantel_plotlist.rds")
# extract plot from the list
fungi1 <- sens_mantel_list[["fungi1"]]

# legend
# across-site correlation
leg_cor_as <- readRDS("FigCode/FigM_correlation_btw_heterogeneity_out/FigM_cor_across_legend.obj")
# Mantel
mantel_leg <- get_legend(fungi1)

# combine plots
figS_sensi <- plot_grid(sens_pcoa+
                          theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                        fungi1 + theme(legend.position = "none",
                                       plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                        mantel_leg,
                        sens_patho_ws +
                          theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                        sens_patho_as,
                        leg_cor_as,
                        nrow = 2,
                        ncol = 3,
                        rel_widths = c(1.3, 1, 0.45),
                        labels = c("a",  "b", NA, "c", "d", NA),
                        label_size = 20
                        )


# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_sensitivity.pdf", width = 14, height = 12)
print(figS_sensi)
dev.off()
# save png
ggsave(filename = "FigCode/FigSupple_out/FigS_sensitivity.png",
       plot = figS_sensi, width = 14, height = 12, bg = "white")





### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigSupple_seisho_re_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


