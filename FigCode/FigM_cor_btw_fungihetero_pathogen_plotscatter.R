####
#### R script for Ohigashi et al (2024)
#### correlation between heterogeneity of fungal community x abundance of pathogen
#### 2024.09.13 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(ggplot2); packageVersion("ggplot2")
library(ggExtra); packageVersion("ggExtra")
library(cowplot); packageVersion("cowplot")


### load data
# for within-site distance
ws_df <- read.table("05_distance_to_centroid_out/distance_to_centroids_withinsite.txt", header = T, sep = "\t")
# for across-site distance
as_df <- read.table("05_distance_to_centroid_out/distance_to_centroids_acrosssite.txt", header = T, sep = "\t")

# relative abundances of pathogen
fungilife <- read.table("02_Function_analysis_out/FungalTraits_specific_functions.txt", header = T, sep = "\t")

# environmental data
envdata <- read.csv("Data/soil_metadata.txt", header = T, sep = "\t")
envcategory <- envdata[,1:4]


### format data
# combine distance data and envcategory data
# relocate rownames to "Sample" column
ws_df <- ws_df |> rownames_to_column(var = "Sample") 
as_df <- as_df |> rownames_to_column(var = "Sample")

# merge envcategory and distance
ws.cat_df <- merge(envcategory, ws_df, by = "Sample", sort = F)
as.cat_df <- merge(envcategory, as_df, by = "Sample", sort = F)

# extract only envcategory and distance to centroids data (fungi)
ws_fungi <- ws.cat_df |> select(Sample, Site, Landuse, fungitaxa_within)
as_fungi <- as.cat_df |> select(Sample, Site, Landuse, fungitaxa_across)

# add a column for relative abundance of pathogenic fungi
ws_fungi <- ws_fungi |> mutate(Pathotroph = fungilife$Pathotroph)
as_fungi <- as_fungi |> mutate(Pathotroph = fungilife$Pathotroph)


### create plots
## within-site correlation
# obtain the result of correlation test
cor_res.w <- read.csv("07_correlation_fungtaxahetero_pathogenabund_out/cor_taxafungi_pathoabun_within.csv")
cor_res.w <- cor_res.w |>
  mutate(pval = case_when(
    p.value < 0.001 ~ "p < 0.001",
    p.value < 0.01 ~ "p < 0.01",
    p.value < 0.05 ~ "p < 0.05",
    TRUE ~ "n.s."
  ))

cor_anno.w <- "Site\n"
for (j in 1:nrow(cor_res.w)) {
  cor_anno.w0 <- paste0(cor_res.w[j, 1], ": ", "r = ", round(cor_res.w[j, 2], 2), ", ", cor_res.w[j, 4], "\n")
  cor_anno.w0 <- sub("-", "\U2212", cor_anno.w0)
  cor_anno.w <- paste0(cor_anno.w, cor_anno.w0)
}

# make a list of sites that exhibit a correlation between the pairs
sig_site <- cor_res.w |>
  filter(pval != "n.s." & X != "All") |>
  pull(X)

# plot
ws_fungi$Landuse <- factor(ws_fungi$Landuse, levels = c("Natural", "Farm"))
ws_g <- ggplot(ws_fungi, aes(x = fungitaxa_within, y = Pathotroph, color = Site, shape = Landuse)) +
  geom_point(size=3.5) + 
  scale_shape_manual(values = c(16, 3))+
  labs(x = "Distance to centroids of fungal communities",
       y = "Relative abundance of pathotroph (%)",
       color = "Site",
       shape = "Land use",
       title = "Fungal heterogeneity × Pathogenic fungi (within-site)"
  ) + 
  geom_smooth(data = subset(ws_fungi, Site %in% sig_site),
              method = "lm", se = FALSE, #inherit.aes = FALSE,
              aes(group = Site)) +
  theme_classic()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        # legend.position = "bottom", # to do marginal plot, do not set legend here
        # legend.direction = "horizontal",
        # legend.box = "vertical",
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold", colour = "black")
  ) +
  annotate("text",
           x = 0.52,
           y = max(ws_fungi$Pathotroph)*0.83,
           label = cor_anno.w,
           hjust = "left",
           size = 4
  )
# check
# ws_g


## across-site correlation
cor_res.a <- read.csv("07_correlation_fungtaxahetero_pathogenabund_out/cor_taxafungi_pathoabun_across.csv")
cor_res.a <- cor_res.a |>
  mutate(pval = case_when(
    p.value < 0.001 ~ "p < 0.001",
    p.value < 0.01 ~ "p < 0.01",
    p.value < 0.05 ~ "p < 0.05",
    TRUE ~ "n.s."
  ))

cor_anno.a <- paste0("r = ", round(cor_res.a[1, 2], 2), ", ", cor_res.a[1, 4])
cor_anno.a <- sub("-", "\U2212", cor_anno.a)

# plot
as_fungi$Landuse <- factor(as_fungi$Landuse, levels = c("Natural", "Farm"))
as_g <- ggplot(as_fungi, aes(x = fungitaxa_across, y = Pathotroph, color = Landuse, shape = Site)) +
  geom_point(size=3.5) + 
  scale_color_manual(values=c("Farm"="tan1", "Natural"="palegreen"))+
  labs(x = "Distance to centroids of fungal communities",
       y = "Relative abundance of pathotroph (%)",
       color = "Land use",
       shape = "Site",
       title = "Fungal heterogeneity × Pathogenic fungi (across-site)"
  ) + 
  theme_linedraw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold", colour = "black"),
        legend.position = "none", # to do marginal plot, do not set legend here
        # legend.direction = "horizontal",
        # legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) +
  annotate("text",
           x = max(as_fungi$fungitaxa_across),
           y = max(as_fungi$Pathotroph),
           label = cor_anno.a,
           hjust = "right",
           size = 4)

as_g_fin <- ggMarginal(as_g,
                       groupColour = TRUE,
                       groupFill = TRUE
                       )
# check
# as_g_fin

### save data
# create a directory
dir.create("FigCode/FigM_cor_btw_fungihetero_pathogen_out")

# save object and image files
# within-site
saveRDS(ws_g, file = "FigCode/FigM_cor_btw_fungihetero_pathogen_out/fungihetero_pathogen_within.obj")
ggsave(ws_g, file = "FigCode/FigM_cor_btw_fungihetero_pathogen_out/fungihetero_pathogen_within.png",
       bg = "white", width = 8, height = 7)

# across-site
saveRDS(as_g_fin, file = "FigCode/FigM_cor_btw_fungihetero_pathogen_out/fungihetero_pathogen_across.obj")
ggsave(as_g_fin, file = "FigCode/FigM_cor_btw_fungihetero_pathogen_out/fungihetero_pathogen_across.png",
       bg = "white", width = 8, height = 7)

### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigM_cor_btw_fungihetero_pathogen_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))



