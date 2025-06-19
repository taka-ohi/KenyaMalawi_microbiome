####
#### R script for Ohigashi et al (2024)
#### check relationship between FungalTraits annotation rates and main results
#### 2025.06.11 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
source("Function/F2_HelperFunctions_for_Visualization.R")
source("Function/F1_HelperFunctions.R")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(tibble); packageVersion("tibble")
library(ggExtra); packageVersion("ggExtra")
library(purrr); packageVersion("purrr")
library(SoDA); packageVersion("SoDA")


### load data
# distance to centroids for each sample within treatments
dists_within <- read.table("05_distance_to_centroid_out/distance_to_centroids_withinsite.txt", header = T, sep = "\t")
dists_across <- read.table("05_distance_to_centroid_out/distance_to_centroids_acrosssite.txt", header = T, sep = "\t")

# Fungal traits annotation rate
anno_rates <- read.csv("FigCode/FigS_FungalTraits_summary_out/assignedASV_percent.csv", header = T)

# environment data
env_data <- read.table("Data/soil_metadata.txt", header = T)
env_category <- env_data[,1:4] # extract category data
landuse <- env_data$Landuse
site <- env_data$Site
# with "plots"
env_data2 <- env_data |>
  mutate(Plot=rep(rep(c("1", "2", "3"), each = 3), times = 10)) # add Plot column for later use

# gps
gpsdata <- read.csv("Data/soil_gpsdata.csv", header = T)

# relative abundance of pathotroph data
flife_aggregated <- read.table("02_Function_analysis_out/FungalTraits_specific_functions.txt", header = T, sep = "\t")

# fungal function data
fungilife <- read.table("02_Function_analysis_out/FungalTraits_primarilifestyle_percent.txt", header = T)
fungilife <- fungilife |> filter(primary_lifestyle != "unassigned") # 2024.09.11
fungilife <- fungilife |> column_to_rownames(var = "primary_lifestyle")
fungilife.t <- t(fungilife)


######### sensitivity analysis #########
# set threshold value to remove bottom 25% annotation rate samples
annorate_thres <- quantile(anno_rates$assigned_percent, 0.25) # 41.7%

### 1. distance to centroids x pathogen abundance
## format data
# within-site data frame
ws_df <- env_category %>%
  mutate(treatment = paste(Site, Landuse, sep = "_")) %>%
  left_join(anno_rates, by = "Sample") %>%
  left_join(dists_within %>%
              rownames_to_column("Sample") %>%
              select(Sample, fungitaxa_within, fungilife_within),
            by = "Sample")

# across-site data frame
as_df <- env_category %>%
  mutate(treatment = Landuse) %>%
  left_join(anno_rates, by = "Sample") %>%
  left_join(dists_across %>%
              rownames_to_column("Sample") %>%
              select(Sample, fungitaxa_across, fungilife_across),
            by = "Sample")

# remove samples with annotation rate < 41.7%
ws_df_filterd <- ws_df %>%
  filter(assigned_percent >= annorate_thres)
as_df_filterd <- as_df %>%
  filter(assigned_percent >= annorate_thres)


## 1-1. correlation between distances to centroids of taxa and relative abundances of pathotrophs (within-site)
ws_filterd_patho <- ws_df_filterd %>%
  left_join(flife_aggregated %>% select(Sample, Pathotroph), by = "Sample")

corr_ws_patho <- CorrInCat(ws_filterd_patho, var1 = "fungitaxa_within", var2 = "Pathotroph",
                           method = "pearson", category = "Site")
corr_ws_patho_is <- cor.test(ws_filterd_patho$fungitaxa_within, ws_filterd_patho$Pathotroph) # ignore site
corr_ws_patho_is_summary <- cbind(corr_ws_patho_is$estimate, corr_ws_patho_is$p.value)
colnames(corr_ws_patho_is_summary) <- c("r", "p.value")
row.names(corr_ws_patho_is_summary) <- "All"
corr_ws_patho <- rbind(corr_ws_patho, corr_ws_patho_is_summary) # combine
corr_ws_patho <- corr_ws_patho %>%
  rownames_to_column("Site")

# plot
sig_site <- corr_ws_patho %>% filter(p.value < 0.05, Site != "All") %>% pull(Site)
cor_res <- corr_ws_patho |>
  mutate(pval = case_when(
    p.value < 0.001 ~ "p < 0.001",
    p.value < 0.01 ~ "p < 0.01",
    p.value < 0.05 ~ "p < 0.05",
    TRUE ~ "n.s."
  ))

cor_anno <- "Site\n"
for (j in 1:nrow(cor_res)) {
  cor_anno0 <- paste0(cor_res[j, 1], ": ", "r = ", round(cor_res[j, 2], 2), ", ", cor_res[j, 4], "\n")
  cor_anno0 <- sub("-", "\U2212", cor_anno0)
  cor_anno <- paste0(cor_anno, cor_anno0)
}

ws_filterd_patho$Landuse <- factor(ws_filterd_patho$Landuse, levels = c("Natural", "Farm"))

p_patho_ws <- ggplot(ws_filterd_patho,
                     aes(x = fungitaxa_within, y = Pathotroph, color = Site, shape = Landuse)) + 
  geom_point(size=3.5) + 
  scale_shape_manual(values=c("Farm"=4, "Natural"=16)) +
  labs(x = "Distance to centroids of fungal communities",
       y = "Relative abundance of pathotroph (%)",
       # shape = "Site", color = "Land use",
       shape = "Land use", color = "Site",
       title = "Fungal heterogeneity × Pathogenic fungi\n(within-site; w/o annotation rate < 41.7%)"
  ) + 
  geom_smooth(data = subset(ws_filterd_patho, Site %in% sig_site),
              method = "lm", se = FALSE, #inherit.aes = FALSE,
              aes(group = Site)) +
  theme_classic()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 13, face = "bold", colour = "black"),
        legend.text = element_text(size = 13, face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", colour = "black")
  ) +
  annotate("text",
           x = 0.52,
           y = max(ws_filterd_patho$Pathotroph)*0.83,
           label = cor_anno,
           hjust = "left",
           size = 4
  )


## 1-2. correlation between distances to centroids of taxa and relative abundances of pathotrophs (across-site)
as_filterd_patho <- as_df_filterd %>%
  left_join(flife_aggregated %>% select(Sample, Pathotroph), by = "Sample")

corr_as_test_patho <- cor.test(as_filterd_patho$fungitaxa_across, as_filterd_patho$Pathotroph)
corr_as_patho <- cbind(corr_as_test_patho$estimate, corr_as_test_patho$p.value)
colnames(corr_as_patho) <- c("r", "p.value")


# plot
cor_res <- corr_as_patho %>%
  as.data.frame() %>%
  mutate(pval = case_when(
    p.value < 0.001 ~ "p < 0.001",
    p.value < 0.01 ~ "p < 0.01",
    p.value < 0.05 ~ "p < 0.05",
    TRUE ~ "n.s."
  ))

cor_anno <- paste0("r = ", round(cor_res[1, 1], 2), ", ", cor_res[1, 3])
cor_anno <- sub("-", "\U2212", cor_anno)

as_filterd_patho$Landuse <- factor(as_filterd_patho$Landuse, levels = c("Natural", "Farm"))

p_patho_as <- ggplot(as_filterd_patho,
                     aes(x = fungitaxa_across, y = Pathotroph, color = Landuse, shape = Site)) +
  geom_point(size=3.5) + 
  scale_color_manual(values=c("Farm"="tan1", "Natural"="darkgreen"))+
  labs(x = "Distance to centroids of fungal communities",
       y = "Relative abundance of pathotroph (%)",
       color = "Land use",
       shape = "Site",
       title = "Fungal heterogeneity × Pathogenic fungi\n(across-site; w/o annotation rate < 41.7%)"
  ) + 
  theme_linedraw()+
  geom_smooth(method = "lm", aes(group = 1), color = "black", se = FALSE) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 13, face = "bold", colour = "black"),
        legend.text = element_text(size = 13, face = "bold", colour = "black"),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold", colour = "black"),
        legend.position = "none", # to do marginal plot, do not set legend here
        # legend.direction = "horizontal",
        # legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) +
  annotate("text",
           x = max(as_filterd_patho$fungitaxa_across),
           y = max(as_filterd_patho$Pathotroph),
           label = cor_anno,
           hjust = "right",
           size = 4
  )

p_patho_as_fin <- ggMarginal(p_patho_as, groupColour = TRUE, groupFill = TRUE)


### 2. PCoA (Ordination plot)
## filter samples by the annotation rate threshold
remain_samples <- anno_rates %>%
  filter(assigned_percent >= annorate_thres) %>%
  pull(Sample)

fungilife.t_filt <- fungilife.t[rownames(fungilife.t) %in% remain_samples, ]
landuse_filt <- env_category %>% filter(Sample %in% remain_samples) %>% pull(Landuse)
site_filt <- env_category %>% filter(Sample %in% remain_samples) %>% pull(Site)

## run PERMANOVA
# set the method to generate random values
set.seed(123)

# fungi
fungi_func_perm <- adonis2(fungilife.t_filt~site_filt*landuse_filt, method = "bray", permutations = 100000) # permanova
fungi_func_permsummary <- as.matrix(fungi_func_perm[1:3,5]) # extract p-values
rownames(fungi_func_permsummary) <- c("Site", "Land use", "Site × Land use")
colnames(fungi_func_permsummary) <- "p-value"
fungi_func_perm_result <- changep(fungi_func_permsummary) # convert the p-values to asterisks

## PCoA ordination
# fungal functions
f_func_dist <- vegdist(fungilife.t_filt, method = "bray")

f_func_pcoa <- cmdscale(f_func_dist, k = 2, eig = TRUE)
f_pcoa_scores <- scores(f_func_pcoa)

# calculate eigen values
f_eig_vals <- f_func_pcoa$eig
f_prop_explained <- f_eig_vals / sum(f_eig_vals[f_eig_vals > 0]) # positive values
f_percent_explained <- round(f_prop_explained[1:2] * 100, 1)

# format data frame
f_pcoa_scores <- f_pcoa_scores %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%# make sample column
  mutate(Site = site_filt, Landuse = landuse_filt)

# set the levels for land use factor
f_pcoa_scores$Landuse <- factor(levels = c("Natural", "Farm"), f_pcoa_scores$Landuse)

# set maximum values to annotate text (significance table)
f_pc_yRoof <- (max(f_pcoa_scores$Dim2))

# set axis titles
f_xtitle <- paste0("Axis 1 (", f_percent_explained[1], "%)")
f_ytitle <- paste0("Axis 2 (", f_percent_explained[2], "%)")

## plot
gf_pcoa <- ggplot() +
  geom_point(data=f_pcoa_scores, aes(x=Dim1, y=Dim2, fill=Landuse, shape=Site), size=4) +
  scale_fill_manual(values=c("Farm"="tan1", "Natural"="darkgreen"), guide = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("D"=21, "E"=22, "F"=23, "G"=24, "H"=25)) +
  theme_bw() +
  theme(panel.border = element_rect(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold", colour = "black"),
        legend.text = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size=15, face = "bold", colour = "black"),
        axis.text=element_text(size=10, face = "bold", colour = "black"),
        axis.title=element_text(size=10,face="bold", colour = "black"))+
  labs(title = "Fungal lifestyles (w/o annotation rate < 41.7%)", fill = "Land use", shape = "Site",
       x = f_xtitle, y = f_ytitle) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN)+
  annotate("text", x=min(f_pcoa_scores$Dim1), y=f_pc_yRoof*0.9,
           label=fungi_func_perm_result,
           hjust="left",
           size = 4)
plot(gf_pcoa)


### 3. Mantel correlogram
# => other file



### save results
# correlation plot (between fungal taxonomic heterogeneity x pathogen abundance)
saveRDS(p_patho_ws, "FigCode/FigS_FungalTraits_summary_out/cor_taxahetero_pathogen_ws.rds")
saveRDS(p_patho_as_fin, "FigCode/FigS_FungalTraits_summary_out/cor_taxahetero_pathogen_as.rds")
# PCoA plot of fungal lifestyles
saveRDS(gf_pcoa, "FigCode/FigS_FungalTraits_summary_out/pcoa_filtered.rds")


