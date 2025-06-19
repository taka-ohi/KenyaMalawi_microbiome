####
#### R script for Ohigashi et al (2024)
#### PERMANOVA and PCoA analysis for the microbial functions
#### 2025.06.10 written by Ohigashi
#### R 4.3.3
#### 


### load packages and functions
source("Function/F2_HelperFunctions_for_Visualization.R")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(ggrepel); packageVersion("ggrepel")
library(tibble); packageVersion("tibble")


### load data
# prokaryotic function data
prok_func_all <- read.table("02_Function_analysis_out/PICRUSt2_full_function.txt", header = T)
# fungal function data
fungilife <- read.table("02_Function_analysis_out/FungalTraits_primarilifestyle_percent.txt", header = T)
fungilife <- fungilife |> filter(primary_lifestyle != "unassigned") # 2024.09.11
fungilife <- fungilife |> column_to_rownames(var = "primary_lifestyle")
fungilife.t <- t(fungilife)

# environment data
env_data <- read.table("Data/soil_metadata.txt", header = T)
landuse <- env_data$Landuse
site <- env_data$Site


### run PERMANOVA
# set the method to generate random values
set.seed(123)

# prokaryotic functions
prok_func_perm <- adonis2(prok_func_all~site*landuse, method = "bray", permutations = 100000) # permanova
prok_func_permsummary <- as.matrix(prok_func_perm[1:3,5]) # extract p-values
rownames(prok_func_permsummary) <- c("Site", "Land use", "Site × Land use")
colnames(prok_func_permsummary) <- "p-value"
prok_func_perm_result <- changep(prok_func_permsummary) # convert the p-values to asterisks

# fungi
fungi_func_perm <- adonis2(fungilife.t~site*landuse, method = "bray", permutations = 100000) # permanova
fungi_func_permsummary <- as.matrix(fungi_func_perm[1:3,5]) # extract p-values
rownames(fungi_func_permsummary) <- c("Site", "Land use", "Site × Land use")
colnames(fungi_func_permsummary) <- "p-value"
fungi_func_perm_result <- changep(fungi_func_permsummary) # convert the p-values to asterisks


### PCoA
# prokaryotic functions
b_func_dist <- vegdist(prok_func_all, method = "bray")
# fungal functions
f_func_dist <- vegdist(fungilife.t, method = "bray")

## PCoA ordination
# prokaryotic functions
b_func_pcoa <- cmdscale(b_func_dist, k = 2, eig = TRUE)
b_pcoa_scores <- scores(b_func_pcoa)

# calculate eigen values
b_eig_vals <- b_func_pcoa$eig
b_prop_explained <- b_eig_vals / sum(b_eig_vals[b_eig_vals > 0]) # positive values
b_percent_explained <- round(b_prop_explained[1:2] * 100, 1)


# fungal lifestyles
f_func_pcoa <- cmdscale(f_func_dist, k = 2, eig = TRUE)
f_pcoa_scores <- scores(f_func_pcoa)

# calculate eigen values
f_eig_vals <- f_func_pcoa$eig
f_prop_explained <- f_eig_vals / sum(f_eig_vals[f_eig_vals > 0]) # positive values
f_percent_explained <- round(f_prop_explained[1:2] * 100, 1)


# format data frame
b_pcoa_scores <- b_pcoa_scores %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%# make sample column
  mutate(Site = site, Landuse = landuse)
f_pcoa_scores <- f_pcoa_scores %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%# make sample column
  mutate(Site = site, Landuse = landuse)

# set the levels for land use factor
b_pcoa_scores$Landuse <- factor(levels = c("Natural", "Farm"), b_pcoa_scores$Landuse)
f_pcoa_scores$Landuse <- factor(levels = c("Natural", "Farm"), f_pcoa_scores$Landuse)

# set maximum values to annotate text (significance table)
b_pc_yRoof <- (max(b_pcoa_scores$Dim2))
f_pc_yRoof <- (max(f_pcoa_scores$Dim2))

# set axis titles
b_xtitle <- paste0("Axis 1 (", b_percent_explained[1], "%)")
b_ytitle <- paste0("Axis 2 (", b_percent_explained[2], "%)")

f_xtitle <- paste0("Axis 1 (", f_percent_explained[1], "%)")
f_ytitle <- paste0("Axis 2 (", f_percent_explained[2], "%)")

## plot
# prokaryotes
gb_pcoa <- ggplot() +
  geom_point(data=b_pcoa_scores, aes(x=Dim1, y=Dim2, fill=Landuse, shape=Site), size=4) +
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
        axis.text=element_text(size=12, colour = "black"),
        axis.title=element_text(size=12, colour = "black"))+
  labs(title = "Prokaryotic functions", fill = "Land use", shape = "Site",
       x = b_xtitle, y = b_ytitle
       ) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN)+
  annotate("text", x=max(b_pcoa_scores$Dim1), y=b_pc_yRoof*0.9,
           label=prok_func_perm_result,
           hjust="right",
           size = 4)
plot(gb_pcoa)


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
        axis.text=element_text(size=12, colour = "black"),
        axis.title=element_text(size=12, colour = "black"))+
  labs(title = "Fungal lifestyles", fill = "Land use", shape = "Site",
       x = f_xtitle, y = f_ytitle) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN)+
  annotate("text", x=min(f_pcoa_scores$Dim1), y=f_pc_yRoof*0.9,
           label=fungi_func_perm_result,
           hjust="left",
           size = 4)
plot(gf_pcoa)


### save results
# permutation results
# permanova result
write.csv(prok_func_perm, file = "04_PCoA_PERMANOVA_out/permanova_result_prokfunc.csv", quote = F, row.names = T)
write.csv(fungi_func_perm, file = "04_PCoA_PERMANOVA_out/permanova_result_fungilife.csv", quote = F, row.names = T)

# PCoA plots
saveRDS(gb_pcoa, "04_PCoA_PERMANOVA_out/PCoA_prok_functions.rds")
ggsave("04_PCoA_PERMANOVA_out/PCoA_prok_functions.png", plot = gb_pcoa,
       width = 8, height = 7, bg = "white")
saveRDS(gf_pcoa, "04_PCoA_PERMANOVA_out/PCoA_fungi_functions.rds")
ggsave("04_PCoA_PERMANOVA_out/PCoA_fungi_functions.png", plot = gf_pcoa,
       width = 8, height = 7, bg = "white")
