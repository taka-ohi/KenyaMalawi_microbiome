####
#### R script for Ohigashi et al (2024)
#### PERMANOVA and PCoA analysis for the microbial community
#### 2025.06.19 written by Ohigashi
#### R 4.3.3
####


### load packages
source("Function/F2_HelperFunctions_for_Visualization.R")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(ggrepel); packageVersion("ggrepel")
library(tibble); packageVersion("tibble")


### load data
# prokaryotes
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
b_ASV <- b_ASV.table[,1:(ncol(b_ASV.table)-7)]
b_ASV.t <- t(b_ASV)

# fungi
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)
f_ASV <- f_ASV.table[,1:(ncol(f_ASV.table)-7)]
f_ASV.t <- t(f_ASV)

# environment data
env_data <- read.table("Data/soil_metadata.txt", header = T)
landuse <- env_data$Landuse
site <- env_data$Site

# prokaryotic function data
b_func <- read.table("02_Function_analysis_out/PICRUSt2_aggregated_CN_function.txt", header = T, sep = "\t")

# fungal function data
f_func <- read.table("02_Function_analysis_out/FungalTraits_specific_functions.txt", header = T, sep = "\t")
f_func <- f_func |> column_to_rownames(var = "Sample") # make the sample names rownames
f_func <- f_func |> dplyr::select(ECM, Pathotroph, Saprotroph) # use only data from "primary_lifestyles"


### PERMANOVA
# set the method to generate random values
set.seed(123)

# prokaryotes
b_perm <- adonis2(b_ASV.t~site*landuse, method = "bray", permutations = 100000) # permanova
b_permsummary <- as.matrix(b_perm[1:3,5]) # extract p-values
rownames(b_permsummary) <- c("Site", "Land use", "Site × Land use")
colnames(b_permsummary) <- "p-value"
b_perm_result <- changep(b_permsummary) # convert the p-values to asterisks

# fungi
f_perm <- adonis2(f_ASV.t~site*landuse, method = "bray", permutations = 100000) # permanova
f_permsummary <- as.matrix(f_perm[1:3,5]) # extract p-values
rownames(f_permsummary) <- c("Site", "Land use", "Site × Land use")
colnames(f_permsummary) <- "p-value"
f_perm_result <- changep(f_permsummary) # convert the p-values to asterisks


### PCoA
# calculate bray curtis dissimilarity
# prokaryotes
b_taxa_dist <- vegdist(b_ASV.t, method = "bray")
# fungai
f_taxa_dist <- vegdist(f_ASV.t, method = "bray")

# PCoA ordination
## prokaryotes
b_taxa_pcoa <- cmdscale(b_taxa_dist, k = 2, eig = TRUE)
b_pcoa_scores <- scores(b_taxa_pcoa)

# calculate eigen values
b_eig_vals <- b_taxa_pcoa$eig
b_prop_explained <- b_eig_vals / sum(b_eig_vals[b_eig_vals > 0]) # positive values
b_percent_explained <- round(b_prop_explained[1:2] * 100, 1)

## fungi
f_taxa_pcoa <- cmdscale(f_taxa_dist, k = 2, eig = TRUE)
f_pcoa_scores <- scores(f_taxa_pcoa)

# calculate eigen values
f_eig_vals <- f_taxa_pcoa$eig
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


### check correlations with microbial functions
## prokaryotes
# fitting functions and testing the correlation with permutation
b_ef_func <- envfit(b_pcoa_scores[,2:3], b_func, permu = 100000) # envfit

# formatting
b_ef_func_pvals <- as.data.frame(b_ef_func$vectors$pvals) # extract p-values
b_ef_func_r <- as.data.frame(b_ef_func$vectors$r) # extract r2
b_ef_func_arrows <- as.data.frame(b_ef_func$vectors$arrows*sqrt(b_ef_func$vectors$r)) # add weight by r
b_ef_func_table <- cbind(b_ef_func_arrows, b_ef_func_pvals); colnames(b_ef_func_table)[3] <- "pvals"
b_ef_func_table_sig <- subset(b_ef_func_table, pvals < 0.05)
b_ef_func_table_sig <- b_ef_func_table_sig |> rownames_to_column(var = "factor")
b_ef_func_table_sig$Type <- c(rep("function", nrow(b_ef_func_table_sig))) # for later use

## fungi
# fitting functions and testing the correlation with permutation
f_ef_func <- envfit(f_pcoa_scores[,2:3], f_func, permu = 100000) # envfit

# formatting
f_ef_func_pvals <- as.data.frame(f_ef_func$vectors$pvals) # extract p-values
f_ef_func_r <- as.data.frame(f_ef_func$vectors$r) # extract r2
f_ef_func_arrows <- as.data.frame(f_ef_func$vectors$arrows*sqrt(f_ef_func$vectors$r)) # add weight by r
f_ef_func_table <- cbind(f_ef_func_arrows, f_ef_func_pvals); colnames(f_ef_func_table)[3] <- "pvals"
f_ef_func_table_sig <- subset(f_ef_func_table, pvals < 0.05)
f_ef_func_table_sig <- f_ef_func_table_sig |> rownames_to_column(var = "factor")
f_ef_func_table_sig$Type <- c(rep("function", nrow(f_ef_func_table_sig))) # for later use

### check correlations with environmental variables
# select environment variables
env_fct <- env |> dplyr::select(Gravimetric.water.content,
                                Carbon,
                                Nitrogen,
                                CN_ratio,
                                pH
)

## prokaryotes
# fitting functions and testing the correlation with permutation
b_ef_env <- envfit(b_pcoa_scores[,2:3], env_fct, permu = 100000) # envfit

# formatting
b_ef_env_pvals <- as.data.frame(b_ef_env$vectors$pvals) # extract p-values
b_ef_env_r <- as.data.frame(b_ef_env$vectors$r) # extract r2
b_ef_env_arrows <- as.data.frame(b_ef_env$vectors$arrows*sqrt(b_ef_env$vectors$r)) # add weight by r
b_ef_env_table <- cbind(b_ef_env_arrows, b_ef_env_pvals); colnames(b_ef_env_table)[3] <- "pvals"
b_ef_env_table_sig <- subset(b_ef_env_table, pvals < 0.05)
b_ef_env_table_sig <- b_ef_env_table_sig |> rownames_to_column(var = "factor")
b_ef_env_table_sig$Type <- c(rep("environment", nrow(b_ef_env_table_sig))) # for later use

## fungi
# fitting functions and testing the correlation with permutation
f_ef_env <- envfit(f_pcoa_scores[,2:3], env_fct, permu = 100000) # envfit

# formatting
f_ef_env_pvals <- as.data.frame(f_ef_env$vectors$pvals) # extract p-values
f_ef_env_r <- as.data.frame(f_ef_env$vectors$r) # extract r2
f_ef_env_arrows <- as.data.frame(f_ef_env$vectors$arrows*sqrt(f_ef_env$vectors$r)) # add weight by r
f_ef_env_table <- cbind(f_ef_env_arrows, f_ef_env_pvals); colnames(f_ef_env_table)[3] <- "pvals"
f_ef_env_table_sig <- subset(f_ef_env_table, pvals < 0.05)
f_ef_env_table_sig <- f_ef_env_table_sig |> rownames_to_column(var = "factor")
f_ef_env_table_sig$Type <- c(rep("environment", nrow(f_ef_env_table_sig))) # for later use

## combine tables for environmental factors and functions
b_ef_table_sig <- rbind(b_ef_func_table_sig, b_ef_env_table_sig)
f_ef_table_sig <- rbind(f_ef_func_table_sig, f_ef_env_table_sig)

# change factor names for plotting
b_ef_table_sig <- b_ef_table_sig |>
  mutate(Guild =
           case_when(
             factor == "Nitrification" ~ "Nitrifying",
             factor == "Cellulose_breakdown" ~ "Cellulolytic",
             factor == "Chitin_breakdown" ~ "Chitinolytic",
             factor == "Autotrophic_C_fixation" ~ "Autotroph",
             factor == "Denitrification" ~ "Denitrifying",
             factor == "Lignin_breakdown" ~ "Lignolytic",
             factor == "Methane_oxidation" ~ "Methanotroph",
             factor == "Methanogenesis" ~ "Methanogen",
             factor == "N_fixation" ~ "N-fixing",
             factor == "Xylan_breakdown" ~ "Xylanolytic",
             factor == "Gravimetric.water.content" ~ "Moisture",
             TRUE ~ factor
           ))
f_ef_table_sig <- f_ef_table_sig |>
  mutate(Guild =
           case_when(
             factor == "Gravimetric.water.content" ~ "Moisture",
             TRUE ~ factor
           ))


# set the levels for land use factor
b_pcoa_scores$Landuse <- factor(levels = c("Natural", "Farm"), b_pcoa_scores$Landuse)
f_pcoa_scores$Landuse <- factor(levels = c("Natural", "Farm"), f_pcoa_scores$Landuse)

# set maximum values to annotate text (significance table)
b_pc_yRoof <- (min(b_pcoa_scores$Dim2))
f_pc_yRoof <- (min(f_pcoa_scores$Dim2))

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
  labs(title = "Prokaryotes", fill = "Land use", shape = "Site",
       x = b_xtitle, y = b_ytitle
  ) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN)+
  annotate("text", x=min(b_pcoa_scores$Dim1)*1.1, y=b_pc_yRoof*0.9,
           label=b_perm_result,
           hjust="left",
           size = 4) +
  geom_segment(data = b_ef_table_sig, aes(x=0, y=0, xend=Dim1*0.3, yend=Dim2*0.3, color= Type),
               size=0.5,arrow=arrow(type = "open", length = unit(0.08, "inches")), show.legend = F) +
  geom_text_repel(data = b_ef_table_sig, aes(x = Dim1*0.4, y = Dim2*0.4, label = Guild, color=Type),
                  size = 4, show.legend = F)+
  scale_color_manual(name = "Type", values = c("function"="skyblue", "environment"="red"))
plot(gb_pcoa)


# fungi
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
  labs(title = "Fungi", fill = "Land use", shape = "Site",
       x = f_xtitle, y = f_ytitle) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN)+
  annotate("text", x=max(f_pcoa_scores$Dim1), y=f_pc_yRoof*0.9,
           label=f_perm_result,
           hjust="right",
           size = 4)+
  geom_segment(data = f_ef_table_sig, aes(x=0, y=0, xend=Dim1*0.3, yend=Dim2*0.3, color= Type),
               size=0.5,arrow=arrow(type = "open", length = unit(0.08, "inches")), show.legend = F) +
  geom_text_repel(data = f_ef_table_sig, aes(x = Dim1*0.4, y = Dim2*0.4, label = Guild, color=Type),
                  size = 4, show.legend = F)+
  scale_color_manual(name = "Type", values = c("function"="skyblue", "environment"="red"))
plot(gf_pcoa)



### save data
dir.create("04_PCoA_PERMANOVA_out")
# permanova result
write.csv(b_perm, file = "04_PCoA_PERMANOVA_out/permanova_result_prok.csv", quote = F, row.names = T)
write.csv(f_perm, file = "04_PCoA_PERMANOVA_out/permanova_result_fungi.csv", quote = F, row.names = T)

# save RDS
saveRDS(gb_pcoa, file = "04_PCoA_PERMANOVA_out/PCoA_prok_taxa.rds")
ggsave("04_PCoA_PERMANOVA_out/PCoA_prok_taxa.png", plot = gb_pcoa,
       width = 8, height = 7, bg = "white")
saveRDS(gf_pcoa, file = "04_PCoA_PERMANOVA_out/PCoA_fungi_taxa.rds")
ggsave("04_PCoA_PERMANOVA_out/PCoA_fungi_taxa.png", plot = gf_pcoa,
       width = 8, height = 7, bg = "white")

# save PCoA data frame
saveRDS(b_pcoa_scores, "04_PCoA_PERMANOVA_out/pcoa_scores_prok_df.rds")
saveRDS(f_pcoa_scores, "04_PCoA_PERMANOVA_out/pcoa_scores_fungi_df.rds")

# save envfit results
# with function
saveRDS(b_ef_func, file = "04_PCoA_PERMANOVA_out/envfit_function_prok.rds")
saveRDS(f_ef_func, file = "04_PCoA_PERMANOVA_out/envfit_function_fungi.rds")
# with environmental variables
saveRDS(b_ef_env, file = "04_PCoA_PERMANOVA_out/envfit_env_prok.rds")
saveRDS(f_ef_env, file = "04_PCoA_PERMANOVA_out/envfit_env_fungi.rds")

