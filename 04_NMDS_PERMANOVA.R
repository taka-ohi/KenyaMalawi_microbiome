####
#### R script for Ohigashi et al (2024)
#### PERMANOVA and NMDS analysis for the microbial community
#### 2024.07.02 written by Ohigashi
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


### NMDS for microbes
# ordination
b_ord <- metaMDS(b_ASV.t)
f_ord <- metaMDS(f_ASV.t)

# format result data
# prokaryotes
b_data.scores <- as.data.frame(scores(b_ord)[1]) # extract the result
b_data.scores <- b_data.scores |> rownames_to_column(var = "Sample") # make sample column
b_data.scores <- b_data.scores |>
  mutate(Site = site, Landuse = landuse) # make Site and Landuse columns
colnames(b_data.scores)[2:3] <- c("NMDS1", "NMDS2")

# fungi
f_data.scores <- as.data.frame(scores(f_ord)[1])  # extract the result
f_data.scores <- f_data.scores |> rownames_to_column(var = "Sample") # make sample column
f_data.scores <- f_data.scores |>
  mutate(Site = site, Landuse = landuse) # make Site and Landuse columns
colnames(f_data.scores)[2:3] <- c("NMDS1", "NMDS2")


### check correlations with microbial functions
## prokaryotes
# fitting functions and testing the correlation with permutation
b_ef_func <- envfit(b_ord, b_func, permu = 100000) # envfit

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
f_ef_func <- envfit(f_ord, f_func, permu = 100000) # envfit

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
b_ef_env <- envfit(b_ord, env_fct, permu = 100000) # envfit

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
f_ef_env <- envfit(f_ord, env_fct, permu = 100000) # envfit

# formatting
f_ef_env_pvals <- as.data.frame(f_ef_env$vectors$pvals) # extract p-values
f_ef_env_r <- as.data.frame(f_ef_env$vectors$r) # extract r2
f_ef_env_arrows <- as.data.frame(f_ef_env$vectors$arrows*sqrt(f_ef_env$vectors$r)) # add weight by r
f_ef_env_table <- cbind(f_ef_env_arrows, f_ef_env_pvals); colnames(f_ef_env_table)[3] <- "pvals"
f_ef_env_table_sig <- subset(f_ef_env_table, pvals < 0.05)
f_ef_env_table_sig <- f_ef_env_table_sig |> rownames_to_column(var = "factor")
f_ef_env_table_sig$Type <- c(rep("environment", nrow(f_ef_env_table_sig))) # for later use


### combine tables for environmental factors and functions
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
             factor == "N_fixation" ~ "Nitrogen-fixing",
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


### making NMDS plots
# set maximum values to annotate text (significance table)
b_yRoof <- (max(b_data.scores$NMDS2))
f_yRoof <- (max(f_data.scores$NMDS2))

# set the levels for land use factor
b_data.scores$Landuse <- factor(levels = c("Natural", "Farm"), b_data.scores$Landuse)
f_data.scores$Landuse <- factor(levels = c("Natural", "Farm"), f_data.scores$Landuse)

# define titles of the plots to include "stress values"
proktitle <- sprintf("Prokaryotes (Stress = %s)", round(b_ord$stress, 2))
fungititle <- sprintf("Fungi (Stress = %s)", round(f_ord$stress, 2))

# plot
# prokaryotes
gb <- ggplot() +
  geom_point(data=b_data.scores, aes(x=NMDS1, y=NMDS2, fill=Landuse, shape=Site), size=4) +
  scale_fill_manual(name = "Landuse", values=c("Farm"="orange", "Natural"="green"), guide = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(name = "Site", values = c("D"=21, "E"=22, "F"=23, "G"=24, "H"=25)) +
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
  labs(title = proktitle) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN)+
  annotate("text", x=min(b_data.scores$NMDS1), y=-b_yRoof*0.95,
           label=b_perm_result,
           hjust="left",
           size = 4)+
  geom_segment(data = b_ef_table_sig, aes(x=0, y=0, xend=NMDS1*1.2, yend=NMDS2*1.2, color= Type),
               size=0.5,arrow=arrow(type = "open", length = unit(0.08, "inches")), show.legend = F) +
  geom_text_repel(data = b_ef_table_sig, aes(x = NMDS1*1.5, y = NMDS2*1.5, label = Guild, color=Type),
                  size = 4, show.legend = F)+
  scale_color_manual(name = "Type", values = c("function"="blue", "environment"="red"))
plot(gb)

# fungi
gf <- ggplot() +
  geom_point(data=f_data.scores, aes(x=NMDS1, y=NMDS2, fill=Landuse, shape=Site),size=4) +
  scale_fill_manual(name = "Landuse", values=c("Farm"="orange", "Natural"="green"), guide = guide_legend(override.aes = list(shape=21))) +
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
  labs(title = fungititle)+
  scale_y_continuous(labels = scaleFUN2) +
  scale_x_continuous(labels = scaleFUN2)+
  annotate("text", x=min(f_data.scores$NMDS1), y=-f_yRoof*0.95,
           label=f_perm_result,
           hjust="left",
           size = 4)+
  geom_segment(data = f_ef_table_sig, aes(x=0, y=0, xend=NMDS1*0.12, yend=NMDS2*0.12, colour= Type),
               size=0.5,arrow=arrow(type = "open", length = unit(0.08, "inches")), show.legend = F) +
  geom_text_repel(data = f_ef_table_sig, aes(x = NMDS1*0.15, y = NMDS2*0.15, label = Guild, colour=Type),
                  size = 4, show.legend = F)+
  scale_color_manual(values = c("function"="blue", "environment"="red"))
plot(gf)



### save data
dir.create("04_NMDS_PERMANOVA_out")
# permanova result
write.csv(b_perm, file = "04_NMDS_PERMANOVA_out/permanova_result_prok.csv", quote = F, row.names = T)
write.csv(f_perm, file = "04_NMDS_PERMANOVA_out/permanova_result_fungi.csv", quote = F, row.names = T)

# ordination result
saveRDS(b_ord, file = "04_NMDS_PERMANOVA_out/ordination_prok.obj")
saveRDS(f_ord, file = "04_NMDS_PERMANOVA_out/ordination_fungi.obj")

# NMDS scores
write.csv(b_data.scores, file = "04_NMDS_PERMANOVA_out/NMDS_scores_prok.csv", quote = F, row.names = F)
write.csv(f_data.scores, file = "04_NMDS_PERMANOVA_out/NMDS_scores_fungi.csv", quote = F, row.names = F)

# envfit results
# with function
saveRDS(b_ef_func, file = "04_NMDS_PERMANOVA_out/envfit_function_prok.obj")
saveRDS(f_ef_func, file = "04_NMDS_PERMANOVA_out/envfit_function_fungi.obj")
# with environmental variables
saveRDS(b_ef_env, file = "04_NMDS_PERMANOVA_out/envfit_env_prok.obj")
saveRDS(f_ef_env, file = "04_NMDS_PERMANOVA_out/envfit_env_fungi.obj")

# overall factor correlation table
b_ef_table <- rbind(b_ef_func_table, b_ef_env_table)
write.csv(b_ef_table, file = "04_NMDS_PERMANOVA_out/envfit_result_prok.csv", quote = F, row.names = T)
f_ef_table <- rbind(f_ef_func_table, f_ef_env_table)
write.csv(f_ef_table, file = "04_NMDS_PERMANOVA_out/envfit_result_fungi.csv", quote = F, row.names = T)

# plot objects
saveRDS(gb, file = "04_NMDS_PERMANOVA_out/NMDS_plot_prok.obj")
saveRDS(gf, file = "04_NMDS_PERMANOVA_out/NMDS_plot_fungi.obj")

# plot images
ggsave(filename = "04_NMDS_PERMANOVA_out/NMDS_plot_prok.png", plot = gb)
ggsave(filename = "04_NMDS_PERMANOVA_out/NMDS_plot_fungi.png", plot = gf)

# Shepard plot (to check if the ordination is well fitted to ovserved data)
png(filename = "04_NMDS_PERMANOVA_out/shepard_prok.png", width = 500, height = 500, res = 100)
stressplot(b_ord, main = "Shepard plot (Prokaryotes)")
dev.off()

png(filename = "04_NMDS_PERMANOVA_out/shepard_fungi.png", width = 500, height = 500, res = 100)
stressplot(f_ord, main = "Shepard plot (Fungi)")
dev.off()


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/04_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
