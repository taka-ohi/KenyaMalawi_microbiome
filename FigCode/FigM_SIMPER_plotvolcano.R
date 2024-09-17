####
#### R script for Ohigashi et al (2024)
#### volcano plot for the valued analyzed in SIMPER
#### 2024.09.13 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(stringr); packageVersion("stringr")
library(ggplot2); packageVersion("ggplot2")
library(ggrepel); packageVersion("ggrepel")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data
simper_res <- read.table("11_SIMPER_out/simper_aggregated_signif_lifestyle.txt", header = T, sep = "\t")


### format data
# capitalize primary lifestyle and remove "_"
simper_res <- simper_res |>
  mutate(primary_lifestyle = str_replace_all(primary_lifestyle, "_", " ") |>
           str_to_title())

# remove "unassigned"
simper_res <- simper_res |> filter(primary_lifestyle != "Unassigned")


### create plot
## considering only land use
simper_res$significant <- factor(simper_res$significant, levels = c("Natural", "Farm", "Not significant"))

p_volcano <- ggplot(simper_res, aes(x = log2FC, y = total_ave_contrib_percent)) +
  geom_point(aes(color = significant), size = 3) + 
  scale_color_manual(values = c("Farm" = "tan1", "Natural" = "lightgreen", "Nat significant" = "grey80")) +
  theme_classic() +
  labs(x = "log2 Fold Change", y = "Contribution to dissimilarity between land uses (%)", color = "> 2-fold larger in:"#,
       # title = "Volcano Plot"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey") + 
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey") +
  scale_x_continuous(label = scaleFUN) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, face = "bold", colour = "black")
  ) +
  geom_text_repel(data = subset(simper_res, total_ave_contrib_percent > 1), 
                  aes(label = primary_lifestyle), 
                  size = 4, 
                  box.padding = 0.35, 
                  point.padding = 0.3, 
                  segment.color = 'black')
p_volcano


## considering land use in each site (added on 2024.09.15)
site <- c("D", "E", "F", "G", "H")
p_vol_sites <- list()

for (i in site) {
  # load data for each site
  simres_site <- read.table(sprintf("11_SIMPER_out/simper_aggregated_signif_lifestyle_site%s.txt", i), header = T, sep = "\t")
  
  # capitalize primary lifestyle and remove "_"
  simres_site <- simres_site |>
    mutate(primary_lifestyle = str_replace_all(primary_lifestyle, "_", " ") |>
             str_to_title())
  
  # remove "unassigned"
  simres_site <- simres_site |> filter(primary_lifestyle != "Unassigned")
  
  
  ### create plot
  ## considering only land use
  simres_site$significant <- factor(simres_site$significant, levels = c("Natural", "Farm", "Not significant"))
  
  p_vol_site <- ggplot(simres_site, aes(x = log2FC, y = total_ave_contrib_percent)) +
    geom_point(aes(color = significant), size = 3) + 
    scale_color_manual(values = c("Farm" = "tan1", "Natural" = "lightgreen", "Nat significant" = "grey80")) +
    theme_classic() +
    labs(x = "log2 Fold Change", y = "Contribution to dissimilarity between land uses (%)", color = "> 2-fold larger in:"#,
         # title = "Volcano Plot"
    ) +
    labs(title = sprintf("Site %s", i)) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey") + 
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey") +
    scale_x_continuous(label = scaleFUN) +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 12, face = "bold", colour = "black"),
      legend.text = element_text(size = 12, face = "bold", colour = "black"),
      title = element_text(size = 12, face = "bold", colour = "black")
    ) +
    geom_text_repel(data = subset(simres_site, total_ave_contrib_percent > 1), 
                    aes(label = primary_lifestyle), 
                    size = 4, 
                    box.padding = 0.35, 
                    point.padding = 0.3, 
                    segment.color = 'black')
  p_vol_sites[[i]] <- p_vol_site
}



### save data
# create a directory
dir.create("FigCode/FigM_SIMPER_out")

# save png and object
# all sites
saveRDS(p_volcano, file = "FigCode/FigM_SIMPER_out/SIMPER_allsites.obj")
ggsave(p_volcano, file = "FigCode/FigM_SIMPER_out/SIMPER_allsites.png", width = 8, height = 7, bg = "white")

# each site
for (i in site) {
  saveRDS(p_vol_sites[[i]], file = sprintf("FigCode/FigM_SIMPER_out/SIMPER_site%s.obj", i))
  ggsave(p_vol_sites[[i]], file = sprintf("FigCode/FigM_SIMPER_out/SIMPER_site%s.png", i),
         bg = "white", width = 8, height = 7)
}


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigM_SIMPER_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))




