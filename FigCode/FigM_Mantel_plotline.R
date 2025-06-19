####
#### R script for Ohigashi et al (2024)
#### visualization of Mantel correlogram
#### 2024.09.13 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(stringr); packageVersion("stringr")
library(ggplot2); packageVersion("ggplot2")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data
# Manterl correlograms for microbial communities
prok_man <- read.table("10_Mantel_correlogram_out/mantel_result_prok.txt", header = T)
fungi_man <- read.table("10_Mantel_correlogram_out/mantel_result_fungi.txt", header = T)
# Manterl correlograms for microbial functions
prokfunc_man <- read.table("10_Mantel_correlogram_out/mantel_result_prokallfunc.txt", header = T)
fungilife_man <- read.table("10_Mantel_correlogram_out/mantel_result_fungilife.txt", header = T)
# distance category
dist_cat <- read.table("10_Mantel_correlogram_out/distance_category.txt", header = T)


### format data
# add "dist.class" & # add "type" columns for later use
prok_man <- prok_man |>
  mutate(dist.class = rep(c("1", "2", "3", "4", "5"), 2),
         type = rep("Taxa", nrow(prok_man)))
fungi_man <- fungi_man |>
  mutate(dist.class = rep(c("1", "2", "3", "4", "5"), 2),
         type = rep("Taxa", nrow(fungi_man)))
prokfunc_man <- prokfunc_man |>
  mutate(dist.class = rep(c("1", "2", "3", "4", "5"), 2),
         type = rep("Functions", nrow(prokfunc_man)))
fungilife_man <- fungilife_man |>
  mutate(dist.class = rep(c("1", "2", "3", "4", "5"), 2),
         type = rep("Functions", nrow(fungilife_man)))

# combine taxa and function data
prokdf <- rbind(prok_man, prokfunc_man)
fungidf <- rbind(fungi_man, fungilife_man)


### create plots
## prokaryotes
# set levels for factors
prokdf$landuse <- factor(prokdf$landuse, levels = c("Natural", "Farm"))
prokdf$type <- factor(prokdf$type, levels = c("Taxa", "Functions"))
prok_p <- ggplot(prokdf, aes(x = dist.class, y = Mantel.cor,
                              color = landuse, shape = type, group = interaction(landuse, type), linetype = type)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c("darkgreen", "tan1")) +
  scale_shape_manual(values = c(16, 2)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(
    x = "Distance Class",
    y = "Mantel's r",
    title = "Correlation by distance class\n(Prokaryotes and Prokaryotic functions)",
    color = "Land use",
    linetype = "Correlation in:",
    shape = "Correlation in:"
  ) +
  scale_y_continuous(labels = scaleFUN) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", colour = "black", hjust = 0.5),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, face = "bold", colour = "black")
  ) +
  geom_text(
    data = prokdf |> filter(Pr.corrected. < 0.05 ),
    aes(label = "*", y = Mantel.cor + 0.02),
    color = "black",
    size = 4
  ) +
  annotate("text", x = 1, y = min(prokdf$Mantel.cor)*1.1, label = "Within-site", size = 4) +  
  annotate("text", x = 2.5, y = min(prokdf$Mantel.cor)*1.1, label = "Across-site\n(within country)", size = 4) +  
  annotate("text", x = 4.5, y = min(prokdf$Mantel.cor)*1.1, label = "Across-site\n(between countries)", size = 4)

prok_p


## Fungi
# set levels for factors
fungidf$landuse <- factor(fungidf$landuse, levels = c("Natural", "Farm"))
fungidf$type <- factor(fungidf$type, levels = c("Taxa", "Functions"))
fungi_p <- ggplot(fungidf, aes(x = dist.class, y = Mantel.cor,
                             color = landuse, shape = type, group = interaction(landuse, type), linetype = type)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c("darkgreen", "tan1")) +
  scale_shape_manual(values = c(16, 2)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(
    x = "Distance Class",
    y = "Mantel's r",
    title = "Correlation by distance class\n(Fungi and Fungal lifestyles)",
    color = "Land use",
    linetype = "Correlation in:",
    shape = "Correlation in:"
  ) +
  scale_y_continuous(labels = scaleFUN) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", colour = "black", hjust = 0.5),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, face = "bold", colour = "black")
  ) +
  geom_text(
    data = fungidf |> filter(Pr.corrected. < 0.05 ),
    aes(label = "*", y = Mantel.cor + 0.02),
    color = "black",
    size = 4
  ) +
  annotate("text", x = 1, y = min(fungidf$Mantel.cor)*1.1, label = "Within-site", size = 4) +  
  annotate("text", x = 2.5, y = min(fungidf$Mantel.cor)*1.1, label = "Across-site\n(within country)", size = 4) +  
  annotate("text", x = 4.5, y = min(fungidf$Mantel.cor)*1.1, label = "Across-site\n(between countries)", size = 4)

fungi_p


### save data
# create a directory
dir.create("FigCode/FigM_Mantel_out")

# save object and image files
# prokaryotes
saveRDS(prok_p, file = "FigCode/FigM_Mantel_out/Mantel_prok_taxafunc.obj")
ggsave(prok_p, file = "FigCode/FigM_Mantel_out/Mantel_prok_taxafunc.png",
       bg = "white", width = 8, height = 7)
# fungi
saveRDS(fungi_p, file = "FigCode/FigM_Mantel_out/Mantel_fungi_taxafunc.obj")
ggsave(fungi_p, file = "FigCode/FigM_Mantel_out/Mantel_fungi_taxafunc.png",
       bg = "white", width = 8, height = 7)


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigM_Mantel_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))



