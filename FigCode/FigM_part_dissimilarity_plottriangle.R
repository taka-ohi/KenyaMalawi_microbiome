####
#### R script for Ohigashi et al (2024)
#### triangle plots for total dissimilarity, replacement, and richness difference
#### 2024.09.12 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(tidyr); packageVersion("tidyr")
library(reshape2); packageVersion("reshape2")
library(ggplot2); packageVersion("ggplot2")
library(ggtern); packageVersion("ggtern")
library(cowplot); packageVersion("cowplot")
library(ggpubr); packageVersion("ggpubr")
source("Function/F1_HelperFunctions.R")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data (created in 09_partitioning_dissimilarity.R)
# prokaryotes dissimilarity
prok_df <- read.table("08_partitioning_dissimilarity_out/paritionied_dissimilarity_prok.txt", header = T, sep = "\t")
prok_df <- prok_df |>
  filter(!landuse %in% c("Natural_Farm", "Farm_Natural")) |> # remove pair containing different land uses.
  select(Similarity, Repl, RichDiff, landuse, relation)

# fungi dissimilarity
fungi_df <- read.table("08_partitioning_dissimilarity_out/paritionied_dissimilarity_fungi.txt", header = T, sep = "\t")
fungi_df <- fungi_df |>
  filter(!landuse %in% c("Natural_Farm", "Farm_Natural")) |> # remove pair containing different land uses.
  select(Similarity, Repl, RichDiff, landuse, relation)

# format pair type
prok_df <- prok_df |>
  mutate(landuse2 = case_when(
    landuse == "Natural_Natural" ~ "Natural vs Natural",
    landuse == "Farm_Farm" ~ "Farm vs Farm",
    TRUE ~ NA
  ))

fungi_df <- fungi_df |>
  mutate(landuse2 = case_when(
    landuse == "Natural_Natural" ~ "Natural vs Natural",
    landuse == "Farm_Farm" ~ "Farm vs Farm",
    TRUE ~ NA
  ))


### creating plots
scale_cat <- c("within_site", "diff_site_within_country", "diff_site_diff_country")

## looping
# prokaryotes
prok_plots <- list()
for (i in scale_cat) {
  # subset
  subdf <- prok_df |> filter(relation == i)
  subdf$landuse2 <- factor(subdf$landuse2, levels = c("Natural vs Natural", "Farm vs Farm"))
  
  mean_vals <- subdf |>
    group_by(landuse2) |>
    summarise(across(1:3, mean))
  
  # set plot titles
  if (i == "within_site") {
    plottitle <- "Prokaryotes\n(within-site)"
  } else if (i == "diff_site_within_country") {
    plottitle <- "Prokaryotes\n(across-site, within the same country)"
  } else {
    plottitle <- "Prokaryotes\n(across-site, between the countries)"
  }
  
  # plot
  p <- ggtern(data=subdf, aes(x=Similarity, y=Repl, z=RichDiff)) +
    geom_point(aes(color=landuse2), shape = 4, alpha = 0.3)+
    geom_crosshair_tern(data=mean_vals, aes(color = landuse2)) +
    geom_point(data = mean_vals, aes(color = landuse2), size = 2)+
    scale_color_manual(values = c("darkgreen", "tan1"))+
    guides(fill = "none", alpha = "none")+
    theme_bw(base_size = 13) +
    theme_showarrows() +
    labs(x="Similarity",
         y="Replacement",
         z="Richness difference",
         title=plottitle,
         color="Pair type"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
          )+
    theme_hidetitles()
  
  prok_plots[[i]] <- p
}


# fungi
fungi_plots <- list()
for (i in scale_cat) {
  # subset
  subdf <- fungi_df |> filter(relation == i)
  subdf$landuse2 <- factor(subdf$landuse2, levels = c("Natural vs Natural", "Farm vs Farm"))
  
  mean_vals <- subdf |>
    group_by(landuse2) |>
    summarise(across(1:3, mean))
  
  # set plot titles
  if (i == "within_site") {
    plottitle <- "Fungi\n(within-site)"
  } else if (i == "diff_site_within_country") {
    plottitle <- "Fungi\n(across-site, within the same country)"
  } else {
    plottitle <- "Fungi\n(across-site, between the countries)"
  }
  
  # plot
  p <- ggtern(data=subdf, aes(x=Similarity, y=Repl, z=RichDiff)) +
    geom_point(aes(color=landuse2), shape = 4, alpha = 0.3)+
    geom_crosshair_tern(data=mean_vals, aes(color = landuse2)) +
    geom_point(data = mean_vals, aes(color = landuse2), size = 2)+
    scale_color_manual(values = c("darkgreen", "tan1"))+
    guides(fill = "none", alpha = "none")+
    theme_bw(base_size = 13) +
    theme_showarrows() +
    labs(x="Similarity",
         y="Replacement",
         z="Richness difference",
         title=plottitle,
         color="Pair type"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )+
    theme_hidetitles()
  
  fungi_plots[[i]] <- p
}


## combine plots for the main article
# get legend
p_legend <- get_legend(prok_plots[[1]])
plot(p_legend)

# make w/o-legend version of figures
prok_plots2 <- list()
fungi_plots2 <- list()
for (i in scale_cat) {
  prok_plots2[[i]] <- prok_plots[[i]] + theme(legend.position = "none")
  fungi_plots2[[i]] <- fungi_plots[[i]] + theme(legend.position = "none")
}
allplots <- c(prok_plots2, fungi_plots2)

# make a combined version (cowplot did not work for ggtern objects)
fig_all <- ggarrange(arrangeGrob(grobs = allplots[1]), arrangeGrob(grobs = allplots[2]), arrangeGrob(grobs = allplots[3]), p_legend,
                  arrangeGrob(grobs = allplots[4]), arrangeGrob(grobs = allplots[5]), arrangeGrob(grobs = allplots[6]),
                  ncol = 4, nrow = 2, labels = c("a", "b", "c", NA, "d", "e", "f"), font.label = list(size = 20)
                  )


### save data
# create a directory
dir.create("FigCode/FigM_part_dissimilarity_out")

# save object and image files
for (i in scale_cat) {
  # prokaryotes
  saveRDS(prok_plots[[i]], file = sprintf("FigCode/FigM_part_dissimilarity_out/%s_triangle_prok.obj", i))
  ggsave(prok_plots[[i]], file = sprintf("FigCode/FigM_part_dissimilarity_out/%s_triangle_prok.png", i),
         bg = "white", width = 8, height = 7)
  # fungi
  saveRDS(fungi_plots[[i]], file = sprintf("FigCode/FigM_part_dissimilarity_out/%s_triangle_fungi.obj", i))
  ggsave(fungi_plots[[i]], file = sprintf("FigCode/FigM_part_dissimilarity_out/%s_triangle_fungi.png", i),
         bg = "white", width = 8, height = 7)
}

ggsave("FigCode/FigM_part_dissimilarity_out/FigM_part_dissimilarity.png", plot = fig_all, width = 14, height = 8.5)


# save PDF for the main article
cairo_pdf("FigCode/FigM_part_dissimilarity_out/FigM_part_dissimilarity.pdf",  width = 14, height = 8.5)
print(fig_all)
dev.off()

# save png for the main article
ggsave("FigCode/FigM_part_dissimilarity_out/FigM_part_dissimilarity.png", plot = fig_all, width = 14, height = 8.5)

### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigM_part_dissimilarity_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


