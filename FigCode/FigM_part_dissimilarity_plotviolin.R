####
#### R script for Ohigashi et al (2024)
#### violin plots for total dissimilarity, replacement, and richness difference
#### 2024.09.12 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(tidyr); packageVersion("tidyr")
library(reshape2); packageVersion("reshape2")
library(ggplot2); packageVersion("ggplot2")
library(ggsignif); packageVersion("ggsignif")
library(cowplot); packageVersion("cowplot")
source("Function/F1_HelperFunctions.R")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data (created in 09_partitioning_dissimilarity.R)
# prokaryotes dissimilarity
prok_df <- read.table("08_partitioning_dissimilarity_out/paritionied_dissimilarity_prok.txt", header = T, sep = "\t")
prok_df <- prok_df |>
  filter(!landuse %in% c("Natural_Farm", "Farm_Natural")) |> # remove pair containing different land uses.
  mutate(Dissim = 1-Similarity) |> # put dissimilarity column
  select(Dissim, Repl, RichDiff, landuse, relation)

# fungi dissimilarity
fungi_df <- read.table("08_partitioning_dissimilarity_out/paritionied_dissimilarity_fungi.txt", header = T, sep = "\t")
fungi_df <- fungi_df |>
  filter(!landuse %in% c("Natural_Farm", "Farm_Natural")) |> # remove pair containing different land uses.
  mutate(Dissim = 1-Similarity) |> # put dissimilarity column
  select(Dissim, Repl, RichDiff, landuse, relation)

# p values obtained from permutation test
prok_p <- read.table("08_partitioning_dissimilarity_out/part_dissim_permpval_prok.txt", header = T, sep = "\t")
fungi_p <- read.table("08_partitioning_dissimilarity_out/part_dissim_permpval_fungi.txt", header = T, sep = "\t")

# change p values to asterisks
prok_p <- prok_p |>
  mutate_all(~ case_when(
    . < 0.001 ~ "***",
    . < 0.01  ~ "**",
    . < 0.05  ~ "*",
    TRUE      ~ "n.s."
  ))

fungi_p <- fungi_p |>
  mutate_all(~ case_when(
    . < 0.001 ~ "***",
    . < 0.01  ~ "**",
    . < 0.05  ~ "*",
    TRUE      ~ "n.s."
  ))



### looping for creating plots
scale_cat <- c("within_site", "diff_site_within_country", "diff_site_diff_country")
rownames(prok_p) <- scale_cat
rownames(fungi_p) <- scale_cat

# prokaryotes
prok_plots <- list()
for (i in scale_cat) {
  # subset
  subdf <- prok_df |> filter(relation == i)
  dissim_max <- max(subdf$Dissim)*1.2
  repl_max <- max(subdf$Repl)*1.2
  rich_max <- max(subdf$RichDiff)*1.2
  
  # melt for ggplot2
  meltdf <- melt(subdf, id.vars = c("landuse", "relation"))
  
  # format category names
  meltdf <- meltdf |>
    mutate(variable = case_when(
      variable == "Dissim" ~ "Total\ndissimilarity",
      variable == "Repl" ~ "Replacement",
      variable == "RichDiff" ~ "Richness\ndifference",
      TRUE ~ NA
    ))
  meltdf$variable <- factor(meltdf$variable, levels = c("Total\ndissimilarity", "Replacement", "Richness\ndifference"))
  meltdf$landuse <- sub("_", " vs ", meltdf$landuse)
  meltdf$landuse <- factor(meltdf$landuse, levels = c("Natural vs Natural", "Farm vs Farm"))
  
  # set plot titles
  if (i == "within_site") {
    plottitle <- "Prokaryotes (within-site)"
  } else if (i == "diff_site_within_country") {
    plottitle <- "Prokaryotes (across-site, within the same country)"
  } else {
    plottitle <- "Prokaryotes (across-site, between the countries)"
  }
  
  # get p values
  anno <- unlist(c(prok_p[i, ]))
  
  # plotting
  p <- ggplot(meltdf, aes_string(x = "variable", y = "value", fill = "landuse")) +
    geom_violin(trim = F) +
    geom_boxplot(width=0.2, position = position_dodge(width = 0.9)) +
    labs(
      title = plottitle,
      y=NULL,
      x=NULL,
      fill = "Pair type"
      )+
    theme_classic()+
    scale_fill_manual(values = c("palegreen", "tan1"))+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.title = element_text(size = 12, face = "bold", colour = "black"),
      legend.text = element_text(size = 12, face = "bold", colour = "black"),
      plot.title = element_text(size=13, face = "bold", colour = "black"),
      axis.text=element_text(size=12, face = "bold", colour = "black"),
      axis.title=element_text(size=15,face="bold", colour = "black"),
      legend.position = "right",
      strip.text.x = element_text(face="bold", size = 12),
      strip.background = element_rect(fill="white", colour="black", size=1)
    ) +
    annotate("text", #dummy
             x=-Inf,
             y = dissim_max,
             label="",
             hjust="left",
             vjust=1,
             size = 3) +
    geom_signif(
      annotations = anno,
      y_position=c(dissim_max*0.95, repl_max*0.95, rich_max*0.95),
      xmin=c(0.75, 1.75, 2.75),
      xmax=c(1.25, 2.25, 3.25),
      tip_length=0.01,
      textsize=7,
      size=1,
      color="black"
    )
  prok_plots[[i]] <- p
}


# fungi
fungi_plots <- list()
for (i in scale_cat) {
  # subset
  subdf <- fungi_df |> filter(relation == i)
  dissim_max <- max(subdf$Dissim)*1.2
  repl_max <- max(subdf$Repl)*1.2
  rich_max <- max(subdf$RichDiff)*1.2
  
  # melt for ggplot2
  meltdf <- melt(subdf, id.vars = c("landuse", "relation"))
  
  # format category names
  meltdf <- meltdf |>
    mutate(variable = case_when(
      variable == "Dissim" ~ "Total\ndissimilarity",
      variable == "Repl" ~ "Replacement",
      variable == "RichDiff" ~ "Richness\ndifference",
      TRUE ~ NA
    ))
  meltdf$variable <- factor(meltdf$variable, levels = c("Total\ndissimilarity", "Replacement", "Richness\ndifference"))
  meltdf$landuse <- sub("_", " vs ", meltdf$landuse)
  meltdf$landuse <- factor(meltdf$landuse, levels = c("Natural vs Natural", "Farm vs Farm"))
  
  # set plot titles
  if (i == "within_site") {
    plottitle <- "Fungi (within-site)"
  } else if (i == "diff_site_within_country") {
    plottitle <- "Fungi (across-site, within the same country)"
  } else {
    plottitle <- "Fungi (across-site, between the countries)"
  }
  
  # get p values
  anno <- unlist(c(fungi_p[i, ]))
  
  # plotting
  p <- ggplot(meltdf, aes_string(x = "variable", y = "value", fill = "landuse")) +
    geom_violin(trim = F) +
    geom_boxplot(width=0.2, position = position_dodge(width = 0.9)) +
    labs(
      title = plottitle,
      y=NULL,
      x=NULL,
      fill = "Pair type"
    )+
    theme_classic()+
    scale_fill_manual(values = c("palegreen", "tan1"))+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.title = element_text(size = 12, face = "bold", colour = "black"),
      legend.text = element_text(size = 12, face = "bold", colour = "black"),
      plot.title = element_text(size=13, face = "bold", colour = "black"),
      axis.text=element_text(size=12, face = "bold", colour = "black"),
      axis.title=element_text(size=15,face="bold", colour = "black"),
      legend.position = "right",
      strip.text.x = element_text(face="bold", size = 12),
      strip.background = element_rect(fill="white", colour="black", size=1)
    ) +
    annotate("text", #dummy
             x=-Inf,
             y = dissim_max,
             label="",
             hjust="left",
             vjust=1,
             size = 3) +
    geom_signif(
      annotations = anno,
      y_position=c(dissim_max*0.95, repl_max*0.95, rich_max*0.95),
      xmin=c(0.75, 1.75, 2.75),
      xmax=c(1.25, 2.25, 3.25),
      tip_length=0.01,
      textsize=7,
      size=1,
      color="black"
    )
  fungi_plots[[i]] <- p
}

### save data
# create a directory
dir.create("FigCode/FigM_part_dissimilarity_out")

# save object and image files
for (i in scale_cat) {
  # prokaryotes
  saveRDS(prok_plots[[i]], file = sprintf("FigCode/FigM_part_dissimilarity_out/%s_dissimprok.obj", i))
  ggsave(prok_plots[[i]], file = sprintf("FigCode/FigM_part_dissimilarity_out/%s_dissimprok.png", i),
         bg = "white", width = 8, height = 7)
  # fungi
  saveRDS(fungi_plots[[i]], file = sprintf("FigCode/FigM_part_dissimilarity_out/%s_dissimfungi.obj", i))
  ggsave(fungi_plots[[i]], file = sprintf("FigCode/FigM_part_dissimilarity_out/%s_dissimfungi.png", i),
         bg = "white", width = 8, height = 7)
}


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigM_part_dissimilarity_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


