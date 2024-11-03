####
#### R script for Ohigashi et al (2024)
#### Plotting the mean contributions of ASV to the dissimilarity in the scale group
#### 2024.09.26 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(stringr); packageVersion("stringr")
library(scales); packageVersion("scales")
source("Function/F1_HelperFunctions.R")
source("Function/F2_HelperFunctions_for_Visualization.R")

### load data
# list of data frames that contain average ASVs' contribution to dissimilarity in each scale
ave_dfs <- readRDS("11_SIMPER_out/ASV_ave.contrib_eachscale.obj")

# p-values obtained by permutation test
perm_p_df <- read.csv("11_SIMPER_out/lifestyle_contrib_pvals_comparing_landuse_mean.csv", header = T)


### format data
## combine Natural and Farm data for each scale
# make a vector for scales
scale_vec <- c("within_D", "within_E", "within_F", "within_G", "within_H",
               "within_Kenya", "within_Malawi", "as_ac")

# looping to make data frames for each scale
plotdfs <- list()
for (cat in scale_vec) {
  # find the list elements that match the pattern for Natural and Farm
  nat_index <- grep(paste0(cat, "_Nat"), names(ave_dfs))
  farm_index <- grep(paste0(cat, "_Farm"), names(ave_dfs))
  
  # extract the corresponding list elements
  df_nat <- ave_dfs[[nat_index]]
  df_farm <- ave_dfs[[farm_index]]
  
  # name the value column to combine the data frames later
  df_nat <- df_nat |>
    filter(mean_rel > 0) |>
    select(ASV, w_Natural=mean_contrib, primary_lifestyle)
  df_farm <- df_farm |>
    filter(mean_rel > 0) |>
    select(ASV, w_Farm=mean_contrib, primary_lifestyle)
  
  # combine the data frames
  plotdf0 <- merge(df_nat, df_farm, by = c("ASV", "primary_lifestyle"))
  
  # melt the combined data frame
  plotdf0.melt <- melt(plotdf0, id.vars = c("ASV", "primary_lifestyle"),
                       variable.name = "group", value.name = "mean_contrib")
  
  # modify the group name
  plotdf0.melt <- plotdf0.melt |>
    mutate(group = case_when(
      group == "w_Natural" ~ "Natural vs Natural",
      group == "w_Farm" ~ "Farm vs Farm",
      TRUE ~ group
    ))
  plotdf0.melt$group <- factor(plotdf0.melt$group, levels = c("Natural vs Natural", "Farm vs Farm"))
  
  plotdfs[[cat]] <- plotdf0.melt
}


## change p values to asterisks
# vectorize the changep0 function
changep0.v <- Vectorize(changep0)
# replace p values with asterisks
perm_p_df <- perm_p_df |>
  mutate(across(-1, changep0.v))


### Wilcoxon Mann Whitney test <- will use only permutation result
# make a data frame to store p-values (containing same column names as permutation test one)
# u_p_df <- perm_p_df
# u_p_df[, 2:9] <- NA # reset it
# 
# # put W-M-W test p-values into the data frame
# for (dfname in names(plotdfs)) {
#   df <- plotdfs[[dfname]]
#   for (ls in u_p_df$primary_lifestyle) {
#     lifedf <- df |> filter(primary_lifestyle == ls)
#     wmw_res <- wilcox.test(mean_contrib ~ group, data = lifedf)
#     u_p <- wmw_res$p.value
#     if (dfname == "as_ac") {
#       u_p_df$as_ac[u_p_df$primary_lifestyle == ls] <- u_p
#     } else if (grepl("^within_[A-Z]$", dfname)) {
#       tar_col <- sub("within_", "ws", dfname)
#       u_p_df[[tar_col]][u_p_df$primary_lifestyle == ls] <- u_p
#     } else if (grepl("^within_[A-Z]+", dfname)) {
#       tar_col <- sub("within_", "w", dfname)
#       u_p_df[[tar_col]][u_p_df$primary_lifestyle == ls] <- u_p
#     }
#   }
# }
# 
# # chnage p-values to asterisks
# u_p_df <- u_p_df |>
#   mutate(across(-1, changep0.v))

### create plots for each scale type
# pick up lifestyles to be plotted
plotted_lifestyle <- c("plant_pathogen", "animal_parasite", "dung_saprotroph", "litter_saprotroph",
                       "soil_saprotroph", "unspecified_saprotroph", "wood_saprotroph")

# extract p values of plotted lifestyles
plot_perm_p <- perm_p_df |>
  filter(primary_lifestyle %in% plotted_lifestyle) |>
  mutate(primary_lifestyle = str_to_title(sub("_", " ", primary_lifestyle)))
colnames(plot_perm_p)[2:ncol(plot_perm_p)] <- names(plotdfs) # name as data frames in plotdfs are named

# looping to create plots
plots <- list()
# l <- names(plotdfs)[1]
for (l in names(plotdfs)) {
  # get data set
  plotdf <- plotdfs[[l]]
  
  # make tiny lifestyles as "others"
  plotdf <- plotdf |>
    mutate(primary_lifestyle = case_when(
      primary_lifestyle %in% plotted_lifestyle ~ str_to_title(sub("_", " ", primary_lifestyle)),
      TRUE ~ "Others/Unassigned"
    ))
  plotdf <- plotdf |> filter(primary_lifestyle != "Others/Unassigned")
  
  # set plot title
  if (grepl("^within_[A-Z]$", l)) {
    # within site
    plottitle <- gsub("within_", "Within Site ", l)
    # across site within the country
  } else if (grepl("^within_[A-Z]+", l)) {
    plottitle <- gsub("within_", "Across sites within ", l)
    plottitle <- gsub("_", " ", plottitle)
    # across site between the countries
  } else if (startsWith(l, "as_ac")) {
    plottitle <- gsub("as_ac", "Across sites between countries", l)
  }
  
  # data frame for p values
  pval.tbl <- plot_perm_p |> select(primary_lifestyle, l)
  
  # calculate max values in each lifestyle
  max_lifes <- plotdf |>
    group_by(primary_lifestyle) |>
    summarise(lifeRoof = max(mean_contrib)*10)
  
  # make a data frame for annotation
  ann.tbl <- merge(pval.tbl, max_lifes, by = "primary_lifestyle")
  ann.tbl$x <- c(0.85, 1.85, 2.85, 3.85, 4.85, 5.85, 6.85)
  ann.tbl$xend <- c(1.15, 2.15, 3.15, 4.15, 5.15, 6.15, 7.15)
  
  # plot
  p0 <- ggplot(plotdf, aes(x = primary_lifestyle, y = mean_contrib,
                           fill = group)) +
    geom_violin(trim = FALSE, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.8), color = "grey90", alpha = 0.2) +
    scale_fill_manual(values = c("darkgreen", "tan1")) +
    theme_classic() +
    labs(
      x = NULL,
      y = "Mean contribution of ASVs to dissimilarity",
      fill = "Group",
      title = plottitle
    ) +
    theme(
      # plot.title = element_text(hjust = 0.5),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(size = 8, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right"
    ) +
    scale_y_continuous(
      labels = label_scientific(),  # Scientific notation for the y-axis labels
      trans = "log10"               # Log10 scale on y-axis
    )
  
  # annotate significance
  p_ann <- p0 +
    geom_text(
      data = ann.tbl,
      mapping = aes(x = primary_lifestyle,
                    y = lifeRoof,
                    # fill = group,
                    label = !!sym(l)
      ),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.75)
    ) +
    geom_segment(data = ann.tbl,
               aes(x = x, xend = xend, y = lifeRoof*0.6),
               inherit.aes = FALSE,
               linetype = "solid", color = "black")

  # put the plot in the list
  plots[[l]] <- p_ann
}


### save data
dir.create("FigCode/FigM_SIMPER_out/contrib_within_grp")
for (name in names(plots)) {
  p <- plots[[name]]
  saveRDS(p, file = sprintf("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_%s.obj", name))
  ggsave(filename = sprintf("FigCode/FigM_SIMPER_out/contrib_within_grp/contrib_plot_%s.png", name), 
         plot = p, width = 7, height = 6, bg = "white")
}

