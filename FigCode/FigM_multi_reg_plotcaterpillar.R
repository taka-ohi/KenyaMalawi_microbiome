####
#### R script for Ohigashi et al (2024)
#### caterpillar plot with the multiple regression for microbial heterogeneity
#### 2024.09.15 written by Ohigashi
#### R 4.3.3
####


### load package and function
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(stringr); packageVersion("stringr")
library(cowplot); packageVersion("cowplot")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data
CI_prok_ws <- read.table("13_multiple_regression_for_heterogeneity_out/CI_prok_within.txt", header = T)
CI_fungi_ws <- read.table("13_multiple_regression_for_heterogeneity_out/CI_fungi_within.txt", header = T)
CI_prok_as <- read.table("13_multiple_regression_for_heterogeneity_out/CI_prok_across.txt", header = T)
CI_fungi_as <- read.table("13_multiple_regression_for_heterogeneity_out/CI_fungi_across.txt", header = T)

# make a list for later use
df_list <- list(CI_prok_ws = CI_prok_ws,
                CI_fungi_ws = CI_fungi_ws,
                CI_prok_as = CI_prok_as,
                CI_fungi_as = CI_fungi_as)


### create plot
plots <- list()
for (i in names(df_list)) {
  df <- df_list[[i]]
    
  # remove "scale(" and ")"
  df <- df |> 
    mutate(term = str_remove_all(term, "scale\\(|\\)")) |>
    mutate(term = str_replace(term, "_absdev_(within|across)", " (dev)"))
  
  df <- df |> 
    transform(term = factor(term, levels = c("pH (dev)", "Nitrogen (dev)", "Carbon (dev)", "Moisture (dev)",
                                           "pH", "Nitrogen", "Carbon", "Moisture", "Farming")))
  # set plot title
  plottitle <- "Regression for heterogeneity of\n"
  if (grepl("prok", i)) {
    plottitle <- paste0(plottitle, "prokaryotes (")
  } else {
    plottitle <- paste0(plottitle, "fungi (")
  }
  if (grepl("ws", i)) {
    plottitle <- paste0(plottitle, "within-site)")
  } else if (grepl("as", i)) {
    plottitle <- paste0(plottitle, "across-site)")
  }
  
  plot_CI <- 
    ggplot(data = df) +
    geom_line(aes(x = term,
                  y = estimate),
              color = "black",
              size = 3) +
    geom_pointrange(aes(x = term, 
                        y = estimate,
                        ymin = lower,
                        ymax = upper),
                    color = "black",
                    size = 0.75) +
    geom_hline(yintercept = 0, 
               linetype = 2, 
               color = "red") +
    geom_text(aes(x = term,
                  y = estimate,
                  label = estimate_signif
    ),
    vjust = -1) +
    labs(x = NULL,
         y = "Estimates",
         title = plottitle
    ) +
    theme_classic()+
    theme(axis.text = element_text(size = 12, color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", colour = "black")
    )+
    # scale_y_continuous(labels = scaleFUN2) +
    scale_y_continuous(breaks = seq(-1, 1, by=0.5), labels = scaleFUN2, limits = c(-1,1)) +
    coord_flip()
  
  plots[[i]] <- plot_CI
}


### create a panel for the main article
fig_all <- plot_grid(plots[[1]] + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                     plots[[2]] + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                     plots[[3]] + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                     plots[[4]] + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                     ncol = 2,
                     rel_widths = c(1, 1),
                     labels = c("a", "b", "c", "d"),
                     label_size = 20
)


### save data
# create a directory
dir.create("FigCode/FigM_multi_reg_out")

# save object and image files
for (i in names(df_list)) {
  saveRDS(plots[[i]], file = sprintf("FigCode/FigM_multi_reg_out/multireg%s.obj", i))
  ggsave(plots[[i]], file = sprintf("FigCode/FigM_multi_reg_out/multireg%s.png", i),
         bg = "white", width = 7, height = 7)
}


# save PDF for the main article
cairo_pdf("FigCode/FigM_multi_reg_out/FigM_multireg.pdf", width = 12, height = 12)
print(fig_all)
dev.off()

# save png
ggsave(filename = "FigCode/FigM_multi_reg_out/FigM_multireg.png",
       plot = fig_all, width = 12, height = 12, bg = "white")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigM_multi_reg_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))







