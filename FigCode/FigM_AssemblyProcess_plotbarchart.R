####
#### R script for Ohigashi et al (2024)
#### barchart for assembly process
#### 2024.11.03 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(ggh4x); packageVersion("ggh4x")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data
w.bNTI.RCbray_16S <- read.csv("13_AssemblyProcess_out/summary_bNTI_RCbray_prok.csv", header = T)
w.bNTI.RCbray_ITS <- read.csv("13_AssemblyProcess_out/summary_bNTI_RCbray_fungi.csv", header = T)


### format data
## calculate ratio of the processes
# prokaryotes
prok_res.df <- w.bNTI.RCbray_16S |> 
  group_by(pair_scale, pair_landu, process) |> 
  summarise(count = n(), .groups = "drop") |> 
  group_by(pair_scale, pair_landu) |> 
  mutate(ratio = count / sum(count)) |> 
  ungroup() |> 
  select(pair_scale, pair_landu, process, ratio)

# fungi
fungi_res.df <- w.bNTI.RCbray_ITS |> 
  group_by(pair_scale, pair_landu, process) |> 
  summarise(count = n(), .groups = "drop") |> 
  group_by(pair_scale, pair_landu) |> 
  mutate(ratio = count / sum(count)) |> 
  ungroup() |> 
  select(pair_scale, pair_landu, process, ratio)

# add a prok or fungi column
prok_res.df <- prok_res.df |> mutate(orgs = "Prokaryotes")
fungi_res.df <- fungi_res.df |> mutate(orgs = "Fungi")

# bind vertically
plot.df <- rbind(prok_res.df, fungi_res.df)

# add columns for plot
plot.df <- plot.df |>
  mutate(pair_scale2 = case_when(
    pair_scale == "within_site" ~ "Within-site",
    pair_scale == "across_site_within_country" ~ "Across-site\n(within the country)",
    pair_scale == "across_site_btw_country" ~ "Across-site\n(between the countries)",
    TRUE ~ NA
  )
  )
plot.df$pair_scale2 <- factor(plot.df$pair_scale2, levels = c("Within-site",
                                                              "Across-site\n(within the country)",
                                                              "Across-site\n(between the countries)"
                                                              ))
plot.df$pair_landu <- factor(plot.df$pair_landu, levels = c("Natural", "Farm", "cross_landuse"))
plot.df$orgs <- factor(plot.df$orgs, levels = c("Prokaryotes", "Fungi"))


### create plot
# remove "cross_landuse"
plot.df.fin <- plot.df |> filter(pair_landu != "cross_landuse")

# set color palette
colpal <- colors[1:length(unique(plot.df.fin$process))]

# plot
p <- ggplot(plot.df.fin, aes(x = pair_landu, y = ratio, fill = process)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_nested(~ orgs + pair_scale2, scales = "free") + # nested 
  scale_fill_manual(values = colpal) +
  theme_classic() +
  labs(x = NULL, y = "Ratio", fill = "Assembly Process") +
  theme(
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )


### save data
dir.create("FigCode/FigM_AssemblyProcess_out")

# plot object
ggsave("FigCode/FigM_AssemblyProcess_out/assemblyprocess_ratio.png", plot = p, width = 11.5, height = 7, bg = "white")
saveRDS(p, "FigCode/FigM_AssemblyProcess_out/assemblyprocess_ratio_plot.obj")

# used data frame
saveRDS(plot.df.fin, "FigCode/FigM_AssemblyProcess_out/assemblyprocess_ratio_df.obj")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigM_AssemblyProcess_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

