####
#### R script for Ohigashi et al (2024)
#### Visualize variation partitioning of community composition and functional composition
#### 2024.11.04 written by Ohigashi
#### R 4.3.3
#### 


### load packages
library(tidyverse); packageVersion("tidyverse")
library(reshape2); packageVersion("reshape2")
library(ggh4x); packageVersion("ggh4x")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data
# create a dataframe for variation partitioning results
vpR2_all.df <- data.frame(X = c("EnvAlone", "SpaceAlone", "Overlap", "Residuals"))

# loop to get the results
for (i in c("prok", "fungi")) {
  for (j in c("com", "func")) {
    for (k in c("nat", "farm")) {
      dfname <- paste0(i, j, ".", k)
      tmp.df <- read.csv(sprintf("14_variation_partitioning_out/vp_R2_%s.csv", dfname))
      vpR2_all.df <- vpR2_all.df |>
        left_join(tmp.df |> select(X, Adj.R.squared), by = "X")
      colnames(vpR2_all.df)[ncol(vpR2_all.df)] <- paste0(dfname, ".r2.adj")
    }
  }
}


### format data
# transform data frame
vpR2_all.ml.df <- melt(vpR2_all.df, id.vars = "X")

# add columns for type and landuse
vpR2_all.ml.df <- vpR2_all.ml.df |>
  mutate(type = case_when(
    grepl("com", variable) ~ "Community",
    TRUE ~ "Function"
  ),
  Landuse = case_when(
    grepl("nat", variable) ~ "Natural",
    TRUE ~ "Farm"
  ),
  orgs = case_when(
    grepl("prok", variable) ~ "Prokaryotes",
    TRUE ~ "Fungi"
  )
  )

# convert negative adjusted R2 values to 0
vpR2_all.ml.df <- vpR2_all.ml.df %>%
  mutate(value = ifelse(value < 0, 0, value))

# adjust R2 values to percentage
plot.df <- vpR2_all.ml.df |>
  group_by(orgs, type, Landuse) |>
  summarise(value.percent = value/sum(value)*100,
            X = X) |>
  ungroup()

# make tidy version of X
plot.df <- plot.df|>
  mutate(X2 = case_when(
    X == "EnvAlone" ~ "Environment",
    X == "Overlap" ~ "Environment + Space",
    X == "SpaceAlone" ~ "Space",
    TRUE ~ "Residuals"
  ))


### create plot
# set color palette
colpal <- colors[1:length(unique(plot.df$X))]

# set levels of factors
plot.df$Landuse <- factor(plot.df$Landuse, levels = c("Natural", "Farm"))
plot.df$X2 <- factor(plot.df$X2, levels = c("Environment", "Space",
                                                          "Environment + Space", "Residuals"))
plot.df$orgs <- factor(plot.df$orgs, levels = c("Prokaryotes", "Fungi"))

# plot
p <- ggplot(plot.df, aes(x = Landuse, y = value.percent, fill = X2)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_nested(~ orgs + type, scales = "free") + # nested 
  scale_fill_manual(values = colpal) +
  theme_classic() +
  labs(x = NULL, y = "Variation Explained (%)", fill = "") +
  theme(
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )
p


### save data
dir.create("FigCode/FigM_VariationPartitioning_out")

# plot object
ggsave("FigCode/FigM_VariationPartitioning_out/variation_explained.png", plot = p, width = 11.5, height = 7, bg = "white")
saveRDS(p, "FigCode/FigM_VariationPartitioning_out/variation_explained_plot.obj")

# used data frame
saveRDS(plot.df, "FigCode/FigM_VariationPartitioning_out/variation_explained_df.obj")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigM_VariationPartitioning_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

