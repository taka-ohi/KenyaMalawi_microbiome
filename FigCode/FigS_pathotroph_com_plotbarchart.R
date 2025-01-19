####
#### R script for Ohigashi et al (2024)
#### barchart for fungal pathogen composition
#### 2024.11.30 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(ggh4x); packageVersion("ggh4x")
library(reshape2); packageVersion("reshape2")
library(tibble); packageVersion("tibble")
library(stringr); packageVersion("stringr")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data
# fungal ASV table with lifestyles
fungilife_w_ASVtbl <- read.table("02_Function_analysis_out/FungalTraits_w_rarefied_ASV_table_fungi.txt", header = T, sep = "\t")

# environmental data
envdata <- read.table("Data/soil_metadata.txt", header = T, sep = "\t")

### format data
# convert the value to relative abundance (percent)
rel_fungilife_ASV <- fungilife_w_ASVtbl %>%
  mutate(across(1:90, ~ . / sum(.) * 100))

# extract only pathotroph (i.e., "patho" in primary lifestyle)
rel_patho_ASV <- rel_fungilife_ASV %>%
  filter(grepl("patho", primary_lifestyle))

## create data frames before plotting
# 1. data frame by primary lifestyle
# summarize by primary lifestyle
rel_patho.life <- rel_patho_ASV |>
  group_by(primary_lifestyle) |>
  summarise(across(1:90, sum))

# primary lifestyle -> rownames, and transpose
rel_patho.life <- rel_patho.life |> column_to_rownames("primary_lifestyle")
rel_patho.life.t <- t(rel_patho.life)
rel_patho.life.t <- as.data.frame(rel_patho.life.t)

# rownames -> sample
rel_patho.life.t <- rel_patho.life.t |> rownames_to_column("Sample")

# combine with envdata
rel_patho.life.df <- merge(envdata[,1:4], rel_patho.life.t, by = "Sample") # just used categorical environmental data

# melt
rel_patho.life.df <- melt(rel_patho.life.df, id.vars = c("Sample", "Country", "Site", "Landuse"))

# 2. data frame by Genus
# summarize by Genus
rel_patho.genus <- rel_patho_ASV |>
  group_by(Genus) |>
  summarise(across(1:90, sum))

# convert genus names to "Others" if it doesn't have 0.1% of abundance
rel_patho.genus2 <- rel_patho.genus %>%
  mutate(Genus = ifelse(rowMeans(select(., 2:91)) < 0.1, "Others (< 0.1%)", Genus))
# summarize by Genus again
rel_patho.genus2 <- rel_patho.genus2 |>
  group_by(Genus) |>
  summarise(across(1:90, sum))

# genus -> rownames, and transpose
rel_patho.G <- rel_patho.genus2 |> column_to_rownames("Genus")
rel_patho.G.t <- t(rel_patho.G)
rel_patho.G.t <- as.data.frame(rel_patho.G.t)

# rownames -> sample
rel_patho.G.t <- rel_patho.G.t |> rownames_to_column("Sample")

# combine with envdata
rel_patho.G.df <- merge(envdata[,1:4], rel_patho.G.t, by = "Sample") # just used categorical environmental data

# melt
rel_patho.G.df <- melt(rel_patho.G.df, id.vars = c("Sample", "Country", "Site", "Landuse"))


### plotting
## 1. primary lifestyle
# capitalize the characters of lifestyle
rel_patho.life.df <- rel_patho.life.df |>
  mutate(variable = str_to_title(sub("_", " ", variable)))

# set order of land use
rel_patho.life.df$Landuse <- factor(rel_patho.life.df$Landuse, levels = c("Natural", "Farm"))

# plot
plife <- ggplot(rel_patho.life.df, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_nested(~ Site + Landuse, scales = "free") + # nested 
  scale_fill_manual(values = colors[1:length(unique(rel_patho.life.df$variable))]) +
  theme_classic() +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Lifestyle") +
  theme(
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 14),
    axis.text.x = element_blank()
  )


## 2. genus
# set order of land use
rel_patho.G.df$Landuse <- factor(rel_patho.G.df$Landuse, levels = c("Natural", "Farm"))

# set order of genus
# get a genus vector arranged by mean relative abundance
genus_means <- rel_patho.genus2 %>%
  mutate(mean_abundance = rowMeans(select(., 2:ncol(.)))) %>%
  arrange(desc(mean_abundance)) %>%
  pull(Genus)
# make "Others" at last
new_levels <- c(setdiff(genus_means, "Others (< 0.1%)"), "Others (< 0.1%)")
# set the levels of genus
rel_patho.G.df$variable <- factor(rel_patho.G.df$variable, levels = new_levels)


# set colors
colpal <- colors[1:length(unique(rel_patho.G.df$variable))]
colpal[length(colpal)] <- "grey" # set the last one (Others) to grey

# plot
pgenus <- ggplot(rel_patho.G.df, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_nested(~ Site + Landuse, scales = "free") + # nested 
  scale_fill_manual(values = colpal) +
  theme_classic() +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Genus") +
  theme(
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 14),
    axis.text.x = element_blank()
  )


### save data
dir.create("FigCode/FigS_pathotroph_subanalysis_out")

# lifestyle
ggsave("FigCode/FigS_pathotroph_subanalysis_out/relabun_lifestyle.png", plot = plife, width = 11.5, height = 7, bg = "white")
saveRDS(plife, "FigCode/FigS_pathotroph_subanalysis_out/relabun_lifestyle.obj")

# genus
ggsave("FigCode/FigS_pathotroph_subanalysis_out/relabun_genus.png", plot = pgenus, width = 11.5, height = 7, bg = "white")
saveRDS(pgenus, "FigCode/FigS_pathotroph_subanalysis_out/relabun_genus.obj")

# used data frames (genus one)
saveRDS(rel_patho.G.df, "FigCode/FigS_pathotroph_subanalysis_out/relabun_pathogen.genus_df.obj")

### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("FigCode/Fig_SessionInfo/FigS_pathotroph_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


