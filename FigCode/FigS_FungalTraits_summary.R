####
#### R script for Ohigashi et al (2024)
#### Proportion of ASVs that was annotated with primary lifestyles
#### 2025.04.18 written by Ohigashi
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

# summarize by primary lifestyle
rel_primarylife <- rel_fungilife_ASV |>
  group_by(primary_lifestyle) |>
  summarise(across(1:90, sum))
# replace NA with "Unassigned"
rel_primarylife <- rel_primarylife %>% 
  mutate(primary_lifestyle = case_when(
    is.na(primary_lifestyle) ~ "Unassigned",
    TRUE ~ str_to_title(sub("_", " ", primary_lifestyle))
  ))

# convert lifestyle names to "Others" if it doesn't have 0.1% of abundance
rel_primarylife2 <- rel_primarylife %>%
  mutate(primary_lifestyle = ifelse(rowMeans(select(., 2:91)) < 1, "Others (< 1%)", primary_lifestyle))
# summarize by Genus again
rel_primarylife2 <- rel_primarylife2 |>
  group_by(primary_lifestyle) |>
  summarise(across(1:90, sum))


# primary lifestyle -> rownames, and transpose
rel_primarylife2 <- rel_primarylife2 |> column_to_rownames("primary_lifestyle")
rel_primarylife.t <- t(rel_primarylife2)
rel_primarylife.t <- as.data.frame(rel_primarylife.t)

# rownames -> sample
rel_primarylife.t <- rel_primarylife.t |> rownames_to_column("Sample")

# combine with envdata
rel_primarylife.df <- merge(envdata[,1:4], rel_primarylife.t, by = "Sample") # just used categorical environmental data

# melt
rel_primarylife.df <- melt(rel_primarylife.df, id.vars = c("Sample", "Country", "Site", "Landuse"))


### plot barchart
# set order of land use
rel_primarylife.df$Landuse <- factor(rel_primarylife.df$Landuse, levels = c("Natural", "Farm"))

# set order of lifestyles
rel_primarylife.df$variable <- factor(rel_primarylife.df$variable, levels = rev(c("Soil Saprotroph", "Litter Saprotroph",
                                                                                  "Wood Saprotroph", "Unspecified Saprotroph",
                                                                                  "Plant Pathogen", "Sooty Mold",
                                                                                  "Others (< 1%)", "Unassigned")))
# rel_primarylife.df$variable <- factor(rel_primarylife.df$variable, levels = c("Soil Saprotroph", "Litter Saprotroph",
#                                                                                   "Wood Saprotroph", "Unspecified Saprotroph",
#                                                                                   "Plant Pathogen", "Sooty Mold",
#                                                                                   "Others (< 1%)", "Unassigned"))


# set colors for plot
plo_cols <- colors
plo_cols[(length(unique(rel_primarylife.df$variable))-1):length(unique(rel_primarylife.df$variable))] <- c("grey30", "grey")

# plot
plife <- ggplot(rel_primarylife.df, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_nested(~ Site + Landuse, scales = "free") + # nested 
  scale_fill_manual(values = plo_cols[length(unique(rel_primarylife.df$variable)):1]) +
  theme_classic() +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Lifestyle") +
  theme(
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 14),
    axis.text.x = element_blank()
  )+
  guides(fill = guide_legend(reverse = TRUE))


### calculate percentage of annotate ASVs
assigned_percent <- rel_primarylife.t %>%
  mutate(assigned_percent = 100 - Unassigned) %>% 
  select(Sample, assigned_percent)

assigned_summary <- assigned_percent %>%
  summarise(max_assign = max(assigned_percent, na.rm = T),
            min_assign = min(assigned_percent, na.rm = T),
            mean_assign = mean(assigned_percent, na.rm = T)
            )

### save data
dir.create("FigCode/FigS_FungalTraits_summary")

# assigned percentage (raw)
write.csv(assigned_percent, "FigCode/FigS_FungalTraits_summary/assignedASV_percent.csv",
          row.names = F, quote = F)
# assigned percentage (summary)
write.csv(assigned_summary, "FigCode/FigS_FungalTraits_summary/assignedASV_percent_summary.csv",
          row.names = F, quote = F)

# figure
# save PDF
cairo_pdf("FigCode/FigS_FungalTraits_summary/primary_lifestyles.pdf", width = 12, height = 6.5)
print(plife)
dev.off()
ggsave("FigCode/FigS_FungalTraits_summary/primary_lifestyles.png", plot = plife, width = 12, height = 6.5, bg = "white")


