####
#### R script for Ohigashi et al (2024)
#### showing correlations between fungal pathogenic genera & pH
#### 2024.11.30 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(tibble); packageVersion("tibble")
library(stringr); packageVersion("stringr")
library(kableExtra); packageVersion("kableExtra")
library(gt); packageVersion("gt")
library(gtsummary); packageVersion("gtsummary")
library(flextable); packageVersion("flextable")
source("Function/F1_HelperFunctions.R")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data
# reuse the dataframe that was used in the relative abundance of pathogenic genera
rel_patho.genus <- readRDS("FigCode/FigS_pathotroph_subanalysis_out/relabun_pathogen.genus_df.obj")

# environmental data
envdata <- read.table("Data/soil_metadata.txt", header = T, sep = "\t")


### format data
# convert dataframe to calculate "all pathotroph"
rel_patho.genus_raw <- dcast(rel_patho.genus, variable ~ Sample)

# variable -> rownames
rel_patho.genus_raw <- rel_patho.genus_raw |>
  column_to_rownames("variable") |>
  t() |>
  as.data.frame()

# add "all pathotroph" column
rel_patho.genus_raw <- rel_patho.genus_raw |>
  mutate(All = rowSums(rel_patho.genus_raw))

# remove Others
rel_patho.genus_raw <- rel_patho.genus_raw |> select(-`Others (< 0.1%)`)

# rownames -> Sample
rel_patho.genus_raw <- rel_patho.genus_raw |> rownames_to_column("Sample")

# melt
rel_patho.genus2 <- melt(rel_patho.genus_raw, id.vars = c("Sample"))

# merge with environemntal data
rel_patho.genus_env <- rel_patho.genus2 %>%
  left_join(envdata %>%
              select(Sample, Site, Landuse, pH)
              , by = "Sample")


### calculation of correlations for each genus
cor_res <- CorrInCat(rel_patho.genus_env, var1 = "pH", var2 = "value",
                     method = "pearson", category = "variable")



### create tidy table
# round r and convert p values
cor_res.tidy <- cor_res |>
  mutate(pval = case_when(
    p.value < 0.001 ~ "< 0.001",
    p.value < 0.01 ~ "< 0.01",
    p.value < 0.05 ~ "< 0.05",
    TRUE ~ "n.s."
  ),
  r = sub("-", "\U2212", round(r, 2))
  )

# select columns
cor_res.tidy <- cor_res.tidy |> select(r, "p" = pval)

### save
# save as pdf
cor_res.tidy %>%
  kbl(caption = "<b>Table S1. Pearson's correlation between soil pH and relative abundances of pathotrophic genera</b>", align = "c") %>%
  kable_classic(full_width = T, html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%   # bold the header
  row_spec(13, italic = FALSE) %>%
  column_spec(1, width = "4cm", italic = T) %>%
  column_spec(2, width = "4cm") %>%
  column_spec(3, width = "4cm") %>%
  save_kable("FigCode/FigS_pathotroph_subanalysis_out/patho_cor_table.pdf", density = 300)


# cor_res <- cor_res |> rownames_to_column("Genus")
# 
# cor_anno <- "\n"
# for (j in 1:nrow(cor_res)) {
#   cor_anno0 <- paste0(cor_res[j, 1], ": ", "r = ", round(cor_res[j, 2], 2), ", ", cor_res[j, 4], "\n")
#   cor_anno0 <- sub("-", "\U2212", cor_anno0)
#   cor_anno <- paste0(cor_anno, cor_anno0)
# }
# 
# # make a list of genera that exhibit a correlation between the pairs
# sig_genera <- cor_res |>
#   filter(pval != "n.s.") |>
#   pull(Genus)
# 
# # set order of land use
# rel_patho.genus_env$Landuse <- factor(rel_patho.genus_env$Landuse, levels = c("Natural", "Farm"))
# 
# # plot
# pcor <- ggplot(rel_patho.genus_env,
#                aes(x = pH, y = value, color = variable, shape = Landuse)) +
#   geom_point(size=3.5) + 
#   scale_shape_manual(values = c(16, 3))+
#   labs(x = "pH",
#        y = "Relative abundance (%)",
#        color = "Genus",
#        shape = "Land use",
#        title = "pH × Pathogenic genus"
#   ) + 
#   geom_smooth(data = subset(rel_patho.genus_env, variable %in% sig_genera),
#               method = "lm", se = FALSE, #inherit.aes = FALSE,
#               aes(group = variable)) +
#   theme_classic()+
#   theme(axis.text = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 12, color = "black"),
#         legend.title = element_text(size = 13, face = "bold", colour = "black"),
#         legend.text = element_text(size = 13, face = "bold", colour = "black"),
#         plot.title = element_text(hjust = 0.5, face = "bold", colour = "black")
#   ) +
#   annotate("text",
#            x = max(rel_patho.genus_env$pH)*0.85,
#            y = max(rel_patho.genus_env$value)*0.85,
#            label = cor_anno,
#            hjust = "left",
#            size = 4
#   )



### save data
# ggsave("FigCode/FigS_pathotroph_subanalysis_out/relabun_genus_pH.png", plot = pcor, width = 9, height = 7, bg = "white")
# saveRDS(pcor, "FigCode/FigS_pathotroph_subanalysis_out/relabun_genus_pH.obj")

