####
#### R script for Ohigashi et al (2024)
#### model selection for multiple regression for heterogeneity of within-site and across-site scale
#### 2025.04.17 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(broom); packageVersion("broom")
library(car); packageVersion("car")
library(MuMIn); packageVersion("MuMIn")
library(stringr); packageVersion("stringr")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data
# distance to centroids
dists_within <- read.table("05_distance_to_centroid_out/distance_to_centroids_withinsite.txt", header = T)
dists_within <- dists_within |> rownames_to_column(var = "Sample")
dists_across <- read.table("05_distance_to_centroid_out/distance_to_centroids_acrosssite.txt", header = T, sep = "\t")
dists_across <- dists_across |> rownames_to_column(var = "Sample")

# environment
env_data <- read.table("Data/soil_metadata.txt", header = T)
colnames(env_data)[5] <- "Moisture"


### calculation of standard deviation of environmental factors
# use only soil physicochemical factors
env_data <- env_data |> dplyr::select(-log_copy16S_persoil, -log_copyITS_persoil, -shannon_16S, -shannon_ITS)

# add a column for location (site x land use)
env_data <- env_data |> mutate(treatment = paste(Site, Landuse, sep = "_"))
env_data_tre <- env_data |> dplyr::select(treatment, Moisture, Carbon, Nitrogen, CN_ratio, pH)

# calculate absolute values of deviation within the location (e.g., Site D-Natural)
env_dev_within <- env_data_tre |>
  group_by(treatment) |>
  mutate(across(c(Moisture, Carbon, Nitrogen, pH), 
                list(deviation_abs = ~ abs(. - mean(.))), 
                .names = "{col}_absdev_within"))

# add a column for sample
row.names(env_dev_within) <- env_data$Sample
env_dev_within <- env_dev_within |> rownames_to_column(var = "Sample")

# calculate absolute values of deviation across sites (e.g., Natural)
env_data_lu <- env_data |> dplyr::select(Sample, Landuse, Moisture, Carbon, Nitrogen, CN_ratio, pH)
env_dev_across <- env_data_lu |>
  group_by(Landuse) |>
  mutate(across(c(Moisture, Carbon, Nitrogen, pH), 
                list(deviation_abs = ~ abs(. - mean(.))), 
                .names = "{col}_absdev_across"))


### multiple regression
## 1. within-site scale
data_within <- env_dev_within
data_within$prokdist_within <- dists_within$proktaxa_within # add dist to cent of prokaryotes
data_within$fungidist_within <- dists_within$fungitaxa_within # add dist to cent of fungi
data_within <- data_within |>
  mutate(Farming = ifelse(grepl("Natural", treatment), 0, 1))

## 1-1. prokaryotes (within-site)
# create model for prokaryotes
model_prok_within <- lm(scale(prokdist_within) ~ scale(Farming) + scale(Moisture) + scale(Carbon) + scale(Nitrogen) + scale(pH) +
                          scale(Moisture_absdev_within) + scale(Carbon_absdev_within) + scale(Nitrogen_absdev_within) + scale(pH_absdev_within),
                        data = data_within, na.action = na.fail)

# calculate BIC using all combinations of variables
model_prok_within_dre <- dredge(model_prok_within, rank = "BIC")
# extract models with delta BIC < 2 => top 5 model for checking purpose
model_prok_within_topsdf <- model_prok_within_dre[1:5, ]
# format the data frame
prok_within_tops <- model_prok_within_topsdf %>%
  as.data.frame() %>%
  mutate(Model = paste0("Prok_within_", seq(1:nrow(.)))) %>%
  select(Model, BIC, delta, `scale(Farming)`, `scale(Moisture)`, `scale(Carbon)`, `scale(Nitrogen)`,
        `scale(pH)`, `scale(Moisture_absdev_within)`, `scale(Carbon_absdev_within)`, `scale(Nitrogen_absdev_within)`,
        `scale(pH_absdev_within)`)
varnames <- colnames(prok_within_tops) 
varnames <- str_remove_all(varnames, "scale\\(|\\)")
varnames <- str_replace(varnames, "_absdev_(within|across)", " (dev)")
colnames(prok_within_tops) <- varnames

# round values and replace minus values and NA with appropriate characters
prok_within_tops <- prok_within_tops %>%
  mutate(across(
    2:ncol(.),
    ~ if_else(is.na(.x), "", scaleFUN2(.x))
  ))


## 1-2. fungi (within-site)
# create model for fungi
model_fungi_within <- lm(scale(fungidist_within) ~ scale(Farming) + scale(Moisture) + scale(Carbon) + scale(Nitrogen) + scale(pH) +
                          scale(Moisture_absdev_within) + scale(Carbon_absdev_within) + scale(Nitrogen_absdev_within) + scale(pH_absdev_within),
                        data = data_within, na.action = na.fail)

# calculate BIC using all combinations of variables
model_fungi_within_dre <- dredge(model_fungi_within, rank = "BIC")
# extract top 5 models
model_fungi_within_topsdf <- model_fungi_within_dre[1:5,]
# format the data frame
fungi_within_tops <- model_fungi_within_topsdf %>%
  as.data.frame() %>%
  mutate(Model = paste0("Fungi_within_", seq(1:nrow(.)))) %>%
  select(Model, BIC, delta, `scale(Farming)`, `scale(Moisture)`, `scale(Carbon)`, `scale(Nitrogen)`,
         `scale(pH)`, `scale(Moisture_absdev_within)`, `scale(Carbon_absdev_within)`, `scale(Nitrogen_absdev_within)`,
         `scale(pH_absdev_within)`)
colnames(fungi_within_tops) <- varnames

# round values and replace minus values and NA with appropriate characters
fungi_within_tops <- fungi_within_tops %>%
  mutate(across(
    2:ncol(.),
    ~ if_else(is.na(.x), "", scaleFUN2(.x))
  ))


## 2. across-site scale
data_across <- env_dev_across
data_across$prokdist_across <- dists_across$proktaxa_across # add dist to cent of prokaryotes
data_across$fungidist_across <- dists_across$fungitaxa_across # add dist to cent of fungi
data_across <- data_across |>
  mutate(Farming = ifelse(grepl("Natural", Landuse), 0, 1))

## 2-1. prokaryotes (across-site)
# create model for prokaryotes
model_prok_across <- lm(scale(prokdist_across) ~ scale(Farming) + scale(Moisture) + scale(Carbon) + scale(Nitrogen) + scale(pH) +
                          scale(Moisture_absdev_across) + scale(Carbon_absdev_across) + scale(Nitrogen_absdev_across) + scale(pH_absdev_across),
                        data = data_across, na.action = na.fail)

# calculate AICc using all combinations of variables
model_prok_across_dre <- dredge(model_prok_across, rank = "BIC")
model_prok_across_topsdf <- model_prok_across_dre[1:5,]
prok_across_tops <- model_prok_across_topsdf %>%
  as.data.frame() %>%
  mutate(Model = paste0("Prok_across_", seq(1:nrow(.)))) %>%
  select(Model, BIC, delta, `scale(Farming)`, `scale(Moisture)`, `scale(Carbon)`, `scale(Nitrogen)`,
         `scale(pH)`, `scale(Moisture_absdev_across)`, `scale(Carbon_absdev_across)`, `scale(Nitrogen_absdev_across)`,
         `scale(pH_absdev_across)`)
colnames(prok_across_tops) <- varnames

# round values and replace minus values and NA with appropriate characters
prok_across_tops <- prok_across_tops %>%
  mutate(across(
    2:ncol(.),
    ~ if_else(is.na(.x), "", scaleFUN2(.x))
  ))


## 2-2. fungi (across-site)
# create model for prokaryotes
model_fungi_across <- lm(scale(fungidist_across) ~ scale(Farming) + scale(Moisture) + scale(Carbon) + scale(Nitrogen) + scale(pH) +
                          scale(Moisture_absdev_across) + scale(Carbon_absdev_across) + scale(Nitrogen_absdev_across) + scale(pH_absdev_across),
                        data = data_across, na.action = na.fail)

# calculate AICc using all combinations of variables
model_fungi_across_dre <- dredge(model_fungi_across, rank = "BIC")
model_fungi_across_topsdf <- model_fungi_across_dre[1:5,]

fungi_across_tops <- model_fungi_across_topsdf %>%
  as.data.frame() %>%
  mutate(Model = paste0("Fungi_across_", seq(1:nrow(.)))) %>%
  select(Model, BIC, delta, `scale(Farming)`, `scale(Moisture)`, `scale(Carbon)`, `scale(Nitrogen)`,
         `scale(pH)`, `scale(Moisture_absdev_across)`, `scale(Carbon_absdev_across)`, `scale(Nitrogen_absdev_across)`,
         `scale(pH_absdev_across)`)
colnames(fungi_across_tops) <- varnames

# round values and replace minus values and NA with appropriate characters
fungi_across_tops <- fungi_across_tops %>%
  mutate(across(
    2:ncol(.),
    ~ if_else(is.na(.x), "", scaleFUN2(.x))
  ))


### save data
write.csv(prok_within_tops, "13_multiple_regression_for_heterogeneity_out/top5bic_models_prok_within.csv",
          row.names = F, quote = F)
write.csv(fungi_within_tops, "13_multiple_regression_for_heterogeneity_out/top5bic_models_fungi_within.csv",
          row.names = F, quote = F)
write.csv(prok_across_tops, "13_multiple_regression_for_heterogeneity_out/top5bic_models_prok_across.csv",
          row.names = F, quote = F)
write.csv(fungi_across_tops, "13_multiple_regression_for_heterogeneity_out/top5bic_models_fungi_across.csv",
          row.names = F, quote = F)

### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/13_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
