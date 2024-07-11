####
#### R script for Ohigashi et al (2024)
#### multiple regression for heterogeneity of within-site and across-site scale
#### 2024.07.11 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(broom); packageVersion("broom")
library(car); packageVersion("car")
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
                        data = data_within)

# check multicollinearity
vif_prok_within <- vif(model_prok_within)
print(vif_prok_within)

# select variables by step-wise method
model_prok_within_step <- step(model_prok_within, direction = "both")
summary(model_prok_within_step)

# calculate 95% confidence interval
CI_prok_within <- model_prok_within |> # for all variables
  tidy() |>
  mutate(lower = estimate - qnorm(0.975) * std.error,
         upper = estimate + qnorm(0.975) * std.error) |>
  filter(term != "(Intercept)") |>
  mutate(estimate_signif = ifelse(p.value < 0.001, paste(scaleFUN2(estimate), "***", sep = ""),
                                  ifelse(p.value < 0.01, paste(scaleFUN2(estimate), "**", sep = ""),
                                         ifelse(p.value < 0.05, paste(scaleFUN2(estimate), "*", sep = ""),
                                                scaleFUN2(estimate)))))

CI_prok_within_step <- model_prok_within_step |> # selected model
  tidy() |>
  mutate(lower = estimate - qnorm(0.975) * std.error,
         upper = estimate + qnorm(0.975) * std.error) |>
  filter(term != "(Intercept)") |>
  mutate(estimate_signif = ifelse(p.value < 0.001, paste(scaleFUN2(estimate), "***", sep = ""),
                                  ifelse(p.value < 0.01, paste(scaleFUN2(estimate), "**", sep = ""),
                                         ifelse(p.value < 0.05, paste(scaleFUN2(estimate), "*", sep = ""),
                                                scaleFUN2(estimate)))))

## 1-2. fungi (within-site)
# create model for fungi
model_fungi_within <- lm(scale(fungidist_within) ~ scale(Farming) + scale(Moisture) + scale(Carbon) + scale(Nitrogen) + scale(pH) +
                          scale(Moisture_absdev_within) + scale(Carbon_absdev_within) + scale(Nitrogen_absdev_within) + scale(pH_absdev_within),
                        data = data_within)

# check multicollinearity
vif_fungi_within <- vif(model_fungi_within)
print(vif_fungi_within)

# select variables by step-wise method
model_fungi_within_step <- step(model_fungi_within, direction = "both")
summary(model_fungi_within_step)

# calculate 95% confidence interval
CI_fungi_within <- model_fungi_within |>
  tidy() |>
  mutate(lower = estimate - qnorm(0.975) * std.error,
         upper = estimate + qnorm(0.975) * std.error) |>
  filter(term != "(Intercept)") |>
  mutate(estimate_signif = ifelse(p.value < 0.001, paste(scaleFUN2(estimate), "***", sep = ""),
                                  ifelse(p.value < 0.01, paste(scaleFUN2(estimate), "**", sep = ""),
                                         ifelse(p.value < 0.05, paste(scaleFUN2(estimate), "*", sep = ""),
                                                scaleFUN2(estimate)))))

CI_fungi_within_step <- model_fungi_within_step |>
  tidy() |> 
  mutate(lower = estimate - qnorm(0.975) * std.error,
         upper = estimate + qnorm(0.975) * std.error) |>
  filter(term != "(Intercept)") |>
  mutate(estimate_signif = ifelse(p.value < 0.001, paste(scaleFUN2(estimate), "***", sep = ""),
                                  ifelse(p.value < 0.01, paste(scaleFUN2(estimate), "**", sep = ""),
                                         ifelse(p.value < 0.05, paste(scaleFUN2(estimate), "*", sep = ""),
                                                scaleFUN2(estimate)))))


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
                        data = data_across)

# check multicollinearity
vif_prok_across <- vif(model_prok_across)
print(vif_prok_across)

# select variables by step-wise method
model_prok_across_step <- step(model_prok_across, direction = "both")
summary(model_prok_across_step)

# calculate 95% confidence interval
CI_prok_across <- model_prok_across |> # for all variables
  tidy() |>
  mutate(lower = estimate - qnorm(0.975) * std.error,
         upper = estimate + qnorm(0.975) * std.error) |>
  filter(term != "(Intercept)") |>
  mutate(estimate_signif = ifelse(p.value < 0.001, paste(scaleFUN2(estimate), "***", sep = ""),
                                  ifelse(p.value < 0.01, paste(scaleFUN2(estimate), "**", sep = ""),
                                         ifelse(p.value < 0.05, paste(scaleFUN2(estimate), "*", sep = ""),
                                                scaleFUN2(estimate)))))

CI_prok_across_step <- model_prok_across_step |> # selected model
  tidy() |>
  mutate(lower = estimate - qnorm(0.975) * std.error,
         upper = estimate + qnorm(0.975) * std.error) |>
  filter(term != "(Intercept)") |>
  mutate(estimate_signif = ifelse(p.value < 0.001, paste(scaleFUN2(estimate), "***", sep = ""),
                                  ifelse(p.value < 0.01, paste(scaleFUN2(estimate), "**", sep = ""),
                                         ifelse(p.value < 0.05, paste(scaleFUN2(estimate), "*", sep = ""),
                                                scaleFUN2(estimate)))))

## 2-2. fungi (across-site)
# create model for prokaryotes
model_fungi_across <- lm(scale(fungidist_across) ~ scale(Farming) + scale(Moisture) + scale(Carbon) + scale(Nitrogen) + scale(pH) +
                          scale(Moisture_absdev_across) + scale(Carbon_absdev_across) + scale(Nitrogen_absdev_across) + scale(pH_absdev_across),
                        data = data_across)

# check multicollinearity
vif_fungi_across <- vif(model_fungi_across)
print(vif_fungi_across)

# select variables by step-wise method
model_fungi_across_step <- step(model_fungi_across, direction = "both")
summary(model_fungi_across_step)

# calculate 95% confidence interval
CI_fungi_across <- model_fungi_across |> # for all variables
  tidy() |>
  mutate(lower = estimate - qnorm(0.975) * std.error,
         upper = estimate + qnorm(0.975) * std.error) |>
  filter(term != "(Intercept)") |>
  mutate(estimate_signif = ifelse(p.value < 0.001, paste(scaleFUN2(estimate), "***", sep = ""),
                                  ifelse(p.value < 0.01, paste(scaleFUN2(estimate), "**", sep = ""),
                                         ifelse(p.value < 0.05, paste(scaleFUN2(estimate), "*", sep = ""),
                                                scaleFUN2(estimate)))))

CI_fungi_across_step <- model_fungi_across_step |> # selected model
  tidy() |>
  mutate(lower = estimate - qnorm(0.975) * std.error,
         upper = estimate + qnorm(0.975) * std.error) |>
  filter(term != "(Intercept)") |>
  mutate(estimate_signif = ifelse(p.value < 0.001, paste(scaleFUN2(estimate), "***", sep = ""),
                                  ifelse(p.value < 0.01, paste(scaleFUN2(estimate), "**", sep = ""),
                                         ifelse(p.value < 0.05, paste(scaleFUN2(estimate), "*", sep = ""),
                                                scaleFUN2(estimate)))))


### save data
dir.create("09_multiple_regression_for_heterogeneity_out")

## within-site regression
# prokaryotes
saveRDS(model_prok_within, "09_multiple_regression_for_heterogeneity_out/model_prok_within.obj")
write.table(CI_prok_within, "09_multiple_regression_for_heterogeneity_out/CI_prok_within.txt", row.names = F, quote = F, sep = "\t")
saveRDS(model_prok_within_step, "09_multiple_regression_for_heterogeneity_out/model_prok_within_selected.obj")
write.table(CI_prok_within_step, "09_multiple_regression_for_heterogeneity_out/CI_prok_within_selected.txt", row.names = F, quote = F, sep = "\t")

# fungi
saveRDS(model_fungi_within, "09_multiple_regression_for_heterogeneity_out/model_fungi_within.obj")
write.table(CI_fungi_within, "09_multiple_regression_for_heterogeneity_out/CI_fungi_within.txt", row.names = F, quote = F, sep = "\t")
saveRDS(model_fungi_within_step, "09_multiple_regression_for_heterogeneity_out/model_fungi_within_selected.obj")
write.table(CI_fungi_within_step, "09_multiple_regression_for_heterogeneity_out/CI_fungi_within_selected.txt", row.names = F, quote = F, sep = "\t")

## across-site regression
# prokaryotes
saveRDS(model_prok_across, "09_multiple_regression_for_heterogeneity_out/model_prok_across.obj")
write.table(CI_prok_across, "09_multiple_regression_for_heterogeneity_out/CI_prok_across.txt", row.names = F, quote = F, sep = "\t")
saveRDS(model_prok_across_step, "09_multiple_regression_for_heterogeneity_out/model_prok_across_selected.obj")
write.table(CI_prok_across_step, "09_multiple_regression_for_heterogeneity_out/CI_prok_across_selected.txt", row.names = F, quote = F, sep = "\t")

# fungi
saveRDS(model_fungi_across, "09_multiple_regression_for_heterogeneity_out/model_fungi_across.obj")
write.table(CI_fungi_across, "09_multiple_regression_for_heterogeneity_out/CI_fungi_across.txt", row.names = F, quote = F, sep = "\t")
saveRDS(model_fungi_across_step, "09_multiple_regression_for_heterogeneity_out/model_fungi_across_selected.obj")
write.table(CI_fungi_across_step, "09_multiple_regression_for_heterogeneity_out/CI_fungi_across_selected.txt", row.names = F, quote = F, sep = "\t")

# collinearity
write.table(vif_prok_within, "09_multiple_regression_for_heterogeneity_out/vif_within.txt", row.names = T, quote = F, sep = "\t") # same as vif_fungi_within
write.table(vif_prok_across, "09_multiple_regression_for_heterogeneity_out/vif_across.txt", row.names = T, quote = F, sep = "\t") # same as vif_fungi_across


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/09_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

