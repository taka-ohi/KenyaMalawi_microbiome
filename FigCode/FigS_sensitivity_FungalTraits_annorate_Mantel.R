####
#### R script for Ohigashi et al (2024)
#### check relationship between FungalTraits annotation rates and main results (Mantel correlogram)
#### 2025.06.12 written by Ohigashi
#### R 4.3.3
####



### load packages and functions
# source("Function/F2_HelperFunctions_for_Visualization.R")
# source("Function/F1_HelperFunctions.R")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(tibble); packageVersion("tibble")
library(purrr); packageVersion("purrr")
library(SoDA); packageVersion("SoDA")



### load data
# Fungal traits annotation rate
anno_rates <- read.csv("FigCode/FigS_FungalTraits_summary_out/assignedASV_percent.csv", header = T)

# environment data
env_data <- read.table("Data/soil_metadata.txt", header = T)
env_category <- env_data[,1:4] # extract category data
landuse <- env_data$Landuse
site <- env_data$Site
# with "plots"
env_data2 <- env_data |>
  mutate(Plot=rep(rep(c("1", "2", "3"), each = 3), times = 10)) # add Plot column for later use

# gps
gpsdata <- read.csv("Data/soil_gpsdata.csv", header = T)

# fungal function data
fungilife <- read.table("02_Function_analysis_out/FungalTraits_primarilifestyle_percent.txt", header = T)
fungilife <- fungilife |> filter(primary_lifestyle != "unassigned") # 2024.09.11
fungilife <- fungilife |> column_to_rownames(var = "primary_lifestyle")
fungilife.t <- t(fungilife)

# prokaryotic function data
prok_func_all <- read.table("02_Function_analysis_out/PICRUSt2_full_function.txt", header = T)

# set function converting Lat/Lon data (e.g., N32° -> 32)
dms_to_deg <- function(dms_string) {
  parts <- as.numeric(strsplit(gsub("[^0-9.]", " ", dms_string), " ")[[1]])
  deg <- parts[2] + parts[3]/60 + parts[4]/3600
  if (substring(dms_string, 1, 1) %in% c("W", "S")) {
    deg <- -deg
  }
  return(deg)
}

## convert Lat/Lon data to XY data
# convert North-South, East-West data to values
gpsdata <- gpsdata |>
  mutate(lon_val = map_dbl(Coord_long, dms_to_deg),
         lat_val = map_dbl(Coord_lat, dms_to_deg))

# get XY coordinate from longitude and latitude values
xy <- geoXY(latitude = gpsdata$lat_val, longitude = gpsdata$lon_val)
xy <- as.data.frame(xy) # convert to data frame

# format XY data (add plot, site, country, and landuse columns)
xy <- xy |>
  mutate(Plot = rep(rep(c("1", "2", "3")), times = 10), 
         Country = gpsdata$Country,
         Site = gpsdata$Site,
         Landuse = gpsdata$Landuse
  )

# combine environmental data and xy data
env_xy <- merge(env_data2, xy, by = c("Country", "Site", "Landuse", "Plot"), all.x = T, sort = F)

# add replication column to set different XY coordinates for the replications
env_xy <- env_xy |> 
  mutate(Reps = rep(rep(c("1", "2", "3")), times = 30))

# add a tiny value to set different XY coordinates for the replications 
env_xy <- env_xy |>
  group_by(Plot) |>
  mutate(X = case_when(
    Reps == 1 ~ X + 0,
    Reps == 2 ~ X + 0.3,
    Reps == 3 ~ X + 0.6,
    TRUE ~ X),
    Y = case_when(
      Reps == 1 ~ Y + 0,
      Reps == 2 ~ Y + 0.3,
      Reps == 3 ~ Y + 0.6,
      TRUE ~ Y)
  ) |>
  ungroup()


## for Mantel correlogram, exclude samples with some of the lowest samples within group (e.g., Site D-Farm)
# define a function to run Mantel correlogram while removing some samples
compute_mantel_correlog_grp <- function(mat,
                                        n_remove = 1,
                                        nperm = 9999,
                                        random = FALSE) {
  library(dplyr)
  library(vegan)
  
  # define extracted samples
  remain_samples2 <- anno_rates %>%
    left_join(env_category, by = "Sample") %>%
    group_by(Site, Landuse) %>%
    { 
      if (random) {
        slice_sample(., prop = 1) %>%
          slice((n_remove + 1):n())
      } else {
        arrange(., assigned_percent) %>%
          slice((n_remove + 1):n())
      }
    } %>%
    ungroup() %>%
    pull(Sample)
  
  # filtering coordinates
  env_xy_filt <- env_xy %>% filter(Sample %in% remain_samples2)
  env_xy_nat <- env_xy_filt %>%
    filter(Landuse == "Natural") %>%
    select(X, Y)
  env_xy_farm <- env_xy_filt %>%
    filter(Landuse == "Farm") %>%
    select(X, Y)
  
  # extract lifestyle data
  mat_filt <- mat[rownames(mat) %in% remain_samples2, ]
  mat_lu <- as.data.frame(mat_filt)
  mat_lu$Landuse <- env_xy_filt$Landuse
  
  mat_nat <- mat_lu %>%
    filter(Landuse == "Natural") %>%
    select(-Landuse)
  mat_farm <- mat_lu %>%
    filter(Landuse == "Farm") %>%
    select(-Landuse)
  
  # calculate Bray-Curtis
  mat.bray_nat <- vegdist(mat_nat, method = "bray")
  mat.bray_farm <- vegdist(mat_farm, method = "bray")
  
  # Mantel correlogram（Natural）
  set.seed(123)
  mat.correlog.bray_nat <- mantel.correlog(
    mat.bray_nat,
    XY = env_xy_nat,
    nperm = nperm,
    cutoff = FALSE,
    break.pts = c(0, 200, 30000, 70000, 1444000, 1500000)
  )
  
  # Mantel correlogram（Farm）
  mat.correlog.bray_farm <- mantel.correlog(
    mat.bray_farm,
    XY = env_xy_farm,
    nperm = nperm,
    cutoff = FALSE,
    break.pts = c(0, 200, 30000, 70000, 1444000, 1500000)
  )
  
  # format data
  df_nat <- as.data.frame(mat.correlog.bray_nat$mantel.res)
  df_nat$landuse <- "Natural"
  df_farm <- as.data.frame(mat.correlog.bray_farm$mantel.res)
  df_farm$landuse <- "Farm"
  
  # class.index
  class.index <- rep(c("0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000"), 2)
  
  # combine
  df_combined <- bind_rows(df_nat, df_farm)
  rownames(df_combined) <- NULL
  df_combined$class.index <- class.index
  
  return(df_combined)
}




##### compute Mantel when randomly removing samples from each group #####
library(pbapply)

# randomly remove 2 samples from each group 
set.seed(123)
mantel_rm_random2_grp <- pblapply(seq_along(1:1000), function(i) {
  # set different rng.seed for each iteration
  set.seed(100 + i)
  # compute mantel correlogram
  compute_mantel_correlog_grp(mat = fungilife.t, n_remove = 2, nperm = 999, random = TRUE)
}, cl = 32)
# make a data frame for the result
mantel_rm_random2_grp_df <- bind_rows(mantel_rm_random2_grp, .id = "replicate")

# remove 23 lowest samples based on FungalTraits annotation rate
mantel_rm_low2_grp_df <- compute_mantel_correlog_grp(mat = fungilife.t, n_remove = 2, nperm = 999, random = FALSE)

# randomly remove 1 sample from each
set.seed(123)
mantel_rm_random1_grp <- pblapply(seq_along(1:1000), function(i) {
  # set different rng.seed for each iteration
  set.seed(100 + i)
  # compute mantel correlogram
  compute_mantel_correlog_grp(mat = fungilife.t, n_remove = 1, nperm = 999, random = TRUE)
}, cl = 32)
# make a data frame for the result
mantel_rm_random1_grp_df <- bind_rows(mantel_rm_random1_grp, .id = "replicate")

# remove 1 lowest samples based on FungalTraits annotation rate
mantel_rm_low1_grp_df <- compute_mantel_correlog_grp(mat = fungilife.t, n_remove = 1, nperm = 999, random = FALSE)




####### confirmation with PICRUSt2 data #######
# randomly remove 2 sample from each group
set.seed(123)
prok_mantel_rm_random2_grp <- pblapply(seq_along(1:1000), function(i) {
  # set different rng.seed for each iteration
  set.seed(100 + i)
  # compute mantel correlogram
  compute_mantel_correlog_grp(mat = prok_func_all, n_remove = 2, nperm = 999, random = TRUE)
}, cl = 32)
# make a data frame for the result
prok_mantel_rm_random2_grp_df <- bind_rows(prok_mantel_rm_random2_grp, .id = "replicate")

# remove 2 lowest samples based on FungalTraits annotation rate
prok_mantel_rm_low2_grp_df <- compute_mantel_correlog_grp(mat = prok_func_all, n_remove = 2, nperm = 999, random = FALSE)


# randomly remove 1 sample from each group
set.seed(123)
prok_mantel_rm_random1_grp <- pblapply(seq_along(1:1000), function(i) {
  # set different rng.seed for each iteration
  set.seed(100 + i)
  # compute mantel correlogram
  compute_mantel_correlog_grp(mat = prok_func_all, n_remove = 1, nperm = 999, random = TRUE)
}, cl = 32)
# make a data frame for the result
prok_mantel_rm_random1_grp_df <- bind_rows(prok_mantel_rm_random1_grp, .id = "replicate")

# remove 1 lowest samples based on FungalTraits annotation rate
prok_mantel_rm_low1_grp_df <- compute_mantel_correlog_grp(mat = prok_func_all, n_remove = 1, nperm = 999, random = FALSE)




###### create plots ######
# format data frame for 95% CI
mantel2_ci_df <- mantel_rm_random2_grp_df %>%
  group_by(class.index, landuse) %>%
  summarise(
    mean_cor = mean(Mantel.cor, na.rm = TRUE),
    ci_lower = quantile(Mantel.cor, probs = 0.025, na.rm = TRUE),
    ci_upper = quantile(Mantel.cor, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

mantel1_ci_df <- mantel_rm_random1_grp_df %>%
  group_by(class.index, landuse) %>%
  summarise(
    mean_cor = mean(Mantel.cor, na.rm = TRUE),
    ci_lower = quantile(Mantel.cor, probs = 0.025, na.rm = TRUE),
    ci_upper = quantile(Mantel.cor, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )


prok_mantel2_ci_df <- prok_mantel_rm_random2_grp_df %>%
  group_by(class.index, landuse) %>%
  summarise(
    mean_cor = mean(Mantel.cor, na.rm = TRUE),
    ci_lower = quantile(Mantel.cor, probs = 0.025, na.rm = TRUE),
    ci_upper = quantile(Mantel.cor, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

prok_mantel1_ci_df <- prok_mantel_rm_random1_grp_df %>%
  group_by(class.index, landuse) %>%
  summarise(
    mean_cor = mean(Mantel.cor, na.rm = TRUE),
    ci_lower = quantile(Mantel.cor, probs = 0.025, na.rm = TRUE),
    ci_upper = quantile(Mantel.cor, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )



# store data frames as a list
analysis_list <- list("fungi2" = list("low" = mantel_rm_low2_grp_df, "random" = mantel2_ci_df),
                      "fungi1" = list("low" = mantel_rm_low1_grp_df, "random" = mantel1_ci_df),
                      "prok2" = list("low" = prok_mantel_rm_low2_grp_df, "random" = prok_mantel2_ci_df),
                      "prok1" = list("low" = prok_mantel_rm_low1_grp_df, "random" = prok_mantel1_ci_df)
)

# create a list to store plots
plot_list <- list()

# function to round the results to show with one decimal place
scaleFUN <- function(x) {
  a <- sprintf("%.1f", x)
  a2 <- sub('^-', '\U2212', format(a))
  a2 <- trimws(a2)
  return(a2)
}

# loop
for (run in names(analysis_list)) {
  # get list of data frames
  ana_list <- analysis_list[[run]]
  # get data frames
  low_df <- ana_list$low
  random_ci_ori <- ana_list$random
  
  ## format data
  low_df <- low_df %>%
    select(class.index, landuse, Mantel.cor) %>% 
    mutate(type = "Low annotation removal")
  random_ci <- random_ci_ori %>%
    select(class.index, landuse, Mantel.cor = mean_cor) %>%
    mutate(type = "Random removal mean")
  # combine
  plot_df <- bind_rows(low_df, random_ci)
  # set factor levels
  plot_df$landuse <- factor(plot_df$landuse, levels = c("Natural", "Farm"))
  random_ci_ori$landuse <- factor(random_ci_ori$landuse, levels = c("Natural", "Farm"))
  
  # define class labels
  plot_df <- plot_df %>%
    mutate(dist.class = case_when(
      class.index == "0_200" ~ "1",
      class.index == "200_30000" ~ "2",
      class.index == "30000_70000" ~ "3",
      class.index == "70000_1444000" ~ "4",
      class.index == "1444000_1500000" ~ "5",
    ))
  random_ci_ori <- random_ci_ori %>%
    mutate(dist.class = case_when(
      class.index == "0_200" ~ "1",
      class.index == "200_30000" ~ "2",
      class.index == "30000_70000" ~ "3",
      class.index == "70000_1444000" ~ "4",
      class.index == "1444000_1500000" ~ "5",
    ))
  
  
  # set title names
  if (grepl("fungi", run)) {
    titlepart <- sprintf("Fungal lifestyles; \U2212%s sample per group", sub("fungi", "", run))
  } else {
    titlepart <- sprintf("Prokaryotic functions; \U2212%s sample per group", sub("prok", "", run))
  }
  
  # create Mantel correlogram
  p_man <- ggplot() +
    # 95% confidence interval
    geom_errorbar(
      data = random_ci_ori,
      aes(x = dist.class, ymin = ci_lower, ymax = ci_upper, color = landuse,  group = landuse),
      width = 0.2,
      position = position_dodge(width = 0.4),
      show.legend = FALSE
    ) +
    # removed low annotation rate samples
    geom_point(
      data = plot_df,
      aes(x = dist.class, y = Mantel.cor, color = landuse, shape = type,  group = landuse),
      size = 3,
      position = position_dodge(width = 0.4),
      show.legend = TRUE
    ) +
    scale_shape_manual(values = c(16, 2)) +
    scale_color_manual(values = c("Farm"="tan1", "Natural"="darkgreen")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_classic() +
    labs(x = "Distance class", y = "Mantel's r",
         title = sprintf("Correlation by distance class\n(%s)", titlepart),
         color = "Land use", shape = "Sample removal method") +
    scale_y_continuous(labels = scaleFUN) +
    theme(
      plot.title = element_text(color = "black", hjust = 0.5),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 12, face = "bold", colour = "black"),
      legend.text = element_text(size = 12, colour = "black")
    )
  
  plot_list[[run]] <- p_man
}




### save data
saveRDS(analysis_list, "FigCode/FigS_FungalTraits_summary_out/sensitivity_Mantel_dflist.rds")
saveRDS(plot_list, "FigCode/FigS_FungalTraits_summary_out/sensitivity_Mantel_plotlist.rds")


