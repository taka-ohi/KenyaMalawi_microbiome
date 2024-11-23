####
#### R script for Ohigashi et al (2024)
#### Mantel correlogram for microbial community x spatial distance
#### 2024.07.12 written by Ohigashi -> 2024.07.14 modified by Ohigashi
#### 2024.09.11 edited by Ohigashi to add Mantel correlogram for functions
#### R 4.3.3
####


### load packages and functions
source("Function/F1_HelperFunctions.R")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(purrr); packageVersion("purrr")
library(tibble); packageVersion("tibble")
# install gfortran (using gfortran-12.2-universal.pkg) following the web site below
# https://cran.r-project.org/bin/macosx/tools/
# download SoDA package from the archive (https://cran.r-project.org/src/contrib/Archive/SoDA/)
# install.packages("~/Downloads/SoDA_1.0-6.1.tar.gz", type = "source")
library(SoDA); packageVersion("SoDA")


### load data
# gps
gpsdata <- read.csv("Data/soil_gpsdata.csv", header = T)

# environment
env_data <- read.table("Data/soil_metadata.txt", header = T)
env_data <- env_data |>
  mutate(Plot=rep(rep(c("1", "2", "3"), each = 3), times = 10)) # add Plot column for later use

# microbial counts
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T) # prokaryotes
b_ASV <- b_ASV.table[,1:(ncol(b_ASV.table)-7)] # remove taxa
b_ASV.t <- t(b_ASV) # transpose

f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T) # fungi
f_ASV <- f_ASV.table[,1:(ncol(f_ASV.table)-7)] # remove taxa
f_ASV.t <- t(f_ASV) # transpose

# microbial functions
b_func.table <- read.table("02_Function_analysis_out/PICRUSt2_full_function.txt", header = T) # prokaryotes

f_func.table <- read.table("02_Function_analysis_out/FungalTraits_primarilifestyle_percent.txt", header = T) # prokaryotes
f_func.table <- f_func.table |> filter(primary_lifestyle != "unassigned")
f_func.table <- f_func.table |> column_to_rownames(var = "primary_lifestyle")
f_func.t <- t(f_func.table) # transpose

### convert Lat/Lon data to XY data
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
env_xy <- merge(env_data, xy, by = c("Country", "Site", "Landuse", "Plot"), all.x = T, sort = F)

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


### Mantel correlogram for microbial communities
## extract Natural and Farm communities
# prokaryotes
b_ASV.t_lu <- as.data.frame(b_ASV.t)
b_ASV.t_lu$Landuse <- env_data$Landuse
b_ASV.t_nat <- b_ASV.t_lu |>
  dplyr::filter(Landuse == "Natural") |> # extract
  dplyr::select(-Landuse)
b_ASV.t_farm <- b_ASV.t_lu |>
  dplyr::filter(Landuse == "Farm") |> # extract
  dplyr::select(-Landuse)

# fungi
f_ASV.t_lu <- as.data.frame(f_ASV.t)
f_ASV.t_lu$Landuse <- env_data$Landuse
f_ASV.t_nat <- f_ASV.t_lu |>
  dplyr::filter(Landuse == "Natural") |> # extract
  dplyr::select(-Landuse)
f_ASV.t_farm <- f_ASV.t_lu |>
  dplyr::filter(Landuse == "Farm") |> # extract
  dplyr::select(-Landuse)

# calculate bray dissimilarity
b_ASV.bray_nat <- vegdist(b_ASV.t_nat, method = "bray")
b_ASV.bray_farm <- vegdist(b_ASV.t_farm, method = "bray")
f_ASV.bray_nat <- vegdist(f_ASV.t_nat, method = "bray")
f_ASV.bray_farm <- vegdist(f_ASV.t_farm, method = "bray")

# extract XY data of Natural and Farm samples
env_xy_nat <- env_xy |> dplyr::filter(Landuse == "Natural")
env_xy_farm <- env_xy |> dplyr::filter(Landuse == "Farm")

## Mantel test (prokaryotes)
set.seed(123)
# natural
b_ASV.correlog.bray_nat <- mantel.correlog(
  b_ASV.bray_nat,
  XY = env_xy_nat[15:16],
  nperm = 99999,
  cutoff = FALSE, # calculate for all classes
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(b_ASV.correlog.bray_nat)
prok_mantel_nat <- b_ASV.correlog.bray_nat$mantel.res

# farm
b_ASV.correlog.bray_farm <- mantel.correlog(
  b_ASV.bray_farm,
  XY = env_xy_farm[15:16],
  nperm = 99999,
  cutoff = FALSE, # calculate for all classes
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(b_ASV.correlog.bray_farm)
prok_mantel_farm <- b_ASV.correlog.bray_farm$mantel.res

# combine Natural and Farm data
prok_mantel_nat <- as.data.frame(prok_mantel_nat)
prok_mantel_nat$landuse <- rep("Natural", 5)
prok_mantel_farm <- as.data.frame(prok_mantel_farm)
prok_mantel_farm$landuse <- rep("Farm", 5)

prok_mantel_df <- rbind(prok_mantel_nat, prok_mantel_farm)
prok_mantel_df$class.index <- c("0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000",
                                "0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000")

## Mantel test (fungi)
# natural
f_ASV.correlog.bray_nat <- mantel.correlog(
  f_ASV.bray_nat,
  XY = env_xy_nat[15:16],
  nperm = 99999,
  cutoff = FALSE, # calculate for all classes
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(f_ASV.correlog.bray_nat)
fungi_mantel_nat <- f_ASV.correlog.bray_nat$mantel.res

# farm
f_ASV.correlog.bray_farm <- mantel.correlog(
  f_ASV.bray_farm,
  XY = env_xy_farm[15:16],
  nperm = 99999,
  cutoff = FALSE, # calculate for all classes
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(f_ASV.correlog.bray_farm)
fungi_mantel_farm <- f_ASV.correlog.bray_farm$mantel.res

# combine Natural and Farm data
fungi_mantel_nat <- as.data.frame(fungi_mantel_nat)
fungi_mantel_nat$landuse <- rep("Natural", 5)
fungi_mantel_farm <- as.data.frame(fungi_mantel_farm)
fungi_mantel_farm$landuse <- rep("Farm", 5)

fungi_mantel_df <- rbind(fungi_mantel_nat, fungi_mantel_farm)
fungi_mantel_df$class.index <- c("0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000",
                                 "0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000")



### Mantel correlogram for microbial functions (added on 2024.09.11)
## extract Natural and Farm functions
# prokaryotic functions
b_func.t_lu <- as.data.frame(b_func.table)
b_func.t_lu$Landuse <- env_data$Landuse
b_func.t_nat <- b_func.t_lu |>
  dplyr::filter(Landuse == "Natural") |> # extract
  dplyr::select(-Landuse)
b_func.t_farm <- b_func.t_lu |>
  dplyr::filter(Landuse == "Farm") |> # extract
  dplyr::select(-Landuse)

# fungal lifestyles
f_func.t_lu <- as.data.frame(f_func.t)
f_func.t_lu$Landuse <- env_data$Landuse
f_func.t_nat <- f_func.t_lu |>
  dplyr::filter(Landuse == "Natural") |> # extract
  dplyr::select(-Landuse)
f_func.t_farm <- f_func.t_lu |>
  dplyr::filter(Landuse == "Farm") |> # extract
  dplyr::select(-Landuse)

# calculate bray dissimilarity 
b_func.bray_nat <- vegdist(b_func.t_nat, method = "bray")
b_func.bray_farm <- vegdist(b_func.t_farm, method = "bray")
f_func.bray_nat <- vegdist(f_func.t_nat, method = "bray")
f_func.bray_farm <- vegdist(f_func.t_farm, method = "bray")

# extract XY data of Natural and Farm samples
env_xy_nat <- env_xy |> dplyr::filter(Landuse == "Natural")
env_xy_farm <- env_xy |> dplyr::filter(Landuse == "Farm")

## Mantel test (prokaryotic functions)
set.seed(123)
# natural
b_func.correlog.bray_nat <- mantel.correlog(
  b_func.bray_nat,
  XY = env_xy_nat[15:16],
  nperm = 99999,
  cutoff = FALSE, # calculate for all classes
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(b_func.correlog.bray_nat)
prokallfunc_mantel_nat <- b_func.correlog.bray_nat$mantel.res

# farm
b_func.correlog.bray_farm <- mantel.correlog(
  b_func.bray_farm,
  XY = env_xy_farm[15:16],
  nperm = 99999,
  cutoff = FALSE, # calculate for all classes
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(b_func.correlog.bray_farm)
prokallfunc_mantel_farm <- b_func.correlog.bray_farm$mantel.res


# combine Natural and Farm data
prokallfunc_mantel_nat <- as.data.frame(prokallfunc_mantel_nat)
prokallfunc_mantel_nat$landuse <- rep("Natural", 5)
prokallfunc_mantel_farm <- as.data.frame(prokallfunc_mantel_farm)
prokallfunc_mantel_farm$landuse <- rep("Farm", 5)

prokallfunc_mantel_df <- rbind(prokallfunc_mantel_nat, prokallfunc_mantel_farm)
prokallfunc_mantel_df$class.index <- c("0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000",
                                "0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000")

## Mantel test (fungal lifestyles)
# natural
f_func.correlog.bray_nat <- mantel.correlog(
  f_func.bray_nat,
  XY = env_xy_nat[15:16],
  nperm = 99999,
  cutoff = FALSE, # calculate for all classes
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(f_func.correlog.bray_nat)
fungilife_mantel_nat <- f_func.correlog.bray_nat$mantel.res

# farm
f_func.correlog.bray_farm <- mantel.correlog(
  f_func.bray_farm,
  XY = env_xy_farm[15:16],
  nperm = 99999,
  cutoff = FALSE, # calculate for all classes
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(f_func.correlog.bray_farm)
fungilife_mantel_farm <- f_func.correlog.bray_farm$mantel.res

# combine Natural and Farm data
fungilife_mantel_nat <- as.data.frame(fungilife_mantel_nat)
fungilife_mantel_nat$landuse <- rep("Natural", 5)
fungilife_mantel_farm <- as.data.frame(fungilife_mantel_farm)
fungilife_mantel_farm$landuse <- rep("Farm", 5)

fungilife_mantel_df <- rbind(fungilife_mantel_nat, fungilife_mantel_farm)
fungilife_mantel_df$class.index <- c("0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000",
                                 "0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000")




######## how to determine the distance class (break.pts) above ☆ #########
# get distance matrix from the XY coordinates of samples
xycoord <- data.frame(X = env_xy$X, Y = env_xy$Y)
row.names(xycoord) <- env_xy$Sample
dist_matrix <- dist(xycoord)

# get the combinations of samples
sample_names <- env_data$Sample
combinations <- expand.grid(sample_names, sample_names, stringsAsFactors = FALSE) # list all pairs
combinations <- combinations[combinations[,1] < combinations[,2],] # extract pairs (removed the replaced pairs)
comb <- apply(combinations, 1, function(x) paste(x, collapse = "_"))
comb <- sort(comb)

# create a data frame of sample-sample distances and pair categories
distdata <- as.data.frame(dist_matrix)
colnames(distdata)[1] <- "dist_m"
distdata$comb <- comb

# add columns for land use, site, country
distdata$landuse <- paste(env_data$Landuse[match(substr(distdata$comb, 1, 3), env_data$Sample)],
                             env_data$Landuse[match(substr(distdata$comb, 5, 7), env_data$Sample)],
                             sep = "_"
                          )
distdata$site <- paste(env_data$Site[match(substr(distdata$comb, 1, 3), env_data$Sample)],
                          env_data$Site[match(substr(distdata$comb, 5, 7), env_data$Sample)],
                          sep = "_"
                       )
distdata$country <- paste(env_data$Country[match(substr(distdata$comb, 1, 3), env_data$Sample)],
                             env_data$Country[match(substr(distdata$comb, 5, 7), env_data$Sample)],
                             sep = "_"
                          )

# summarize the min and max distances for each site-site and landuse-landuse pairs
distdata_summary <- distdata |>
  group_by(site, landuse, country) |>
  summarise(Min_dist = min(dist_m),
            Max_dist = max(dist_m)
            )

# remove different-landuse pairs (e.g., Natural_Farm)
distdata_summary <- distdata_summary |>
  rowwise() |>
  dplyr::filter(strsplit(landuse, "_")[[1]][1] == strsplit(landuse, "_")[[1]][2])

# sort by Max distances
distdata_sorted <- distdata_summary |>
  arrange(Max_dist)

# categorize the pairs and add a column for the category
distdata_sorted <- distdata_sorted |>
  rowwise() |>
  mutate(relation = case_when(
    strsplit(site, "_")[[1]][1] == strsplit(site, "_")[[1]][2] ~ "within_site",  # if the 1st and 3rd characters are the same.
    strsplit(site, "_")[[1]][1] != strsplit(site, "_")[[1]][2] && strsplit(country, "_")[[1]][1] == strsplit(country, "_")[[1]][2] ~ "diff_site_within_country",
    TRUE ~ "diff_site_diff_country"  # else
  )) |>
  ungroup()

# the Max distances (m) were as follows
distdata_sorted |>
  group_by(relation) |>
  summarise(Max = max(Max_dist))
# relation                      Max
# <chr>                       <dbl>
# 1 diff_site_diff_country   1466717.
# 2 diff_site_within_country   65572.
# 3 within_site                  133.

# -> so, I set the distance class at 0~200 m, 200~70000 m, and 70000~1500000 m
# -> on 2024.07.14, increased the number of distance class 0~200, 200~30000, 30000~700000, 70000~1444000, 1444000~1500000
# 30000 is a distance that includes any pairs of 2 sites but without Site F.
# 1444000 is a distance that includes "Site F to G or H", but not "Site D or E to G or H"


### save data
dir.create("10_Mantel_correlogram_out")
write.table(prok_mantel_df, "10_Mantel_correlogram_out/mantel_result_prok.txt", row.names = F, quote = F, sep = "\t")
write.table(fungi_mantel_df, "10_Mantel_correlogram_out/mantel_result_fungi.txt", row.names = F, quote = F, sep = "\t")
write.table(prokallfunc_mantel_df, "10_Mantel_correlogram_out/mantel_result_prokallfunc.txt", row.names = F, quote = F, sep = "\t")
write.table(fungilife_mantel_df, "10_Mantel_correlogram_out/mantel_result_fungilife.txt", row.names = F, quote = F, sep = "\t")
write.table(distdata_sorted, "10_Mantel_correlogram_out/distance_category.txt", row.names = F, quote = F, sep = "\t")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/10_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

