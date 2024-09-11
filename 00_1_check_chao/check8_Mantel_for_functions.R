

### load packages and functions
source("Function/F1_HelperFunctions.R")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(tibble)
library(purrr); packageVersion("purrr")
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

# prokaryotic functions
b_func.table <- read.table("02_Function_analysis_out/PICRUSt2_full_function.txt", header = T) # prokaryotes
# b_func.t <- t(b_func.table) # transpose

# fungal functions
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


### Mantel correlogram
## extract Natural and Farm communities
# prokaryotic functions
b_func.t_lu <- as.data.frame(b_func.table)
b_func.t_lu$Landuse <- env_data$Landuse
b_func.t_nat <- b_func.t_lu |>
  dplyr::filter(Landuse == "Natural") |> # extract
  dplyr::select(-Landuse)
b_func.t_farm <- b_func.t_lu |>
  dplyr::filter(Landuse == "Farm") |> # extract
  dplyr::select(-Landuse)

# fungi
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
  cutoff = FALSE, # set as FALSE otherwise it generates NA
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(b_func.correlog.bray_nat)
prok_mantel_nat <- b_func.correlog.bray_nat$mantel.res

# farm
b_func.correlog.bray_farm <- mantel.correlog(
  b_func.bray_farm,
  XY = env_xy_farm[15:16],
  nperm = 99999,
  cutoff = FALSE, # set as FALSE otherwise it generates NA
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(b_func.correlog.bray_farm)
prok_mantel_farm <- b_func.correlog.bray_farm$mantel.res

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
f_func.correlog.bray_nat <- mantel.correlog(
  f_func.bray_nat,
  XY = env_xy_nat[15:16],
  nperm = 99999,
  cutoff = FALSE, # set as FALSE otherwise it generates NA
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(f_func.correlog.bray_nat)
fungi_mantel_nat <- f_func.correlog.bray_nat$mantel.res

# farm
f_func.correlog.bray_farm <- mantel.correlog(
  f_func.bray_farm,
  XY = env_xy_farm[15:16],
  nperm = 99999,
  cutoff = FALSE, # set as FALSE otherwise it generates NA
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(f_func.correlog.bray_farm)
fungi_mantel_farm <- f_func.correlog.bray_farm$mantel.res

# combine Natural and Farm data
fungi_mantel_nat <- as.data.frame(fungi_mantel_nat)
fungi_mantel_nat$landuse <- rep("Natural", 5)
fungi_mantel_farm <- as.data.frame(fungi_mantel_farm)
fungi_mantel_farm$landuse <- rep("Farm", 5)

fungi_mantel_df <- rbind(fungi_mantel_nat, fungi_mantel_farm)
fungi_mantel_df$class.index <- c("0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000",
                                 "0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000")

### save data
write.table(prok_mantel_df, "00_1_check_chao/mantel_result_prokfuncall.txt", row.names = F, quote = F, sep = "\t")
write.table(fungi_mantel_df, "00_1_check_chao/mantel_result_fungilife.txt", row.names = F, quote = F, sep = "\t")


