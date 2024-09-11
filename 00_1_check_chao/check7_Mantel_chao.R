
####
#### R script for Ohigashi et al (2024)
#### Mantel correlogram for microbial community x spatial distance
#### 2024.07.12 written by Ohigashi -> 2024.07.14 modified by Ohigashi
#### R 4.3.3
####


### load packages and functions
source("Function/F1_HelperFunctions.R")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
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

# microbial counts
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T) # prokaryotes
b_ASV <- b_ASV.table[,1:(ncol(b_ASV.table)-7)] # remove taxa
b_ASV.t <- t(b_ASV) # transpose

f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T) # fungi
f_ASV <- f_ASV.table[,1:(ncol(f_ASV.table)-7)] # remove taxa
f_ASV.t <- t(f_ASV) # transpose

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

# calculate chao dissimilarity 
b_ASV.chao_nat <- vegdist(b_ASV.t_nat, method = "chao")
b_ASV.chao_farm <- vegdist(b_ASV.t_farm, method = "chao")
f_ASV.chao_nat <- vegdist(f_ASV.t_nat, method = "chao")
f_ASV.chao_farm <- vegdist(f_ASV.t_farm, method = "chao")

# extract XY data of Natural and Farm samples
env_xy_nat <- env_xy |> dplyr::filter(Landuse == "Natural")
env_xy_farm <- env_xy |> dplyr::filter(Landuse == "Farm")

## Mantel test (prokaryotes)
set.seed(123)
# natural
b_ASV.correlog.chao_nat <- mantel.correlog(
  b_ASV.chao_nat,
  XY = env_xy_nat[15:16],
  nperm = 99999,
  cutoff = FALSE, # set as FALSE otherwise it generates NA
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(b_ASV.correlog.chao_nat)
prok_mantel_nat <- b_ASV.correlog.chao_nat$mantel.res

# farm
b_ASV.correlog.chao_farm <- mantel.correlog(
  b_ASV.chao_farm,
  XY = env_xy_farm[15:16],
  nperm = 99999,
  cutoff = FALSE, # set as FALSE otherwise it generates NA
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(b_ASV.correlog.chao_farm)
prok_mantel_farm <- b_ASV.correlog.chao_farm$mantel.res

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
f_ASV.correlog.chao_nat <- mantel.correlog(
  f_ASV.chao_nat,
  XY = env_xy_nat[15:16],
  nperm = 99999,
  cutoff = FALSE, # set as FALSE otherwise it generates NA
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(f_ASV.correlog.chao_nat)
fungi_mantel_nat <- f_ASV.correlog.chao_nat$mantel.res

# farm
f_ASV.correlog.chao_farm <- mantel.correlog(
  f_ASV.chao_farm,
  XY = env_xy_farm[15:16],
  nperm = 99999,
  cutoff = FALSE, # set as FALSE otherwise it generates NA
  break.pts = c(0, 200, 30000, 70000, 1444000, 1500000) # how I determined the values is written later☆
)
plot(f_ASV.correlog.chao_farm)
fungi_mantel_farm <- f_ASV.correlog.chao_farm$mantel.res

# combine Natural and Farm data
fungi_mantel_nat <- as.data.frame(fungi_mantel_nat)
fungi_mantel_nat$landuse <- rep("Natural", 5)
fungi_mantel_farm <- as.data.frame(fungi_mantel_farm)
fungi_mantel_farm$landuse <- rep("Farm", 5)

fungi_mantel_df <- rbind(fungi_mantel_nat, fungi_mantel_farm)
fungi_mantel_df$class.index <- c("0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000",
                                 "0_200", "200_30000", "30000_70000", "70000_1444000", "1444000_1500000")

### save data
write.table(prok_mantel_df, "00_1_check_chao/mantel_result_prok_chao.txt", row.names = F, quote = F, sep = "\t")
write.table(fungi_mantel_df, "00_1_check_chao/mantel_result_fungi_chao.txt", row.names = F, quote = F, sep = "\t")
