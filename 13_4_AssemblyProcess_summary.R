####
#### R script for Ohigashi et al (2024)
#### Summarization of bNTI and RCbray
#### 2024.11.03 written by Ohigashi
#### R 4.3.3



### load package
library(dplyr); packageVersion("dplyr")
library(reshape2); packageVersion("reshape2")


### load data
# beta NTI (saved as data frame but the shape of distance matrix)
weighted.bNTI.16S <- read.csv("13_AssemblyProcess_out/weighted_bNTI_16S.csv", row.names = 1)
weighted.bNTI.ITS <- read.csv("13_AssemblyProcess_out/weighted_bNTI_ITS.csv", row.names = 1)
# convert them to dist object
w.bNTI.16S.d <- as.dist(weighted.bNTI.16S)
w.bNTI.ITS.d <- as.dist(weighted.bNTI.ITS)

# RCbray (dist object)
b_rc_res <- readRDS("13_AssemblyProcess_out/RCbray_prok.obj")
f_rc_res <- readRDS("13_AssemblyProcess_out/RCbray_fungi.obj")

# environmental data (country, site, land use)
env <- read.table("Data/soil_metadata.txt", header = T)


### format data
# transform dist to data frame
# bNTI
w.bNTI.16S.df <- melt(as.matrix(w.bNTI.16S.d), varnames = c("row", "col"), value.name = "bNTI")
w.bNTI.ITS.df <- melt(as.matrix(w.bNTI.ITS.d), varnames = c("row", "col"), value.name = "bNTI")
# RCbray
b_RC.df <- melt(as.matrix(b_rc_res), varnames = c("row", "col"), value.name = "RCbray")
f_RC.df <- melt(as.matrix(f_rc_res), varnames = c("row", "col"), value.name = "RCbray")


# prepare colomns of each pair
sample_names <- row.names(weighted.bNTI.16S) # get sample names
combinations <- expand.grid(sample_names, sample_names, stringsAsFactors = FALSE) # list all pairs
combinations <- combinations[combinations[,1] < combinations[,2],] # extract pairs (use only unique pairs)
colnames(combinations) <- c("row", "col")

# merge and use only unique pairs
# 16S
w.bNTI.16S.df <- merge(combinations, w.bNTI.16S.df, by = c("row", "col"), all.x = T, sort = T) # first, bNTI
w.bNTI.RCbray_16S <- merge(w.bNTI.16S.df, b_RC.df, by = c("row", "col"), all.x = T, sort = T) # next, RCbray

# ITS
w.bNTI.ITS.df <- merge(combinations, w.bNTI.ITS.df, by = c("row", "col"), all.x = T, sort = T) # first, bNTI
w.bNTI.RCbray_ITS <- merge(w.bNTI.ITS.df, b_RC.df, by = c("row", "col"), all.x = T, sort = T) # next, RCbray


### calculate the assembly process based on Stegen et al. and 
w.bNTI.RCbray_16S <- w.bNTI.RCbray_16S |>
  mutate(process = case_when(
    bNTI > 2 ~ "Heterogeneous Selection", # significantly more phylogenetic turnover than expected
    bNTI < -2 ~ "Homogeneous Selection", # significantly less phylogenetic turnover than expected
    abs(bNTI) < 2 & RCbray > 0.95 ~ "Dispersal Limitation", 
    abs(bNTI) < 2 & RCbray < -0.95 ~ "Homogenizing Dispersal",  
    abs(bNTI) < 2 & abs(RCbray) < 0.95 ~ "Undominated",
    TRUE ~ "Undefined" 
  ))


w.bNTI.RCbray_ITS <- w.bNTI.RCbray_ITS |>
  mutate(process = case_when(
    bNTI > 2 ~ "Heterogeneous Selection", # significantly more phylogenetic turnover than expected
    bNTI < -2 ~ "Homogeneous Selection", # significantly less phylogenetic turnover than expected
    abs(bNTI) < 2 & RCbray > 0.95 ~ "Dispersal Limitation", 
    abs(bNTI) < 2 & RCbray < -0.95 ~ "Homogenizing Dispersal",  
    abs(bNTI) < 2 & abs(RCbray) <= 0.95 ~ "Undominated",
    TRUE ~ "Undefined" 
  ))



### add "pair type" column
## annotate country, site, land use columns for the pairs
# prokaryotes
w.bNTI.RCbray_16S <- w.bNTI.RCbray_16S |>
  left_join(env %>% select(Sample, Country, Site, Landuse), by = c("row" = "Sample")) |>
  rename(row.country = Country,
         row.site = Site,
         row.landuse = Landuse
         ) |>
  left_join(env %>% select(Sample, Country, Site, Landuse), by = c("col" = "Sample")) |>
  rename(col.country = Country,
         col.site = Site,
         col.landuse = Landuse
  )

# fungi
w.bNTI.RCbray_ITS <- w.bNTI.RCbray_ITS |>
  left_join(env %>% select(Sample, Country, Site, Landuse), by = c("row" = "Sample")) |>
  rename(row.country = Country,
         row.site = Site,
         row.landuse = Landuse
  ) |>
  left_join(env %>% select(Sample, Country, Site, Landuse), by = c("col" = "Sample")) |>
  rename(col.country = Country,
         col.site = Site,
         col.landuse = Landuse
  )

## pair type (scale & land use)
# prokaryotes
w.bNTI.RCbray_16S <- w.bNTI.RCbray_16S |> 
  mutate(pair_scale = case_when(
    row.site == col.site ~ "within_site",
    row.site != col.site & row.country == col.country ~ "across_site_within_country",
    row.site != col.site & row.country != col.country ~ "across_site_btw_country",
    TRUE ~ NA
  ),
  pair_landu = case_when(
    row.landuse == col.landuse ~ row.landuse,
    row.landuse != col.landuse ~ "cross_landuse",
    TRUE ~ NA
  )
  )

# fungi
w.bNTI.RCbray_ITS <- w.bNTI.RCbray_ITS |> 
  mutate(pair_scale = case_when(
    row.site == col.site ~ "within_site",
    row.site != col.site & row.country == col.country ~ "across_site_within_country",
    row.site != col.site & row.country != col.country ~ "across_site_btw_country",
    TRUE ~ NA
  ),
  pair_landu = case_when(
    row.landuse == col.landuse ~ row.landuse,
    row.landuse != col.landuse ~ "cross_landuse",
    TRUE ~ NA
  )
  )




### save data
write.csv(w.bNTI.RCbray_16S, "13_AssemblyProcess_out/summary_bNTI_RCbray_prok.csv")
write.csv(w.bNTI.RCbray_ITS, "13_AssemblyProcess_out/summary_bNTI_RCbray_fungi.csv")

