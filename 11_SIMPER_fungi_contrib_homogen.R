####
#### R script for Ohigashi et al (2024)
#### Analysis on which ASVs do not contribute to heterogeneity within the group
#### 2024.09.20 written by Ohigashi
#### 2024.11.23 edited by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(tibble); packageVersion("tibble")
library(stringr); packageVersion("stringr")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")
source("Function/F1_HelperFunctions.R")
source("Function/F2_HelperFunctions_for_Visualization.R")


### load data
# fungal ASV
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)
f_ASV <- f_ASV.table[,1:(ncol(f_ASV.table)-7)]
f_taxonomy <- f_ASV.table[,(ncol(f_ASV.table)-6):ncol(f_ASV.table)]

# lifestyle list
lifestyle_list <- read.table("02_Function_analysis_out/FungalTraits_w_rarefied_ASV_table_fungi.txt", header = T, row.names = 1, sep = "\t")
lifestyle_list <- lifestyle_list |>
  rownames_to_column(var = "ASV")
lifestyle_list <- lifestyle_list |>
  dplyr::select(ASV, primary_lifestyle)

# soil metadata
DESIGN <- read.table("Data/soil_metadata.txt", header = T)
DESIGN <- DESIGN |>
  select(Sample, Country, Site, Landuse)


### format data
# fungal data
f_percent <- f_ASV/mean(colSums(f_ASV))*100 # convert to percentage
f_percent.t <- t(f_percent) # transpose


### analysis on contribution to heterogeneity within the group (i.e., Natural or Farm in a certain scale)
## across-site between the countries scale (done with the workstation in the Ushiolab)

# set category
asac_cat <- unique(DESIGN$Landuse)
res_dfs_asac <- list()

# looping for simper
for (i in asac_cat) {
  f_percent.t.df <- as.data.frame(f_percent.t)
  f_percent.t.df <- f_percent.t.df |> mutate(cat = DESIGN$Landuse)
  subdf <- f_percent.t.df |>
    filter(cat == i) |>
    select(-cat)
  subDESIGN <- DESIGN |> filter(Landuse == i)
  
  # simper
  sim_fungi <- with(subDESIGN, simper(subdf, Sample))
  sim_fungi_summary <- summary(sim_fungi)
  res_df <- data.frame(ASV = colnames(f_percent.t))
  for (j in 1:length(sim_fungi_summary)) {
    res_df0 <- data.frame(ASV = row.names(sim_fungi_summary[[j]]), X = sim_fungi_summary[[j]]$average)
    colnames(res_df0)[2] <- names(sim_fungi_summary)[j]
    res_df <- merge(res_df, res_df0, by = "ASV", sort = F)
  }
  
  res_dfs_asac[[i]] <- res_df
}


### categorize each pair of samples
combinations <- expand.grid(DESIGN$Sample, DESIGN$Sample, stringsAsFactors = FALSE) # list all pairs
combinations <- combinations[combinations[,1] < combinations[,2],] # extract pairs (removed the replaced pairs)
comb <- apply(combinations, 1, function(x) paste(x, collapse = "_"))
comb.df <- data.frame(comb = sort(comb))

# recognize land use for each pair
comb.df$landuse <- paste(DESIGN$Landuse[match(substr(comb.df$comb, 1, 3), DESIGN$Sample)],
                         DESIGN$Landuse[match(substr(comb.df$comb, 5, 7), DESIGN$Sample)],
                         sep = "_"
)
# recognize site for each pair
comb.df$site <- paste(DESIGN$Site[match(substr(comb.df$comb, 1, 3), DESIGN$Sample)],
                      DESIGN$Site[match(substr(comb.df$comb, 5, 7), DESIGN$Sample)],
                      sep = "_"
)
# recognize country for each pair
comb.df$country <- paste(DESIGN$Country[match(substr(comb.df$comb, 1, 3), DESIGN$Sample)],
                         DESIGN$Country[match(substr(comb.df$comb, 5, 7), DESIGN$Sample)],
                         sep = "_"
)

# set the categories
comb.df <- comb.df |>
  rowwise() |>
  mutate(scale_type = case_when(
    # different land use
    strsplit(landuse, "_")[[1]][1] != strsplit(landuse, "_")[[1]][2] ~ "diff_landuse",
    
    # same land use & same site
    strsplit(landuse, "_")[[1]][1] == strsplit(landuse, "_")[[1]][2] & strsplit(site, "_")[[1]][1] == strsplit(site, "_")[[1]][2] ~ 
      paste0("within_", strsplit(site, "_")[[1]][1], "_", strsplit(landuse, "_")[[1]][1]),
    
    # same land use & different site & same country
    strsplit(landuse, "_")[[1]][1] == strsplit(landuse, "_")[[1]][2] & strsplit(site, "_")[[1]][1] != strsplit(site, "_")[[1]][2] &
      strsplit(country, "_")[[1]][1] == strsplit(country, "_")[[1]][2] ~ 
      paste0("as_within_", strsplit(country, "_")[[1]][1], "_", strsplit(landuse, "_")[[1]][1]),
    
    # same land use & different site & between different countries
    strsplit(landuse, "_")[[1]][1] == strsplit(landuse, "_")[[1]][2] & strsplit(site, "_")[[1]][1] != strsplit(site, "_")[[1]][2] &
      strsplit(country, "_")[[1]][1] != strsplit(country, "_")[[1]][2] ~  
      paste0("as_ac_", str_split(landuse, "_")[[1]][1])
  ))

# remove "different landuse" rows
comb.df_samelu <- comb.df |> filter(scale_type != "diff_landuse")


### average the contribution of each ASV to the dissimilarity of the samples within the group
# using the category above
scale_types <- unique(comb.df_samelu$scale_type)
ave_dfs <- list()

# looping to get average of contribution within each scale type
for (k in scale_types){
  # calculate mean contribution
  comb_vec <- comb.df_samelu |>
    filter(scale_type == k) |>
    select(comb)
  if (str_detect(k, "Natural")) {
    subdf <- res_dfs_asac[["Natural"]] |> select(ASV, comb_vec$comb)
  } else {
    subdf <- res_dfs_asac[["Farm"]] |> select(ASV, comb_vec$comb)
  }
  meandf <- subdf |>
    mutate(mean_contrib = rowMeans(subdf[,-1], na.rm = TRUE)) |>
    select(ASV, mean_contrib)
  
  # calculate mean relative abundance of each ASV within each scale type
  sample_list <- unique(unlist(strsplit(comb_vec$comb, "_")))
  f_percent_sub <- f_percent |> select(sample_list)
  ASV_mean_rel <- f_percent_sub |>
    mutate(ASV = row.names(f_percent_sub),
           mean_rel = rowMeans(f_percent_sub, na.rm = TRUE)) |>
    select(ASV, mean_rel)
  
  # merge mean contribution, mean relative abundance, and primary lifestyle columns
  ave_df0 <- merge2(list(meandf, ASV_mean_rel, lifestyle_list), by = "ASV", sort = F)
  ave_df0 <- ave_df0 |> 
    mutate(primary_lifestyle = replace_na(primary_lifestyle, "unassigned"))
  
  # put it into the data frame list
  ave_dfs[[k]] <- ave_df0
}


# average ASV's contribution to dissimilarity in each scale category
wsD_Nat <- ave_dfs[["within_D_Natural"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsD_Nat_mean = mean(mean_contrib),
            wsD_Nat_med = median(mean_contrib))
wsD_Farm <- ave_dfs[["within_D_Farm"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsD_Farm_mean = mean(mean_contrib),
            wsD_Farm_med = median(mean_contrib))
wsE_Nat <- ave_dfs[["within_E_Natural"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsE_Nat_mean = mean(mean_contrib),
            wsE_Nat_med = median(mean_contrib))
wsE_Farm <- ave_dfs[["within_E_Farm"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsE_Farm_mean = mean(mean_contrib),
            wsE_Farm_med = median(mean_contrib))
wsF_Nat <- ave_dfs[["within_F_Natural"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsF_Nat_mean = mean(mean_contrib),
            wsF_Nat_med = median(mean_contrib))
wsF_Farm <- ave_dfs[["within_F_Farm"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsF_Farm_mean = mean(mean_contrib),
            wsF_Farm_med = median(mean_contrib))
wsG_Nat <- ave_dfs[["within_G_Natural"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsG_Nat_mean = mean(mean_contrib),
            wsG_Nat_med = median(mean_contrib))
wsG_Farm <- ave_dfs[["within_G_Farm"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsG_Farm_mean = mean(mean_contrib),
            wsG_Farm_med = median(mean_contrib))
wsH_Nat <- ave_dfs[["within_H_Natural"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsH_Nat_mean = mean(mean_contrib),
            wsH_Nat_med = median(mean_contrib))
wsH_Farm <- ave_dfs[["within_H_Farm"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(wsH_Farm_mean = mean(mean_contrib),
            wsH_Farm_med = median(mean_contrib))
as_wKenya_Natural <- ave_dfs[["as_within_Kenya_Natural"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(as_wKenya_Nat_mean = mean(mean_contrib),
            as_wKenya_Nat_med = median(mean_contrib))
as_wKenya_Farm <- ave_dfs[["as_within_Kenya_Farm"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(as_wKenya_Farm_mean = mean(mean_contrib),
            as_wKenya_Farm_med = median(mean_contrib))
as_wMalawi_Natural <- ave_dfs[["as_within_Malawi_Natural"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(as_wMalawi_Nat_mean = mean(mean_contrib),
            as_wMalawi_Nat_med = median(mean_contrib))
as_wMalawi_Farm <- ave_dfs[["as_within_Malawi_Farm"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(as_wMalawi_Farm_mean = mean(mean_contrib),
            as_wMalawi_Farm_med = median(mean_contrib))
asac_Natural <- ave_dfs[["as_ac_Natural"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(as_ac_Nat_mean = mean(mean_contrib),
            as_ac_Nat_med = median(mean_contrib))
asac_Farm <- ave_dfs[["as_ac_Farm"]] |>
  filter(mean_rel != 0) |>
  group_by(primary_lifestyle) |>
  summarize(as_ac_Farm_mean = mean(mean_contrib),
            as_ac_Farm_med = median(mean_contrib))

# make a data frame for observed averages of ASV contributions within each category
lifeave <- merge2(list(wsD_Nat, wsD_Farm, wsE_Nat, wsE_Farm,
                       wsF_Nat, wsF_Farm, wsG_Nat, wsG_Farm, wsH_Nat, wsH_Farm,
                       as_wKenya_Natural, as_wKenya_Farm, as_wMalawi_Natural, as_wMalawi_Farm,
                       asac_Natural, asac_Farm), by = "primary_lifestyle", sort = F, all = F)




### permutation test comparing the contributions to dissimilarity summed within fungal lifestyles between Natural & Farm
### heavy calculation!!! conducted with the workstation in the Ushiolab
## 1000-time permutation to get average contributions of each ASV when randomly shuffled
# number of cores to be used
# num_cores <- 96 # maybe too much
num_cores <- 32

# create cluster
cl <- makeCluster(num_cores)
# export packages, functions, and data to the cluster
clusterEvalQ(cl, {
  library(dplyr)
  library(vegan)
  library(tibble)
  library(stringr)
  library(tidyr)
})
clusterExport(cl, c("merge2", "DESIGN", "f_percent.t", "f_percent", "comb.df", "lifestyle_list"))

# pseudf <- pblapply(1:1000, function(x) {
# pseudf2 <- pblapply(1:50, function(x) { # trial
pseudf2 <- pblapply(1:1000, function(x) {
  set.seed(123 + x)
  # make pseudo "land use" under the fixed "site" as original
  DESIGN <- DESIGN |>
    group_by(Site) |>
    mutate(pse_landuse = sample(c(rep("pseNatural", 9), rep("pseFarm", 9)), replace = FALSE)) |>
    ungroup()
  
  # set category
  ps_cat <- unique(DESIGN$pse_landuse)
  res_dfs_ps <- list()
  
  # looping for simper
  for (j in ps_cat) {
    f_percent.t.df <- as.data.frame(f_percent.t)
    f_percent.t.df <- f_percent.t.df |> mutate(cat = DESIGN$pse_landuse)
    subdf <- f_percent.t.df |>
      filter(cat == j) |>
      select(-cat)
    subDESIGN <- DESIGN |> filter(pse_landuse == j)
    
    # simper
    sim_fungi <- with(subDESIGN, simper(subdf, Sample))
    sim_fungi_summary <- summary(sim_fungi)
    res_df <- data.frame(ASV = colnames(f_percent.t))
    for (k in 1:length(sim_fungi_summary)) {
      res_df0 <- data.frame(ASV = row.names(sim_fungi_summary[[k]]), X = sim_fungi_summary[[k]]$average)
      colnames(res_df0)[2] <- names(sim_fungi_summary)[k]
      res_df <- merge(res_df, res_df0, by = "ASV", sort = F)
    }
    
    res_dfs_ps[[j]] <- res_df
  }
  
  # categorize the scale type (pseudo ver.)
  # recognize land use for each pair
  comb.df_ps <- comb.df |>
    mutate(pse_landuse = paste(DESIGN$pse_landuse[match(substr(comb, 1, 3), DESIGN$Sample)],
                               DESIGN$pse_landuse[match(substr(comb, 5, 7), DESIGN$Sample)],
                               sep = "_"
    ))
  # set the scale types (pseudo)
  comb.df_ps <- comb.df_ps |>
    rowwise() |>
    mutate(pse_scale_type = case_when(
      # different land use
      strsplit(pse_landuse, "_")[[1]][1] != strsplit(pse_landuse, "_")[[1]][2] ~ "diff_pse_landuse",
      
      # same land use & same site
      strsplit(pse_landuse, "_")[[1]][1] == strsplit(pse_landuse, "_")[[1]][2] & strsplit(site, "_")[[1]][1] == strsplit(site, "_")[[1]][2] ~ 
        paste0("within_", strsplit(site, "_")[[1]][1], "_", strsplit(pse_landuse, "_")[[1]][1]),
      
      # same land use & different site & same country
      strsplit(pse_landuse, "_")[[1]][1] == strsplit(pse_landuse, "_")[[1]][2] & strsplit(site, "_")[[1]][1] != strsplit(site, "_")[[1]][2] &
        strsplit(country, "_")[[1]][1] == strsplit(country, "_")[[1]][2] ~ 
        paste0("as_within_", strsplit(country, "_")[[1]][1], "_", strsplit(pse_landuse, "_")[[1]][1]),
      
      # same land use & different site & between different countries
      strsplit(pse_landuse, "_")[[1]][1] == strsplit(pse_landuse, "_")[[1]][2] & strsplit(site, "_")[[1]][1] != strsplit(site, "_")[[1]][2] &
        strsplit(country, "_")[[1]][1] != strsplit(country, "_")[[1]][2] ~  
        paste0("as_ac_", str_split(pse_landuse, "_")[[1]][1])
    ))
  
  # remove "diff_pse_landuse" rows
  comb.df_ps_samelu <- comb.df_ps |> filter(pse_scale_type != "diff_pse_landuse")
  
  # looping to get average of contribution within each scale type (pseudo)
  ps_scale_types <- unique(comb.df_ps_samelu$pse_scale_type)
  ps_ave_dfs <- list()
  for (l in ps_scale_types){
    # calculate mean contribution
    comb_vec <- comb.df_ps_samelu |>
      filter(pse_scale_type == l) |>
      select(comb)
    if (str_detect(l, "Natural")) {
      subdf <- res_dfs_ps[["pseNatural"]] |> select(ASV, comb_vec$comb)
    } else {
      subdf <- res_dfs_ps[["pseFarm"]] |> select(ASV, comb_vec$comb)
    }
    meandf <- subdf |>
      mutate(mean_contrib = rowMeans(subdf[,-1], na.rm = TRUE)) |>
      select(ASV, mean_contrib)
    
    # calculate mean relative abundance of each ASV within each scale type
    sample_list <- unique(unlist(strsplit(comb_vec$comb, "_")))
    f_percent_sub <- f_percent |> select(sample_list)
    ASV_mean_rel <- f_percent_sub |>
      mutate(ASV = row.names(f_percent_sub),
             mean_rel = rowMeans(f_percent_sub, na.rm = TRUE)) |>
      select(ASV, mean_rel)
    
    # merge mean contribution and primary lifestyle columns
    ave_df0 <- merge2(list(meandf, ASV_mean_rel, lifestyle_list), by = "ASV", sort = F)
    # ave_df0 <- merge(meandf, lifestyle_list, by = "ASV", sort = F) # deleted to also consider abundance
    ave_df0 <- ave_df0 |> 
      mutate(primary_lifestyle = replace_na(primary_lifestyle, "unassigned"))
    
    # put it into the data frame list
    ps_ave_dfs[[l]] <- ave_df0
  }
  
  # calculate mean values of contributions grouped by primary lifestyles
  wsD_ps_Nat <- ps_ave_dfs[["within_D_pseNatural"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsD_ps_Nat_mean = mean(mean_contrib),
              wsD_ps_Nat_med = median(mean_contrib))
  wsD_ps_Farm <- ps_ave_dfs[["within_D_pseFarm"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsD_ps_Farm_mean = mean(mean_contrib),
              wsD_ps_Farm_med = median(mean_contrib))
  wsE_ps_Nat <- ps_ave_dfs[["within_E_pseNatural"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsE_ps_Nat_mean = mean(mean_contrib),
              wsE_ps_Nat_med = median(mean_contrib))
  wsE_ps_Farm <- ps_ave_dfs[["within_E_pseFarm"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsE_ps_Farm_mean = mean(mean_contrib),
              wsE_ps_Farm_med = median(mean_contrib))
  wsF_ps_Nat <- ps_ave_dfs[["within_F_pseNatural"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsF_ps_Nat_mean = mean(mean_contrib),
              wsF_ps_Nat_med = median(mean_contrib))
  wsF_ps_Farm <- ps_ave_dfs[["within_F_pseFarm"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsF_ps_Farm_mean = mean(mean_contrib),
              wsF_ps_Farm_med = median(mean_contrib))
  wsG_ps_Nat <- ps_ave_dfs[["within_G_pseNatural"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsG_ps_Nat_mean = mean(mean_contrib),
              wsG_ps_Nat_med = median(mean_contrib))
  wsG_ps_Farm <- ps_ave_dfs[["within_G_pseFarm"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsG_ps_Farm_mean = mean(mean_contrib),
              wsG_ps_Farm_med = median(mean_contrib))
  wsH_ps_Nat <- ps_ave_dfs[["within_H_pseNatural"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsH_ps_Nat_mean = mean(mean_contrib),
              wsH_ps_Nat_med = median(mean_contrib))
  wsH_ps_Farm <- ps_ave_dfs[["within_H_pseFarm"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(wsH_ps_Farm_mean = mean(mean_contrib),
              wsH_ps_Farm_med = median(mean_contrib))
  as_wKenya_ps_Natural <- ps_ave_dfs[["as_within_Kenya_pseNatural"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(as_wKenya_ps_Nat_mean = mean(mean_contrib),
              as_wKenya_ps_Nat_med = median(mean_contrib))
  as_wKenya_ps_Farm <- ps_ave_dfs[["as_within_Kenya_pseFarm"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(as_wKenya_ps_Farm_mean = mean(mean_contrib),
              as_wKenya_ps_Farm_med = median(mean_contrib))
  as_wMalawi_ps_Natural <- ps_ave_dfs[["as_within_Malawi_pseNatural"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(as_wMalawi_ps_Nat_mean = mean(mean_contrib),
              as_wMalawi_ps_Nat_med = median(mean_contrib))
  as_wMalawi_ps_Farm <- ps_ave_dfs[["as_within_Malawi_pseFarm"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(as_wMalawi_ps_Farm_mean = mean(mean_contrib),
              as_wMalawi_ps_Farm_med = median(mean_contrib))
  asac_ps_Natural <- ps_ave_dfs[["as_ac_pseNatural"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(as_ac_ps_Nat_mean = mean(mean_contrib),
              as_ac_ps_Nat_med = median(mean_contrib))
  asac_ps_Farm <- ps_ave_dfs[["as_ac_pseFarm"]] |>
    filter(mean_rel != 0) |>
    group_by(primary_lifestyle) |>
    summarize(as_ac_ps_Farm_mean = mean(mean_contrib),
              as_ac_ps_Farm_med = median(mean_contrib))
  
  
  pseudf0 <- merge2(list(wsD_ps_Nat, wsD_ps_Farm, wsE_ps_Nat, wsE_ps_Farm,
                         wsF_ps_Nat, wsF_ps_Farm, wsG_ps_Nat, wsG_ps_Farm, wsH_ps_Nat, wsH_ps_Farm,
                         as_wKenya_ps_Natural, as_wKenya_ps_Farm, as_wMalawi_ps_Natural, as_wMalawi_ps_Farm,
                         asac_ps_Natural, asac_ps_Farm), by = "primary_lifestyle", sort = F)
  for (i in 1:20) {
    gc()
  }
  
  return(pseudf0)
},

cl = cl
)
# stop cluster
stopCluster(cl)


### calculate the observed difference between Nat vs Farm
# make a vector for scales
scale_cat <- c("wsD", "wsE", "wsF", "wsG", "wsH", "wKenya", "wMalawi", "as_ac")

# make a data frame to store the result
obs.diff <- data.frame(matrix(ncol = length(scale_cat), nrow = length(lifeave$primary_lifestyle)))
rownames(obs.diff) <- lifeave$primary_lifestyle
colnames(obs.diff) <- scale_cat
  
# loop through each scale category
for (cat in scale_cat) {
  # extract the Natural and Farm columns for the current scale category
  col_nat <- lifeave |> select(contains(paste0(cat, "_Nat_mean"))) # mean
  # col_nat <- lifeave |> select(contains(paste0(cat, "_Nat_med"))) # median
  col_farm <- lifeave |> select(contains(paste0(cat, "_Farm_mean"))) # mean
  # col_farm <- lifeave |> select(contains(paste0(cat, "_Farm_med"))) # median
  # calculate the absolute difference and store in the list
  obs.diff[[cat]] <- abs(col_nat[[1]] - col_farm[[1]])
}

obs.diff <- obs.diff |> rownames_to_column("primary_lifestyle")


### calculate the pseudo difference between pseudo Nat vs pseudo Farm looping for each data frame
ps.dfs <- list()

for (i in 1:1000) {
  # extract one iteration from the pseudo result
  ori_ps <- pseudf2[[i]]
  # ori_ps <- pseudf[[i]]
  
  # make a data frame to store the result
  ps.diff <- data.frame(matrix(ncol = length(scale_cat), nrow = length(ori_ps$primary_lifestyle)))
  rownames(ps.diff) <- ori_ps$primary_lifestyle
  colnames(ps.diff) <- scale_cat
  
  # loop through each scale category
  for (cat in scale_cat) {
    # extract the Natural and Farm columns for the current scale category
    col_nat <- ori_ps %>% select(contains(paste0(cat, "_ps_Nat_mean"))) # mean
    # col_nat <- ori_ps %>% select(contains(paste0(cat, "_ps_Nat_med"))) # median
    col_farm <- ori_ps %>% select(contains(paste0(cat, "_ps_Farm_mean"))) # mean
    # col_farm <- ori_ps %>% select(contains(paste0(cat, "_ps_Farm_med"))) # median
    # calculate the absolute difference and store in the list
    ps.diff[[cat]] <- abs(col_nat[[1]] - col_farm[[1]])
  }
  
  ps.diff <- ps.diff |> rownames_to_column("primary_lifestyle")
  
  ps.dfs[[i]] <- ps.diff
}

# find common primary lifestyles that appear in all of 1000 permutations
common_lifestyles <- Reduce(intersect, lapply(ps.dfs, function(df) unique(df$primary_lifestyle)))

# filter data frames in ps.dfs and make a new list
ps.dfs.f <- lapply(ps.dfs, function(df) {
  df |>
    filter(primary_lifestyle %in% common_lifestyles)
})


### calculate p-values
# remove lifestyles not common in all 1000 permutations
obs.diff <- obs.diff |> filter(primary_lifestyle %in% common_lifestyles)
# data frame that has the same structure as obs.diff
perm_p <- obs.diff
perm_p[,-1] <- 0 # initialize the value

# looping to count how many times psuedo diff. was over observed diff.
for (ps_df in ps.dfs.f) {
  perm_p[,-1] <- perm_p[,-1] + (ps_df[,-1] >= obs.diff[,-1])
}

# calculate p-values by dividing the count by the number of permutation
perm_p_df2 <- cbind(perm_p[,1, drop = FALSE], perm_p[,-1] / length(ps.dfs.f))

# calculate adjusted p-values (because it was multiple testing) edited on 2024.11.23
perm_p_df2 <- read.csv("11_SIMPER_out/lifestyle_contrib_pvals_comparing_landuse_mean.csv")
perm_p_df3 <- perm_p_df2 |>
  mutate(across(-primary_lifestyle, 
              .fns = ~ p.adjust(., method = "BH"),
              .names = "{.col}.adj"))


### save data
saveRDS(res_dfs_acwc, file = "11_SIMPER_out/ASV_contrib.to.hetero_acrosssite_acrosscountry.obj") # downloaded from the workstation
saveRDS(ave_dfs, file = "11_SIMPER_out/ASV_ave.contrib_eachscale.obj")
saveRDS(lifeave, file = "11_SIMPER_out/lifestyle_mean_ave.contrib_observed.obj")
saveRDS(pseudf2, file = "11_SIMPER_out/lifestyle_mean_ave.contrib_pseudo2.obj") # downloaded from the workstation
# saveRDS(pseudf, file = "11_SIMPER_out/lifestyle_mean_ave.contrib_pseudo.obj") # downloaded from the workstation
# write.csv(perm_p_df, file = "11_SIMPER_out/lifestyle_contrib_pvals_comparing_landuse_median.csv", quote = F, row.names = F)
write.csv(perm_p_df2, file = "11_SIMPER_out/lifestyle_contrib_pvals_comparing_landuse_mean.csv", quote = F, row.names = F)
write.csv(perm_p_df3, file = "11_SIMPER_out/lifestyle_contrib_pvals.adj_comparing_landuse_mean.csv", quote = F, row.names = F)


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/11_2_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))




