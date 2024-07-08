####
#### R script for Ohigashi et al (2024)
#### Calculate distance to centroids for ...
#### 2024.07.04 written by Ohigashi
#### R 4.3.3
#### 


### load packages and functions
source("Function/F1_HelperFunctions.R")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")


### load data and prepare for use
# environtment table
env_data <- read.table("Data/soil_metadata.txt", header = T)
env_category <- env_data[,1:4] # extract category data
env_vars <- env_data[, 5:9] # extract soil env data
scale_values <- function(x){(x-min(x))/(max(x)-min(x))} # scale the variables so that they distribute from 0 to 1
s_env_vars <- scale_values(env_vars)

# prokaryotes ASV table
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
b_ASV <- b_ASV.table[,1:(ncol(b_ASV.table)-7)]
b_ASV.t <- t(b_ASV)

# fungi ASV table
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)
f_ASV <- f_ASV.table[,1:(ncol(f_ASV.table)-7)]
f_ASV.t <- t(f_ASV)

# prokaryotes all functions
prok_func_all <- read.table("02_Function_analysis_out/PICRUSt2_full_function.txt", header = T)

# prokaryotic C and N functions
prok_cn_cycle <- read.table("02_Function_analysis_out/PICRUSt2_CN_function_table.txt", header = T, sep = "\t")
# subset C cycle table
c_cycle <- prok_cn_cycle |>
  dplyr::filter(Cycle == "Carbon") |> 
  na.omit() |>
  dplyr::select(-KO, -Symbol, -Name, -Cycle, -Function, -Pathway)
c_cycle.t <- t(c_cycle)
# subset N cycle table
n_cycle <- prok_cn_cycle |>
  dplyr::filter(Cycle == "Nitrogen") |> 
  na.omit() |>
  dplyr::select(-KO, -Symbol, -Name, -Cycle, -Function, -Pathway)
n_cycle.t <- t(n_cycle)

# fungal lifestyles
fungilife <- read.table("02_Function_analysis_out/FungalTraits_primarilifestyle_percent.txt", header = T)
fungilife <- fungilife |> column_to_rownames(var = "primary_lifestyle")
fungilife.t <- t(fungilife)



### calculate distance to centroids of taxa, functions, and environments for "across-site" scale (just by land use)
# prokaryotes taxa
proktaxa_across <- DistToCent(b_ASV.t, method = "bray", group = env_category$Landuse,
                              name = "proktaxa_across")

# fungi taxa
fungitaxa_across <- DistToCent(f_ASV.t, method = "bray", group = env_category$Landuse,
                               name = "fungitaxa_across")

# prokaryotes all functions
prokallfunc_across <- DistToCent(prok_func_all, method = "bray", group = env_category$Landuse,
                                 name = "prokallfunc_across")

# prokaryotes C cycle
prok_cfunc_across <- DistToCent(c_cycle.t, method = "bray", group = env_category$Landuse,
                                name = "prok_cfunc_across")

# prokaryotes N cycle
prok_nfunc_across <- DistToCent(n_cycle.t, method = "bray", group = env_category$Landuse,
                                name = "prok_nfunc_across")

# fungi lifestyle
fungilife_across <- DistToCent(fungilife.t, method = "bray", group = env_category$Landuse,
                               name = "fungilife_across")

# environmental data
env_across <- DistToCent(s_env_vars, method = "bray", group = env_category$Landuse,
                         name = "env_across")

# combine the data
dists_across <- cbind(proktaxa_across,
                      fungitaxa_across,
                      prokallfunc_across,
                      prok_cfunc_across,
                      prok_nfunc_across,
                      fungilife_across,
                      env_across
                      )


### calculate distance to centroids of taxa, functions, and environments for "within-site" scale (by land use in each site)
# add "Treatment" column to the env_category table
env_category <- env_category |>
  mutate(Treatment = paste(Site, Landuse, sep = "_"))

# prokaryotes taxa
proktaxa_within <- DistToCent(b_ASV.t, method = "bray", group = env_category$Treatment,
                              name = "proktaxa_within")

# fungi taxa
fungitaxa_within <- DistToCent(f_ASV.t, method = "bray", group = env_category$Treatment,
                               name = "fungitaxa_within")

# prokaryotes all functions
prokallfunc_within <- DistToCent(prok_func_all, method = "bray", group = env_category$Treatment,
                                 name = "prokallfunc_within")

# prokaryotes C cycle
prok_cfunc_within <- DistToCent(c_cycle.t, method = "bray", group = env_category$Treatment,
                                name = "prok_cfunc_within")

# prokaryotes N cycle
prok_nfunc_within <- DistToCent(n_cycle.t, method = "bray", group = env_category$Treatment,
                                name = "prok_nfunc_within")

# fungi lifestyle
fungilife_within <- DistToCent(fungilife.t, method = "bray", group = env_category$Treatment,
                               name = "fungilife_within")

# environmental data
env_within <- DistToCent(s_env_vars, method = "bray", group = env_category$Treatment,
                         name = "env_within")

# combine the data
dists_within <- cbind(proktaxa_within,
                      fungitaxa_within,
                      prokallfunc_within,
                      prok_cfunc_within,
                      prok_nfunc_within,
                      fungilife_within,
                      env_within
                      )

### save the data
dir.create("05_distance_to_centroid_out")
write.table(dists_within, file = "05_distance_to_centroid_out/distance_to_centroids_withinsite.txt", row.names = T, quote = F, sep = "\t")
write.table(dists_across, file = "05_distance_to_centroid_out/distance_to_centroids_acrosssite.txt", row.names = T, quote = F, sep = "\t")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/05_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
