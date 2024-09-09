

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

# prokaryotes ASV table
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
b_ASV <- b_ASV.table[,1:(ncol(b_ASV.table)-7)]
b_ASV.t <- t(b_ASV)

# fungi ASV table
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)
f_ASV <- f_ASV.table[,1:(ncol(f_ASV.table)-7)]
f_ASV.t <- t(f_ASV)


# fungal lifestyles
fungilife <- read.table("00_1_check_chao/FungalTraits_primarilifestyle_aggregated.txt", header = T)
fungilife <- fungilife |> column_to_rownames(var = "primary_lifestyle")
fungilife.t <- t(fungilife)



### calculate distance to centroids of taxa, functions, and environments for "across-site" scale (just by land use)
# prokaryotes taxa
proktaxa_across <- DistToCent(b_ASV.t, method = "chao", group = env_category$Landuse,
                              name = "proktaxa_across")

# fungi taxa
fungitaxa_across <- DistToCent(f_ASV.t, method = "chao", group = env_category$Landuse,
                               name = "fungitaxa_across")


# fungi lifestyle
fungilife_across <- DistToCent(fungilife.t, method = "bray", group = env_category$Landuse,
                               name = "fungilife_across")

# combine the data
dists_across <- cbind(proktaxa_across,
                      fungitaxa_across,
                      fungilife_across
)


### calculate distance to centroids of taxa, functions, and environments for "within-site" scale (by land use in each site)
# add "Treatment" column to the env_category table
env_category <- env_category |>
  mutate(Treatment = paste(Site, Landuse, sep = "_"))

# prokaryotes taxa
proktaxa_within <- DistToCent(b_ASV.t, method = "chao", group = env_category$Treatment,
                              name = "proktaxa_within")

# fungi taxa
fungitaxa_within <- DistToCent(f_ASV.t, method = "chao", group = env_category$Treatment,
                               name = "fungitaxa_within")


# fungi lifestyle
fungilife_within <- DistToCent(fungilife.t, method = "bray", group = env_category$Treatment,
                               name = "fungilife_within")


# combine the data
dists_within <- cbind(proktaxa_within,
                      fungitaxa_within,
                      fungilife_within
)

### save the data
write.table(dists_within, file = "00_1_check_chao/distance_to_centroids_withinsite_re.txt", row.names = T, quote = F, sep = "\t")
write.table(dists_across, file = "00_1_check_chao/distance_to_centroids_acrosssite_re.txt", row.names = T, quote = F, sep = "\t")


