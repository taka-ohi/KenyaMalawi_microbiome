####
#### R script for Ohigashi et al (2024)
#### Calculate correlation between distances to centroids of fungi and relative abundance of pathogenic fungi
#### 2024.07.10 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
source("Function/F1_HelperFunctions.R")

### load data
dists_within <- read.table("05_distance_to_centroid_out/distance_to_centroids_withinsite.txt", header = T, sep = "\t")
dists_across <- read.table("05_distance_to_centroid_out/distance_to_centroids_acrosssite.txt", header = T, sep = "\t")
env <- read.table("Data/soil_metadata.txt", header = T)
fungilife <- read.table("02_Function_analysis_out/FungalTraits_specific_functions.txt", header = T)


### calculation of correlation within site
# add column for environmental category
dists_within$Site <- env$Site
dists_within$Landuse <- env$Landuse
# add column for pathogen abundance
dists_within$Pathotroph <- fungilife$Pathotroph

# fungal taxa x fungal pathogen abundance
cor_sites <- CorrInCat(dists_within, var1 = "fungitaxa_within", var2 = "Pathotroph",
                       method = "pearson", category = "Site")
cor_is <- cor.test(dists_within$fungitaxa_within, dists_within$Pathotroph) # ignore site
cor_is_summary <- cbind(cor_is$estimate, cor_is$p.value)
colnames(cor_is_summary) <- c("r", "p.value")
row.names(cor_is_summary) <- "All"
cor_sites <- rbind(cor_sites, cor_is_summary) # combine


### calculation of correlation across site
# add column for pathogen abundance
dists_across$Pathotroph <- fungilife$Pathotroph

# fungal taxa x fungal pathogen abundance
cor_abun <- cor.test(dists_across$fungitaxa_across, dists_across$Pathotroph) # ignore site
cor_abun_summary <- cbind(cor_abun$estimate, cor_abun$p.value)
colnames(cor_abun_summary) <- c("r", "p.value")


### save the data
dir.create("07_correlation_fungtaxahetero_pathogenabund_out")
write.csv(cor_sites, "07_correlation_fungtaxahetero_pathogenabund_out/cor_taxafungi_pathoabun_within.csv", row.names = T, quote = F)
write.csv(cor_abun_summary, "07_correlation_fungtaxahetero_pathogenabund_out/cor_taxafungi_pathoabun_across.csv", row.names = T, quote = F)


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/07_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

