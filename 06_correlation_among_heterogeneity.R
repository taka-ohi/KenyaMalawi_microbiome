####
#### R script for Ohigashi et al (2024)
#### Calculate correlation between distances to centroids of ...
#### 2024.07.08 written by Ohigashi
#### R 4.3.3
#### 


### load packages and functions
source("Function/F1_HelperFunctions.R")
library(dplyr); packageVersion("dplyr")


### load data
dists_within <- read.table("05_distance_to_centroid_out/distance_to_centroids_withinsite.txt", header = T, sep = "\t")
dists_across <- read.table("05_distance_to_centroid_out/distance_to_centroids_acrosssite.txt", header = T, sep = "\t")
env <- read.table("Data/soil_metadata.txt", header = T)

##### NOTE #####
# due to too many pairs of variables, I label the pairs as follows to reduce the efforts in labeling "meaningfully" in the scripts
# [NAMES]
# wcorr01: prokaryotic taxa x prokaryotic C cycle within site
# wcorr02: prokaryotic taxa x prokaryotic N cycle within site
# wcorr03: prokaryotic taxa x prokaryotic All functions within site
# wcorr04: fungal taxa x fungal lifestyles within site
# wcorr05: environment x prokaryotic taxa within site
# wcorr06: environment x fungal taxa within site
# wcoor07: environment x prokaryotic C cycle within site
# wcorr08: environment x prokaryotic N cycle within site
# wcorr09: environment x prokaryotic All functions within site
# wcorr10: environment x fungal lifestyles within site

# acorr01: prokaryotic taxa x prokaryotic C cycle across site
# acorr02: prokaryotic taxa x prokaryotic N cycle across site
# acorr03: prokaryotic taxa x prokaryotic All functions across site
# acorr04: fungal taxa x fungal lifestyles across site
# acorr05: environment x prokaryotic taxa across site
# acorr06: environment x fungal taxa across site
# acoor07: environment x prokaryotic C cycle across site
# acorr08: environment x prokaryotic N cycle across site
# acorr09: environment x prokaryotic All functions across site
# acorr10: environment x fungal lifestyles across site


### 1. correlation analysis for "within site" distance to centroids
# add column for environmental category
dists_within$Site <- env$Site
dists_within$Landuse <- env$Landuse

# prokaryotic taxa x prokaryotic C cycle
wcorr01 <- CorrInCat(dists_within, var1 = "proktaxa_within", var2 = "prok_cfunc_within",
                     method = "pearson", category = "Site")
wcorr01_is <- cor.test(dists_within$proktaxa_within, dists_within$prok_cfunc_within) # ignore site
wcorr01_is_summary <- cbind(wcorr01_is$estimate, wcorr01_is$p.value)
colnames(wcorr01_is_summary) <- c("r", "p.value")
row.names(wcorr01_is_summary) <- "All"
wcorr01 <- rbind(wcorr01, wcorr01_is_summary) # combine

# prokaryotic taxa x prokaryotic N cycle
wcorr02 <- CorrInCat(dists_within, var1 = "proktaxa_within", var2 = "prok_nfunc_within",
                     method = "pearson", category = "Site")
wcorr02_is <- cor.test(dists_within$proktaxa_within, dists_within$prok_nfunc_within) # ignore site
wcorr02_is_summary <- cbind(wcorr02_is$estimate, wcorr02_is$p.value)
colnames(wcorr02_is_summary) <- c("r", "p.value")
row.names(wcorr02_is_summary) <- "All"
wcorr02 <- rbind(wcorr02, wcorr02_is_summary) # combine

# prokaryotic taxa x prokaryotic All functions
wcorr03 <- CorrInCat(dists_within, var1 = "proktaxa_within", var2 = "prokallfunc_within",
                     method = "pearson", category = "Site")
wcorr03_is <- cor.test(dists_within$proktaxa_within, dists_within$prokallfunc_within) # ignore site
wcorr03_is_summary <- cbind(wcorr03_is$estimate, wcorr03_is$p.value)
colnames(wcorr03_is_summary) <- c("r", "p.value")
row.names(wcorr03_is_summary) <- "All"
wcorr03 <- rbind(wcorr03, wcorr03_is_summary) # combine

# fungal taxa x fungal lifestyles
wcorr04 <- CorrInCat(dists_within, var1 = "fungitaxa_within", var2 = "fungilife_within",
                     method = "pearson", category = "Site")
wcorr04_is <- cor.test(dists_within$fungitaxa_within, dists_within$fungilife_within) # ignore site
wcorr04_is_summary <- cbind(wcorr04_is$estimate, wcorr04_is$p.value)
colnames(wcorr04_is_summary) <- c("r", "p.value")
row.names(wcorr04_is_summary) <- "All"
wcorr04 <- rbind(wcorr04, wcorr04_is_summary) # combine

# environment x prokaryotic taxa
wcorr05 <- CorrInCat(dists_within, var1 = "env_within", var2 = "proktaxa_within",
                     method = "pearson", category = "Site")
wcorr05_is <- cor.test(dists_within$env_within, dists_within$proktaxa_within) # ignore site
wcorr05_is_summary <- cbind(wcorr05_is$estimate, wcorr05_is$p.value)
colnames(wcorr05_is_summary) <- c("r", "p.value")
row.names(wcorr05_is_summary) <- "All"
wcorr05 <- rbind(wcorr05, wcorr05_is_summary) # combine

# environment x fungal taxa
wcorr06 <- CorrInCat(dists_within, var1 = "env_within", var2 = "fungitaxa_within",
                     method = "pearson", category = "Site")
wcorr06_is <- cor.test(dists_within$env_within, dists_within$fungitaxa_within) # ignore site
wcorr06_is_summary <- cbind(wcorr06_is$estimate, wcorr06_is$p.value)
colnames(wcorr06_is_summary) <- c("r", "p.value")
row.names(wcorr06_is_summary) <- "All"
wcorr06 <- rbind(wcorr06, wcorr06_is_summary) # combine

# environment x prokaryotic C cycle
wcorr07 <- CorrInCat(dists_within, var1 = "env_within", var2 = "prok_cfunc_within",
                     method = "pearson", category = "Site")
wcorr07_is <- cor.test(dists_within$env_within, dists_within$prok_cfunc_within) # ignore site
wcorr07_is_summary <- cbind(wcorr07_is$estimate, wcorr07_is$p.value)
colnames(wcorr07_is_summary) <- c("r", "p.value")
row.names(wcorr07_is_summary) <- "All"
wcorr07 <- rbind(wcorr07, wcorr07_is_summary) # combine

# environment x prokaryotic N cycle
wcorr08 <- CorrInCat(dists_within, var1 = "env_within", var2 = "prok_nfunc_within",
                     method = "pearson", category = "Site")
wcorr08_is <- cor.test(dists_within$env_within, dists_within$prok_nfunc_within) # ignore site
wcorr08_is_summary <- cbind(wcorr08_is$estimate, wcorr08_is$p.value)
colnames(wcorr08_is_summary) <- c("r", "p.value")
row.names(wcorr08_is_summary) <- "All"
wcorr08 <- rbind(wcorr08, wcorr08_is_summary) # combine

# environment x prokaryotic All functions
wcorr09 <- CorrInCat(dists_within, var1 = "env_within", var2 = "prokallfunc_within",
                     method = "pearson", category = "Site")
wcorr09_is <- cor.test(dists_within$env_within, dists_within$prokallfunc_within) # ignore site
wcorr09_is_summary <- cbind(wcorr09_is$estimate, wcorr09_is$p.value)
colnames(wcorr09_is_summary) <- c("r", "p.value")
row.names(wcorr09_is_summary) <- "All"
wcorr09 <- rbind(wcorr09, wcorr09_is_summary) # combine

# environment x fungal lifestyles
wcorr10 <- CorrInCat(dists_within, var1 = "env_within", var2 = "fungilife_within",
                     method = "pearson", category = "Site")
wcorr10_is <- cor.test(dists_within$env_within, dists_within$fungilife_within) # ignore site
wcorr10_is_summary <- cbind(wcorr10_is$estimate, wcorr10_is$p.value)
colnames(wcorr10_is_summary) <- c("r", "p.value")
row.names(wcorr10_is_summary) <- "All"
wcorr10 <- rbind(wcorr10, wcorr10_is_summary) # combine



### 2. correlation analysis for "across site" distance to centroids
# prokaryotic taxa x prokaryotic C cycle
acorr01_test <- cor.test(dists_across$proktaxa_across, dists_across$prok_cfunc_across)
acorr01 <- cbind(acorr01_test$estimate, acorr01_test$p.value)
colnames(acorr01) <- c("r", "p.value")

# prokaryotic taxa x prokaryotic N cycle
acorr02_test <- cor.test(dists_across$proktaxa_across, dists_across$prok_nfunc_across)
acorr02 <- cbind(acorr02_test$estimate, acorr02_test$p.value)
colnames(acorr02) <- c("r", "p.value")

# prokaryotic taxa x prokaryotic all functions
acorr03_test <- cor.test(dists_across$proktaxa_across, dists_across$prokallfunc_across)
acorr03 <- cbind(acorr03_test$estimate, acorr03_test$p.value)
colnames(acorr03) <- c("r", "p.value")

# fungal taxa x fungal lifestyles
acorr04_test <- cor.test(dists_across$fungitaxa_across, dists_across$fungilife_across)
acorr04 <- cbind(acorr04_test$estimate, acorr04_test$p.value)
colnames(acorr04) <- c("r", "p.value")

# environment x prokaryotic taxa
acorr05_test <- cor.test(dists_across$env_across, dists_across$proktaxa_across)
acorr05 <- cbind(acorr05_test$estimate, acorr05_test$p.value)
colnames(acorr05) <- c("r", "p.value")

# environment x fungal taxa
acorr06_test <- cor.test(dists_across$env_across, dists_across$fungitaxa_across)
acorr06 <- cbind(acorr06_test$estimate, acorr06_test$p.value)
colnames(acorr06) <- c("r", "p.value")

# environment x prokaryotic C cycle
acorr07_test <- cor.test(dists_across$env_across, dists_across$prok_cfunc_across)
acorr07 <- cbind(acorr07_test$estimate, acorr07_test$p.value)
colnames(acorr07) <- c("r", "p.value")

# environment x prokaryotic N cycle
acorr08_test <- cor.test(dists_across$env_across, dists_across$prok_nfunc_across)
acorr08 <- cbind(acorr08_test$estimate, acorr08_test$p.value)
colnames(acorr08) <- c("r", "p.value")

# environment x prokaryotic all functions
acorr09_test <- cor.test(dists_across$env_across, dists_across$prokallfunc_across)
acorr09 <- cbind(acorr09_test$estimate, acorr09_test$p.value)
colnames(acorr09) <- c("r", "p.value")

# environment x fungal lifestyles
acorr10_test <- cor.test(dists_across$env_across, dists_across$fungilife_across)
acorr10 <- cbind(acorr10_test$estimate, acorr10_test$p.value)
colnames(acorr10) <- c("r", "p.value")


### save the data
dir.create("06_correlation_among_heterogeneity_out")
# create a list of results
res_list <- list(wcorr01, wcorr02, wcorr03, wcorr04, wcorr05,
                 wcorr06, wcorr07, wcorr08, wcorr09, wcorr10,
                 acorr01, acorr02, acorr03, acorr04, acorr05,
                 acorr06, acorr07, acorr08, acorr09, acorr10
                 )
# create a vector of file names
name_list <- c("cor_taxaprok_funcc_within",
               "cor_taxaprok_funcn_within",
               "cor_taxaprok_funcall_within",
               "cor_taxafungi_fungilife_within",
               "cor_env_taxaprok_within",
               "cor_env_taxafungi_within",
               "cor_env_funcc_within",
               "cor_env_funcn_within",
               "cor_env_funcall_within",
               "cor_env_fungilife_within",
               "cor_taxaprok_funcc_across",
               "cor_taxaprok_funcn_across",
               "cor_taxaprok_funcall_across",
               "cor_taxafungi_fungilife_across",
               "cor_env_taxaprok_across",
               "cor_env_taxafungi_across",
               "cor_env_funcc_across",
               "cor_env_funcn_across",
               "cor_env_funcall_across",
               "cor_env_fungilife_across"
               )
# save
for (i in 1:length(res_list)) {
  write.csv(res_list[[i]], sprintf("06_correlation_among_heterogeneity_out/%s.csv", name_list[i]), row.names = T, quote = F)
}


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/06_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

