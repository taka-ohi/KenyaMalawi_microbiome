####
#### R script for Ohigashi et al (2024)
#### calculate alpha diversity indeces for both prokaryotes and fungi
#### 2024.07.02 written by Ohigashi
#### R 4.3.3
####

### load packages
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")


### load rarefied ASV count data
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)


### format data
# prokaryotes
b_ASV <- b_ASV.table[, 1:(ncol(b_ASV.table)-7)] # extract the count part
b_ASV.t <- t(b_ASV) # transpose

# fungi
f_ASV <- f_ASV.table[, 1:(ncol(f_ASV.table)-7)] # extract the count part
f_ASV.t <- t(f_ASV) # transpose


### calculate diversity indices
# prokaryotes
b_div <- as.data.frame(cbind(
  shannon = diversity(b_ASV.t, index = "shannon", base = 2), 
  simpson = diversity(b_ASV.t, index = "simpson"), 
  invsimpson = diversity(b_ASV.t, index = "invsimpson"), 
  fisher = fisher.alpha(b_ASV.t)
))

# fungi
f_div <- as.data.frame(cbind(
  shannon = diversity(f_ASV.t, index = "shannon", base = 2), 
  simpson = diversity(f_ASV.t, index = "simpson"), 
  invsimpson = diversity(f_ASV.t, index = "invsimpson"), 
  fisher = fisher.alpha(f_ASV.t)
))


### save data
write.table(b_div, "01_DADA2_out/alpha_diversity_16S.txt", quote = F, row.names = T, sep = "\t")
write.table(f_div, "01_DADA2_out/alpha_diversity_ITS.txt", quote = F, row.names = T, sep = "\t")

# you can add these data to the environmental data ("Data/soil_metadata.txt")

