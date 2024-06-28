####
#### R script for Ohigashi et al (2024)
#### combining 16S and ITS rarefied ASV data
#### 2024.06.28 written by Ohigashi
#### R 4.3.3
####

### load packages
library(dplyr)


### load data
setwd("~/Desktop/analysis/R_kenmal/KenyaMalawi_microbiome/")
b_ASV <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
f_ASV <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)

combined_table <- rbind(b_ASV, f_ASV)

### save the table under the "Data" directory for use
write.table(combined_table, "Data/rarefied_ASV_table_combined.txt", quote = F, row.names = T, sep = "\t")
