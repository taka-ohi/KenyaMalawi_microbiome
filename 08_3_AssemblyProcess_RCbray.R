####
#### R script for Ohigashi et al (2024)
#### Calculate RCbray
#### 2024.10.20 written by Ohigashi
#### R 4.3.3
#### ref https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r


### load package and functions
source("Function/F3_HelperFunctions_RC.R")
library(picante); packageVersion("picante")
library(pbapply); packageVersion("pbapply")


### load data
# count table
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)

# use only abundance data
b_ASV <- b_ASV.table[,-((ncol(b_ASV.table)-6):ncol(b_ASV.table))]
f_ASV <- f_ASV.table[,-((ncol(f_ASV.table)-6):ncol(f_ASV.table))]

# transpose
b_ASV.t <- t(b_ASV)
f_ASV.t <- t(f_ASV)


### calculate Raup-Crick bray
## 16S
b_rc_res <- raup_crick_abundance_parallel(spXsite = b_ASV.t, 
                                          plot_names_in_col1 = FALSE,
                                          reps = 9999,
                                          cl = 32)
saveRDS(b_rc_res, "08_AssemblyProcess_out/RCbray_prok.obj")
write.csv(b_rc_res, "08_AssemblyProcess_out/RCbray_prok.csv", quote = F)


## ITS
f_rc_res <- raup_crick_abundance_parallel(spXsite = f_ASV.t, 
                                          plot_names_in_col1 = FALSE,
                                          reps = 9999,
                                          cl = 32)
saveRDS(f_rc_res, "08_AssemblyProcess_out/RCbray_fungi.obj")
write.csv(f_rc_res, "08_AssemblyProcess_out/RCbray_fungi.csv", quote = F)


