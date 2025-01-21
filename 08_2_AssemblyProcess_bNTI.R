
#### R script for Ohigashi et al (2024)
#### calculation of bNTI
#### 2024.10.15 written by Ohigashi
#### R 4.3.3
#### ref https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r 
#### NOTE: using workstation in the Ushio Lab


### load package and functions
library(ape); packageVersion("ape")
library(picante); packageVersion("picante")
library(pbapply); packageVersion("pbapply")

### load data
# phylogenetic tree
treeNJ16s <- read.tree("08_AssemblyProcess_out/treeNJ16s.nwk") # 16S
treeNJITS <- read.tree("08_AssemblyProcess_out/treeNJITS.nwk") # ITS

# count table
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)

b_ASV <- b_ASV.table[,-((ncol(b_ASV.table)-6):ncol(b_ASV.table))]
f_ASV <- f_ASV.table[,-((ncol(f_ASV.table)-6):ncol(f_ASV.table))]


### format data
match.phylo.asv16S <- match.phylo.data(treeNJ16s, b_ASV)
match.phylo.asvITS <- match.phylo.data(treeNJITS, f_ASV)


### 1. calculate empirical (observed) betaMNTD
## 1-1. 16S
beta.mntd.weighted16S <-  as.matrix(comdistnt(t(match.phylo.asv16S$data),cophenetic(match.phylo.asv16S$phy),abundance.weighted=T))
dim(beta.mntd.weighted16S)
write.csv(beta.mntd.weighted16S,"08_AssemblyProcess_out/betaMNTD_weighted16S.csv", quote=F)

identical(colnames(match.phylo.asv16S$data),colnames(beta.mntd.weighted16S)) # just a check, should be TRUE
identical(colnames(match.phylo.asv16S$data),rownames(beta.mntd.weighted16S)) # just a check, should be TRUE


## 1-2. ITS
beta.mntd.weightedITS <-  as.matrix(comdistnt(t(match.phylo.asvITS$data),cophenetic(match.phylo.asvITS$phy),abundance.weighted=T))
dim(beta.mntd.weightedITS)
write.csv(beta.mntd.weightedITS,"08_AssemblyProcess_out/betaMNTD_weightedITS.csv", quote=F)

identical(colnames(match.phylo.asvITS$data),colnames(beta.mntd.weightedITS)) # just a check, should be TRUE
identical(colnames(match.phylo.asvITS$data),rownames(beta.mntd.weightedITS)) # just a check, should be TRUE


### 2. calculate randomized betaMNTD
# set number of randomizations
beta.reps = 999

## 2-1. 16S
# parallel computation for randomized weighted beta MNTD
rand.w.bMNTD.comp16S.tmp <- pblapply(1:beta.reps, function(rep) {
  # set different ways to generate random values for each thread
  set.seed(123 + rep)
  
  # calculate 
  as.matrix(comdistnt(
    t(match.phylo.asv16S$data),
    taxaShuffle(cophenetic(match.phylo.asv16S$phy)),
    abundance.weighted = TRUE,
    exclude.conspecifics = FALSE
  ))
}, cl = 32)

# convert the list to 3-dimension array
rand.w.bMNTD.comp16S <- array(unlist(rand.w.bMNTD.comp16S.tmp),
                              dim = c(ncol(match.phylo.asv16S$data),
                                      ncol(match.phylo.asv16S$data),
                                      beta.reps))

# calculate
weighted.bNTI.16S <- matrix(c(NA),nrow=ncol(match.phylo.asv16S$data),ncol=ncol(match.phylo.asv16S$data))
dim(weighted.bNTI.16S)

for (columns in 1:(ncol(match.phylo.asv16S$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.asv16S$data)) {
    rand.vals <- rand.w.bMNTD.comp16S[rows,columns,]
    weighted.bNTI.16S[rows,columns] = (beta.mntd.weighted16S[rows,columns] - mean(rand.vals)) / sd(rand.vals)
    rm("rand.vals")
  }
}

rownames(weighted.bNTI.16S) <- colnames(match.phylo.asv16S$data)
colnames(weighted.bNTI.16S) <- colnames(match.phylo.asv16S$data)

write.csv(weighted.bNTI.16S,"08_AssemblyProcess_out/weighted_bNTI_16S.csv",quote=F)



## 2-2. ITS
# parallel computation for randomized weighted beta MNTD
rand.w.bMNTD.compITS.tmp <- pblapply(1:beta.reps, function(rep) {
  # set different ways to generate random values for each thread
  set.seed(123 + rep)
  
  # calculate 
  as.matrix(comdistnt(
    t(match.phylo.asvITS$data),
    taxaShuffle(cophenetic(match.phylo.asvITS$phy)),
    abundance.weighted = TRUE,
    exclude.conspecifics = FALSE
  ))
}, cl = 32)

# convert the list to 3-dimension array
rand.w.bMNTD.compITS <- array(unlist(rand.w.bMNTD.compITS.tmp),
                              dim = c(ncol(match.phylo.asvITS$data),
                                      ncol(match.phylo.asvITS$data),
                                      beta.reps))

# calculate
weighted.bNTI.ITS <- matrix(c(NA),nrow=ncol(match.phylo.asvITS$data),ncol=ncol(match.phylo.asvITS$data))
dim(weighted.bNTI.ITS)

for (columns in 1:(ncol(match.phylo.asvITS$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.asvITS$data)) {
    rand.vals <- rand.w.bMNTD.compITS[rows,columns,]
    weighted.bNTI.ITS[rows,columns] = (beta.mntd.weightedITS[rows,columns] - mean(rand.vals)) / sd(rand.vals)
    rm("rand.vals")
  }
}

rownames(weighted.bNTI.ITS) <- colnames(match.phylo.asvITS$data)
colnames(weighted.bNTI.ITS) <- colnames(match.phylo.asvITS$data)

write.csv(weighted.bNTI.ITS,"08_AssemblyProcess_out/weighted_bNTI_ITS.csv",quote=F)
