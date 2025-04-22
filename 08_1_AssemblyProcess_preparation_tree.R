####
#### R script for Ohigashi et al (2024)
#### Creating phylogenetic trees for 16S and ITS
#### 2024.10.08 written by Ohigashi
#### R 4.3.3
####


### load packages
library(seqinr); packageVersion("seqinr")
library(msa); packageVersion("msa")
library(phangorn); packageVersion("phangorn")


### load data
seq16s <- readDNAStringSet("02_Function_analysis_out/rarefied_ASV_seqs_16S.fasta") # 16S
seqITS <- readDNAStringSet("02_Function_analysis_out/rarefied_ASV_seqs_ITS.fasta") # ITS


### 1. 16S analysis
# multiple alignment
mult16s <- msa(seq16s, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)

# create phylogenetic tree (neighbor-joining)
phang.align16s <- as.phyDat(mult16s, type="DNA", names=getSequence(seq16s))
dm16s <- dist.ml(phang.align16s)
treeNJ16s <- NJ(dm16s) # Note, tip order != sequence order
fit16s <- pml(treeNJ16s, data=phang.align16s)

# negative edges length changed to 0!
fitGTR16s <- update(fit16s, k=4, inv=0.2)
fitGTR16s <- optim.pml(fitGTR16s, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))



### 2. ITS analysis
# multiple alignment
multITS <- msa(seqITS, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)

# create phylogenetic tree (neighbor-joining)
phang.alignITS <- as.phyDat(multITS, type="DNA", names=getSequence(seqITS))
dmITS <- dist.ml(phang.alignITS)
treeNJITS <- NJ(dmITS) # Note, tip order != sequence order
fitITS <- pml(treeNJITS, data=phang.alignITS)

# negative edges length changed to 0!
fitGTRITS <- update(fitITS, k=4, inv=0.2)
fitGTRITS <- optim.pml(fitGTRITS, model="GTR", optInv=TRUE, optGamma=TRUE,
                       rearrangement = "stochastic", control = pml.control(trace = 0))


### save data
# create a directory
dir.create("08_AssemblyProcess_out")

## 16S tree
saveRDS(mult16s, "08_AssemblyProcess_out/msa16s.obj") # alignment
saveRDS(treeNJ16s, "08_AssemblyProcess_out/treeNJ16s.obj") # NJ tree
saveRDS(fitGTR16s, "08_AssemblyProcess_out/treeGTR16s.obj") # GTR tree
write.tree(treeNJ16s, "08_AssemblyProcess_out/treeNJ16s.nwk") # NJ tree
write.tree(fitGTR16s$tree, "08_AssemblyProcess_out/treeGTR16s.nwk") # GTR tree

## ITS tree
saveRDS(multITS, "08_AssemblyProcess_out/msaITS.obj") # alignment
saveRDS(treeNJITS, "08_AssemblyProcess_out/treeNJITS.obj") # NJ tree
write.tree(treeNJITS, "08_AssemblyProcess_out/treeNJITS.nwk") # NJ tree
saveRDS(fitGTRITS, "08_AssemblyProcess_out/treeGTRITS.obj") # GTR tree
write.tree(fitGTRITS$tree, "08_AssemblyProcess_out/treeGTRITS.nwk") # GTR tree


detach("package:phangorn", unload=TRUE)



