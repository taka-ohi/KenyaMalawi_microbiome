####
#### R script for Ohigashi et al (2024)
#### preparation for a fasta file of prokaryotes
#### 2024.06.29 written by Ohigashi
#### R 4.3.3
####

### load library and functions
library(tibble)
library(dplyr)

### load data
# rarefied ASV table
rared_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
# original sequence and taxa table (before rarefaction)
unrared_seqs <- read.table("01_DADA2_out/taxonomy_16S.txt", header = T)


### obtain sequence data of the rarefied ASV table
# include sequence data which is in the rowname in the seqs table
unrared_seqs <- unrared_seqs |> rownames_to_column(var = "Sequence")

# name ASVs
asvname <- paste0("Prok_ASV", seq(nrow(unrared_seqs))) # numbering
row.names(unrared_seqs) <- asvname

# include ASV names in the seqs table
unrared_seqs <- unrared_seqs |> rownames_to_column(var = "ASV")

# get a list of ASVs that remain in the rarefied ASV table
rared_asvname <- row.names(rared_ASV.table)

# filter the sequence table 
rared_seqs <- unrared_seqs |>
  dplyr::filter(ASV %in% rared_asvname)

# output fasta file
for (i in 1:nrow(rared_seqs)) {
  write(paste0(">", rared_seqs[i,1]), file = "01_DADA2_out/rarefied_ASV_seqs_16S.fasta", append = TRUE)
  write(rared_seqs[i,2], file= "01_DADA2_out/rarefied_ASV_seqs_16S.fasta", append = TRUE)
}

### save session info
setwd("~/Desktop/analysis/R_kenmal/KenyaMalawi_microbiome/")
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/02_1_SessionInfo_%s.txt", substr(Sys.time(), 1, 10))) 
