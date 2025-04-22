####
#### R script for Ohigashi et al (2024)
#### quality filtering, denoizing, and creating ASV table for ITS community using dada2
#### 2024.06.28 written by Ohigashi; 2025.04.15 edited by Ohigashi for rarefaction curves
#### R 4.3.3
#### README: be sure that ITS fastq files are in the "Data/2020kenya_malawi_ITS" directory


### load library and functions
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(tibble); packageVersion("tibble")

### load database (UNITE ver 8.3; for ITS)
# obtained from here https://doi.plutof.ut.ee/doi/10.15156/BIO/1280089
DATABASE <- "~/Desktop/analysis/Database/sh_general_release_s_10.05.2021/sh_general_release_dynamic_s_10.05.2021.fasta" # set your UNITE directory


### performing DADA2
# set working directory and list files
setwd("Data/2020kenya_malawi_ITS/")  ## CHANGE ME to the directory containing the fastq files.
filez <- list.files()

# prepare directory for output files
# dir.create("../../01_DADA2_out") # no need if you did in "01_1_DADA2_16S.R"

# record path
fnFs <- sort(list.files(getwd(), pattern = ".fastq", full.names = TRUE))

# primers sequence
FWD <- "TCCGTAGGTGAACCTGCGG" ## ITS1
REV <- "GCTGCGTTCTTCATCGATGC" ## ITS2

# create a function to know complement and reverse complement
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# get sample names
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)

# set path for samples to be filtered
fnFs.filtN <- file.path(getwd(), "filtN", paste0(sample.names, ".fastq.gz"))

# “pre-filter” the sequences just to remove those with ambiguous bases (Ns)
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

# Identifying and counting the primers on FASTQ files
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]) )

# set the directory of cutadapt command
# see https://cutadapt.readthedocs.io/en/stable/installation.html 
cutadapt <- "/opt/homebrew/Caskroom/miniconda/base/envs/cutadapt/bin/cutadapt" # CHANGE to the cutadapt path on your machine
system2(cutadapt, args = "--version")

# set paths for cutadapted files
path.cut <- file.path(getwd(), "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnFs)) # dummy file

# set reverse-oriented complement for the primers
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
# R2.flags <- paste("-G", REV, "-A", FWD.RC) # no need for single-end sequences

# Run cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2,
                             "-o", fnFs.cut[i],  # output file                             
                             fnFs.filtN[i])) # input files
}

# check if all the adapters are trimmed
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))

# list cutadapted files
cutFs <- sort(list.files(path.cut, pattern = ".fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs)) # set path
out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2,
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

# calculate the error rates and save an image
errF <- learnErrors(filtFs, multithread = TRUE)
p <- plotErrors(errF, nominalQ = TRUE)
ggsave(plot=p, file="../../01_DADA2_out/errors_ITS.png")

# Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

# Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE) 

# Construct Sequence Table
# see https://github.com/benjjneb/dada2/issues/384
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))

# Track reads through this workflow
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim)) # If processing a single sample, remove the sapply calls: e.g. replace
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,file="../../01_DADA2_out/track_ITS.txt", quote=F)


### Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, DATABASE, multithread = TRUE, tryRC = TRUE)
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) # check
write.table(taxa, file="../../01_DADA2_out/taxonomy_ITS.txt", quote=F)
write.table(seqtab.nochim, file="../../01_DADA2_out/seqtabnochim_ITS.txt", quote=F)

### Create ASV table
samples.out<-rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("Fung_ASV", seq(ntaxa(ps)))
ps

otu_table.t<-t(ps@otu_table)
ps.t<-cbind(otu_table.t,ps@tax_table)
write.table(ps.t,  file="../../01_DADA2_out/ASV_table_ITS.txt", quote=F, sep = "\t")

### create rarefied ASV table
# load the raw ASV table
ps.t <- read.table("../../01_DADA2_out/ASV_table_ITS.txt", header = T)

# remove "g__" or "p__" things from the taxa
# Define a function to remove the prefix
remove_prefix <- function(x) {
  sub("^._{2}", "", x)
}
# Apply the function to each column
ps.t_cleaned <- ps.t %>%
  mutate(across((ncol(ps.t)-6):(ncol(ps.t)), remove_prefix))

ps.t_wo_others <- ps.t_cleaned |>
  dplyr::filter(Kingdom == "Fungi") # removed the ones which were not assigned as "Fungi"
  
# get ASV count table
f_ASV <- ps.t_wo_others[,1:(ncol(ps.t_wo_others)-7)]
f_ASV.t <- t(f_ASV)

# rarefy by the lowest number of reads in a sample
rared_f_ASV.t <- rrarefy(f_ASV.t, min(rowSums(f_ASV.t)))

rared_f_ASV <- t(rared_f_ASV.t)
rared_f_taxa <- ps.t_wo_others[,(ncol(ps.t_wo_others)-6):(ncol(ps.t_wo_others))] # get taxa info

# remove rows that contains only 0 counts after rarefaction
rared_f_ASV_wo_0 <- rared_f_ASV[apply(rared_f_ASV, 1, function(row) !all(row == 0)), ]
rared_f_ASV_wo_0 <- as.data.frame(rared_f_ASV_wo_0)

# combine the count data and taxa data
rared_f_ASV_wo_0 <- rared_f_ASV_wo_0 |> rownames_to_column(var = "ASV")
rared_f_taxa <- rared_f_taxa |> rownames_to_column(var = "ASV")

rared_f_ASV.table <- merge(rared_f_ASV_wo_0, rared_f_taxa, by = "ASV", all.x = T, all.y = F, sort = F) # combine
rared_f_ASV.table <- rared_f_ASV.table |> column_to_rownames(var = "ASV")

write.table(rared_f_ASV.table,  file="../../01_DADA2_out/rarefied_ASV_table_ITS.txt", quote=F, sep = "\t")


### create rarefaction curves ###
# get raw ASV data
ps.t <- read.table("01_DADA2_out/ASV_table_ITS.txt", header = T)
# remove "g__" or "p__" things from the taxa
# Define a function to remove the prefix
remove_prefix <- function(x) {
  sub("^._{2}", "", x)
}
# Apply the function to each column
ps.t_cleaned <- ps.t %>%
  mutate(across((ncol(ps.t)-6):(ncol(ps.t)), remove_prefix))

ps.t_wo_others <- ps.t_cleaned |>
  dplyr::filter(Kingdom == "Fungi") # removed the ones which were not assigned as "Fungi"

# get ASV count table
f_ASV <- ps.t_wo_others[,1:(ncol(ps.t_wo_others)-7)]
f_ASV.t <- t(f_ASV)

# get rarafaction curve data
f_rarecurve <- rarecurve(f_ASV.t, step = 100, sample = min(rowSums(f_ASV.t)),
                         # xlab = "Number of reads", ylab = "Number of ASVs",
                         # col = "blue", cex = 0.6, label = FALSE,
                         tidy = T)
# plot rarefaction curve
p_f_rarecurve <- ggplot(data = f_rarecurve, aes(x = Sample, y = Species, group = Site)) +
  geom_line(color = "blue") +
  geom_vline(xintercept = min(rowSums(f_ASV.t))) +
  theme_classic() +
  labs(x = "Number of reads sampled", y = "Observed ASV richness", title = "Rarefaction curves for fungal communities")
# p_f_rarecurve

# save the curve
write.csv(f_rarecurve, file = "01_DADA2_out/rarefaction_curve_fungi.csv", quote = F, row.names = F)
ggsave("01_DADA2_out/rarefaction_curve_fungi.png", plot = p_f_rarecurve,
       width = 8, height = 7, bg = "white")
saveRDS(p_f_rarecurve, "01_DADA2_out/rarefaction_curve_fungi.rds")


### save session info
setwd("~/Desktop/analysis/R_kenmal/KenyaMalawi_microbiome/")
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/01_SessionInfo_%s.txt", substr(Sys.time(), 1, 10))) 

                    