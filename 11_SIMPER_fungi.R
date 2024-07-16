####
#### R script for Ohigashi et al (2024)
#### SIMPER analysis for fungal ASVs accumulated by fungal lifestyles
#### 2024.07.16 written by Ohigashi
#### R 4.3.3
####


### load packages
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(tibble); packageVersion("tibble")
library(stringr); packageVersion("stringr")


### load data
# fungal ASV
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)
f_ASV <- f_ASV.table[,1:(ncol(f_ASV.table)-7)]
f_taxonomy <- f_ASV.table[,(ncol(f_ASV.table)-6):ncol(f_ASV.table)]

# lifestyle list
lifestyle_list <- read.table("02_Function_analysis_out/FungalTraits_w_rarefied_ASV_table_fungi.txt", header = T, row.names = 1, sep = "\t")
lifestyle_list <- lifestyle_list |>
  rownames_to_column(var = "ASV")
lifestyle_list <- lifestyle_list |>
  dplyr::select(ASV, primary_lifestyle)

# soil metadata
DESIGN <- read.table("Data/soil_metadata.txt", header = T)
DESIGN <- DESIGN |>
  dplyr::select(Sample, Country, Site, Landuse) |>
  mutate(treatment = paste(Site, Landuse, sep = "_"))


### format data
# fungal data
f_percent <- f_ASV/mean(colSums(f_ASV))*100 # convert to percentage
f_percent.t <- t(f_percent) # transpose


### SIMPER analysis
set.seed(123)
# simper
sim_fungi <- with(DESIGN, simper(f_percent.t, Landuse, permutations = 999))
sim_fungi_summary <- summary(sim_fungi)
sim_fungi_summary <- sim_fungi_summary$Natural_Farm

# pick up ASVs whose contriburtion to the difference in land uses
sim_fungi_signif <- sim_fungi_summary |>
  dplyr::filter(p < 0.05)
sim_fungi_signif <- sim_fungi_signif |>
  rownames_to_column(var = "ASV") # make ASV column


### summarize SIMPER result by fungal lifestyle
# combine simper result and fungal lifestyles
sim_signif_lifestyle <- merge(sim_fungi_signif, lifestyle_list, by = "ASV", sort = F, all.x = T)

# convert NA to "undetermined" in lifestyle
sim_signif_lifestyle <- sim_signif_lifestyle |>
  mutate(primary_lifestyle = ifelse(is.na(primary_lifestyle), "unassigned", primary_lifestyle))

# remove ASVs whose contribution = 0
sim_signif_lifestyle <- sim_signif_lifestyle |>
  dplyr::filter(average != 0)

# add a bigger category (same way as NMDS plot)
sim_signif_lifestyle <- sim_signif_lifestyle |>
  mutate(lifestyle = case_when(
    grepl("patho", primary_lifestyle, ignore.case = TRUE) ~ "Pathotroph",
    grepl("saprotroph", primary_lifestyle, ignore.case = TRUE) ~ "Saprotroph",
    grepl("ectomycorrhizal", primary_lifestyle, ignore.case = TRUE) ~ "EcM",
    TRUE ~ "Others"
  ))

# sum by primary lifestyle
summed_sim_life <- sim_signif_lifestyle |>
  group_by(primary_lifestyle) |>
  summarise(
    total_ave_contrib_percent = sum(average, na.rm = TRUE) * 100,
    total_ave_relabun_nat = sum(ava, na.rm = TRUE),
    total_ave_relabun_farm = sum(avb, na.rm = TRUE)
  )

# calculate log2 fold change from natural to farm
summed_sim_life <- summed_sim_life |>
  mutate(log2FC = log2((total_ave_relabun_farm+1e-4)/(total_ave_relabun_nat+1e-4)))

# categorize lifestyles by significance (contribution > 1% & over 1-fold change)
summed_sim_life <- summed_sim_life |>
  mutate(significant = case_when(
    total_ave_contrib_percent > 1 & log2FC > 0.5 ~ "Farm",
    total_ave_contrib_percent > 1 & log2FC < -0.5 ~ "Natural",
    TRUE ~ "Not significant"
  ))


### save data
dir.create("11_SIMPER_out")
write.table(sim_fungi_summary, "11_SIMPER_out/simper_allASV.txt", row.names = T, quote = F, sep = "\t")
write.table(sim_signif_lifestyle, "11_SIMPER_out/simper_signif_lifestyle.txt", row.names = F, quote = F, sep = "\t")
write.table(summed_sim_life, "11_SIMPER_out/simper_aggregated_signif_lifestyle.txt", row.names = F, quote = F, sep = "\t")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/11_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

