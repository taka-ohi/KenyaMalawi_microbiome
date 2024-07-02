####
#### R script for Ohigashi et al (2024)
#### FungalTraits analysis and summarization
#### 2024.07.01 written by Ohigashi
#### R 4.3.3
####

### load library
library(dplyr); packageVersion("dplyr")
library(openxlsx); packageVersion("openxlsx")
library(tibble); packageVersion("tibble")
library(tidyr); packageVersion("tidyr")
library(stringr); packageVersion("stringr")


### load fungaltraits database
# get the xlsx file from here (https://link.springer.com/article/10.1007/s13225-020-00466-2#Sec10 )
fungaltraits_DB <- read.xlsx("../../Database/13225_2020_466_MOESM4_ESM.xlsx", sheet = 1)


### annotate fungal traits based on the database
# load my fungi rarefiled ASV table
ASV.table <- read.table(file="01_DADA2_out/rarefied_ASV_table_ITS.txt", header=T)

# divide into ASV table and taxa info
ASV <- ASV.table[,1:(ncol(ASV.table)-7)] 
taxa <- ASV.table[,(ncol(ASV.table)-6):ncol(ASV.table)]

# add ID columns
taxa <- taxa |> rownames_to_column("ASV")
ASV_levels <- factor(taxa$ASV, levels = taxa$ASV) # for later use

# merge taxa and DB
merged_table <- merge(taxa, fungaltraits_DB, by.x = c("Phylum", "Class", "Order", "Family", "Genus"),
                      by.y = c("Phylum", "Class", "Order", "Family", "GENUS"), all.x = TRUE, sort = F)

# arrange the dataframe based on the ASV column
merged_table$ASV <- factor(merged_table$ASV, levels = ASV_levels)
merged_table <- merged_table |> arrange(ASV)

# change the order of columns for use
col_names <- names(merged_table)
merged_table <- merged_table |>
  dplyr::select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, col_names[9:27])

# move the row names to one of the columns in the ASV count table
ASV <- ASV |> rownames_to_column(var = "ASV")

# merge ASV count table and fungal trait info table
rarefied_fungaltrait <- merge(ASV, merged_table, by = "ASV", sort = F)
rarefied_fungaltrait <- rarefied_fungaltrait |> column_to_rownames(var = "ASV")


### aggregation by primary lifestyles
# ref -> https://github.com/Chikae-Tatsumi/UNE/blob/main/1_Data_construction/Fungi03_Aggregate_FungalTratit_UNE.r
# convert the ASV counts to percentage
percent_fungaltrait <- rarefied_fungaltrait
percent_fungaltrait[,1:90] <- percent_fungaltrait[,1:90]/mean(colSums(percent_fungaltrait[,1:90]))*100

# aggregate by primary lifestyle
agr_fungaltrait <- percent_fungaltrait |>
  group_by(primary_lifestyle) |>
  summarise(across(1:90, sum))

# replace NA with "unassigned"
agr_fungaltrait <- agr_fungaltrait |>
  mutate(primary_lifestyle = replace_na(primary_lifestyle, "unassigned"))


### aggregation by specified characteristics
# categorize the lifestyles (Pathotroph, Saprotroph, ECM)
agr_fungaltrait <- agr_fungaltrait |>
  mutate(category = case_when(
    grepl("patho", primary_lifestyle) ~ "Pathotroph",
    grepl("saprotroph", primary_lifestyle) ~ "Saprotroph",
    primary_lifestyle == "ectomycorrhizal" ~ "ECM",
    TRUE ~ "Others"
  ))

# get a table for the category
big_category <- agr_fungaltrait |>
  group_by(category) |>
  summarise(across(2:91, sum))
big_category <- big_category |>
  dplyr::filter(category != "Others") |> # remove Others
  column_to_rownames(var = "category")
big_category.t <- t(big_category) # transpose

# categorize the lifestyles in the detail (including capacity)
# plant pathogenic capacity
plant_pat_cap <- percent_fungaltrait |>
  dplyr::filter(!is.na(Plant_pathogenic_capacity_template)) # pick up the capacity was not NA
plant_pat_cap_table <- data.frame(row.names = colnames(plant_pat_cap)[1:90], plant_pat_cap = colSums(plant_pat_cap[,1:90]))

# animal pathogenic capacity
animal_pat_cap <- percent_fungaltrait |>
  dplyr::filter(str_detect(Animal_biotrophic_capacity_template, "parasite")) # pick up the "parasite" capacity 
animal_pat_cap_table <- data.frame(row.names = colnames(animal_pat_cap)[1:90], animal_pat_cap = colSums(animal_pat_cap[,1:90]))

# saprotrophs
saprotrophs <- agr_fungaltrait |>
  dplyr::filter(primary_lifestyle %in% c("soil_saprotroph", "litter_saprotroph", "dung_saprotroph", "wood_saprotroph"))
saprotrophs <- saprotrophs |>
  column_to_rownames(var = "primary_lifestyle") |>
  dplyr::select(-category)
saprotrophs.t <- t(saprotrophs)

# create a table for the all extracted characteristics
functable <- cbind(big_category.t, plant_pat_cap_table, animal_pat_cap_table, saprotrophs.t)
functable <- functable |> rownames_to_column(var = "Sample")


### save data
write.table(rarefied_fungaltrait, file = "02_Function_analysis_out/FungalTraits_w_rarefied_ASV_table_fungi.txt", row.names = T, quote = F, sep = "\t")
write.table(agr_fungaltrait, file = "02_Function_analysis_out/FungalTraits_primarilifestyle_percent.txt", row.names = T, quote = F, sep = "\t")
write.table(functable, file = "02_Function_analysis_out/FungalTraits_specific_functions.txt", row.names = F, quote = F, sep = "\t")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/02_4_SessionInfo_%s.txt", substr(Sys.time(), 1, 10))) 

