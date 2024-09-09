
### load library
library(dplyr); packageVersion("dplyr")
# library(openxlsx); packageVersion("openxlsx")
library(tibble); packageVersion("tibble")
library(tidyr); packageVersion("tidyr")
library(stringr); packageVersion("stringr")


### load rarefied fungal traits data
rarefied_fungaltrait <- read.table("02_Function_analysis_out/FungalTraits_w_rarefied_ASV_table_fungi.txt",
                                   header = T, sep = "\t")


### aggregation by primary lifestyles
# aggregate by primary lifestyle
agr_fungaltrait <- rarefied_fungaltrait |>
  group_by(primary_lifestyle) |>
  summarise(across(1:90, sum))

# replace NA with "unassigned"
agr_fungaltrait <- agr_fungaltrait |>
  mutate(primary_lifestyle = replace_na(primary_lifestyle, "unassigned"))


### aggregation by specified characteristics
# categorize the lifestyles (Pathotroph, Saprotroph, ECM)
spc_fungaltrait <- agr_fungaltrait |>
  mutate(category = case_when(
    grepl("patho", primary_lifestyle) ~ "Pathotroph",
    grepl("saprotroph", primary_lifestyle) ~ "Saprotroph",
    primary_lifestyle == "ectomycorrhizal" ~ "ECM",
    TRUE ~ "Others"
  ))

# get a table for the category
big_category <- spc_fungaltrait |>
  group_by(category) |>
  summarise(across(2:91, sum))
big_category <- big_category |>
  dplyr::filter(category != "Others") |> # remove Others
  column_to_rownames(var = "category")
big_category.t <- t(big_category) # transpose

# categorize the lifestyles in the detail (including capacity)
# plant pathogenic capacity
plant_pat_cap <- rarefied_fungaltrait |>
  dplyr::filter(!is.na(Plant_pathogenic_capacity_template)) # pick up the capacity was not NA
plant_pat_cap_table <- data.frame(row.names = colnames(plant_pat_cap)[1:90], plant_pat_cap = colSums(plant_pat_cap[,1:90]))

# animal pathogenic capacity
animal_pat_cap <- rarefied_fungaltrait |>
  dplyr::filter(str_detect(Animal_biotrophic_capacity_template, "parasite")) # pick up the "parasite" capacity 
animal_pat_cap_table <- data.frame(row.names = colnames(animal_pat_cap)[1:90], animal_pat_cap = colSums(animal_pat_cap[,1:90]))

# saprotrophs
saprotrophs <- spc_fungaltrait |>
  dplyr::filter(primary_lifestyle %in% c("soil_saprotroph", "litter_saprotroph", "dung_saprotroph", "wood_saprotroph"))
saprotrophs <- saprotrophs |>
  column_to_rownames(var = "primary_lifestyle") |>
  dplyr::select(-category)
saprotrophs.t <- t(saprotrophs)

# create a table for the all extracted characteristics
functable <- cbind(big_category.t, plant_pat_cap_table, animal_pat_cap_table, saprotrophs.t)
functable <- functable |> rownames_to_column(var = "Sample")


### save data
# write.table(rarefied_fungaltrait, file = "02_Function_analysis_out/FungalTraits_w_rarefied_ASV_table_fungi.txt", row.names = T, quote = F, sep = "\t")
write.table(agr_fungaltrait, file = "00_1_check_chao/FungalTraits_primarilifestyle_aggregated.txt", row.names = F, quote = F, sep = "\t")
write.table(functable, file = "00_1_check_chao/FungalTraits_specific_functions.txt", row.names = F, quote = F, sep = "\t")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/02_4_SessionInfo_%s.txt", substr(Sys.time(), 1, 10))) 

