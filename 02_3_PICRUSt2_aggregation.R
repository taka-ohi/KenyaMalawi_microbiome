####
#### R script for Ohigashi et al (2024)
#### summarize the picrust2 data for C and N cycling
#### 2024.06.30 written by Ohigashi
#### R 4.3.3
####

### load library and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")

### load data
# picrust full data
picrust_data <- read.table("02_Function_analysis_out/picrust2_out_pipeline_rarefied/KO_metagenome_out/pred_metagenome_unstrat.tsv", header = T, row.names = 1)
# function list
func_database <- read.table("Data/KO_list_Kaiser2016_plus_methanogen.txt", header = T, sep = "\t")


### summarizing
# add a column for KO number
picrust_data <- picrust_data |> rownames_to_column(var = "KO")

# merge picrust count data and KO list for C and N cycling, filtering by the list
merged <- merge(func_database, picrust_data, by = "KO", all.x = T, sort = T)
merged <- merged |>
  select(-KO, -Symbol, -Name, -Cycle, -Pathway) |>
  na.omit()

# sum by function category
summarised_data <- merged %>%
  group_by(Function) %>%
  summarise(across(everything(), sum))

# function names as row names
summarised_data <- summarised_data |> column_to_rownames(var = "Function")

# transpose
summarised_data.t <- t(summarised_data)


### save data
# aggregated CN function
write.csv(summarised_data.t, "02_Function_analysis_out/PICRUSt2_aggregated_CN_function.csv", row.names = T, quote = F)

# non-aggregated CN function
non_agr <- merge(func_database, picrust_data, by = "KO", all.x = T, sort = T)
write.csv(non_agr, "02_Function_analysis_out/PICRUSt2_CN_function_table.csv", row.names = T, quote = F)

# full function
picrust_data_wo_KO <- picrust_data |> dplyr::select(-KO)
picrust_data.t <- t(picrust_data_wo_KO)
write.csv(picrust_data.t, "02_Function_analysis_out/PICRUSt2_full_function.csv", row.names = T, quote = F)


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/02_3_SessionInfo_%s.txt", substr(Sys.time(), 1, 10))) 



