####
#### R script for Ohigashi et al (2024)
#### statistics for environmental data
#### 2024.06.25 written by Ohigashi
#### R 4.3.3
#### README: be sure that "soil_metadata.txt" is in the "Data" directory

#### load library and functions
source("Function/F1_HelperFunctions.R")
library(dplyr)
library(forcats)
library(agricolae)
library(multcomp)
library(gt)
library(gtExtras)
library(tibble)
library(flextable)

#### load data
env <- read.table("Data/soil_metadata.txt", header = T)
env$Gravimetric.water.content <- env$Gravimetric.water.content*100 # make it percentage

#### statistics
# perform 2-way ANOVA
env2way_summary <- Do2wayANOVA(data = env,
                               factor1 = "Site",
                               factor2 = "Landuse",
                               showingstyle = "asterisk",
                               transformation = "sqrt",
                               start_col = 5
                               )

# perform Tukey's test (especially for the vars that had interaction between site and land use)
envTukey <- DoTukeyTest(data = env,
                        factor1 = "Site",
                        factor2 = "Landuse",
                        transformation = "sqrt",
                        start_col = 5
                        )
# reorder the envTukey dataframe
treat_order <- c("D_Natural", "D_Farm", "E_Natural", "E_Farm", "F_Natural", "F_Farm",
                 "G_Natural", "G_Farm", "H_Natural", "H_Farm")
envTukey <- envTukey |>
  mutate(treat = fct_relevel(treat, treat_order)) |>
  arrange(treat)

# calculate mean & sd for each variable
env_meansd <- env |>
  dplyr::select(-Sample, -Country) |>
  group_by(Site, Landuse) |>
  summarise(across(everything(), list(mean = ~mean(.x, na.rm = TRUE), 
                                      sd = ~sd(.x, na.rm = TRUE)))) |>
  # set the level of land use
  mutate(Landuse = fct_relevel(Landuse, "Natural", "Farm")) |>
  # reorder
  arrange(Site, Landuse)

# format mean & sd to make a table
varnames <- names(env[,5:ncol(env)]) # extract variable names

env_meansd_txt <- env_meansd
for (i in 1:length(varnames)) {
  mean_col <- paste(varnames[i], "mean", sep = "_")
  sd_col <- paste(varnames[i], "sd", sep = "_")
  new_col <- varnames[i]
  env_meansd_txt <- env_meansd_txt |>
    # add columns with mean ± sd
    mutate(!!new_col := paste0(sprintf("%.1f", .data[[mean_col]]), " ± ", sprintf("%.1f", .data[[sd_col]])))
}
env_meansd_txt <- env_meansd_txt |> dplyr::select(Site, Landuse, all_of(varnames))

# annotate letters next to values if there is an interaction between site and land use on the variable
# env_meansd_formatted <- env_meansd_txt
# for (var in varnames) {
#   interaction_value <- env2way_summary[[var]][3]
#   if (interaction_value != "n.s.") {
#     env_meansd_formatted[[var]] <- paste(env_meansd_formatted[[var]], envTukey[[var]], sep = " ")
#   } else {
#     env_meansd_formatted[[var]] <- paste(env_meansd_formatted[[var]], "", sep = "")
#   }
# }
# ↑did not use this part

#### make table
# format table of "env_meansd_txt"
env_meansd_txt <- env_meansd_txt |>
  mutate(Site := paste("Site", Site))
tab_meansd <- env_meansd_txt |>
  gt(groupname_col = "Site",
     rowname_col = "Landuse") |>
  tab_header(
    title = "Table. 1 Soil physicochemical properties"
  ) |>
  # tab_stubhead("Property") |>
  tab_style(
    style = list(cell_text(align = "left")),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = list(cell_text(align = "left")),
    locations = cells_body()
  ) |>
  tab_spanner(
    label = "Shannon Diversity",
    columns = c(shannon_16S, shannon_ITS)
  ) |>
  tab_spanner(
    label = "Quantity (log copies g\U207B\U00B9soil)",
    columns = c(log_copy16S_persoil, log_copyITS_persoil)
  ) |>
  # the variables affected by land use w/o interaction were made bold. Pls refer to "env2way_summary"
  tab_style(
    locations = cells_body(
      columns = c(log_copy16S_persoil, CN_ratio) 
    ),
    cell_text(weight = "bold")
  ) |>
  # for the variables that showed interaction, some pairs were made if they showed differences between land uses.
  tab_style(
    locations = cells_body(
      columns = c(Gravimetric.water.content, Carbon),
      rows = Site == "Site D" | Site == "Site H"
    ),
    cell_text(weight = "bold")
  ) |>
  tab_style(
    locations = cells_body(
      columns = c(pH),
      rows = Site != "Site E"
    ),
    cell_text(weight = "bold")
  ) |>
  tab_style(
    locations = cells_body(
      columns = c(log_copyITS_persoil),
      rows = Site == "Site G"
    ),
    cell_text(weight = "bold")
  ) |>
  tab_style(
    locations = cells_body(
      columns = c(shannon_16S),
      rows = Site == "Site F"
    ),
    cell_text(weight = "bold")
  ) |>
  tab_style(
    locations = cells_body(
      columns = c(shannon_ITS),
      rows = Site == "Site D"
    ),
    cell_text(weight = "bold")
  ) |>
  text_transform(
    locations = cells_column_labels(),
    fn=function(x){
      dplyr::case_when(
        x == "shannon_16S" ~ "16S rRNA",
        x == "shannon_ITS" ~ "ITS",
        x == "log_copy16S_persoil" ~ "16S rRNA",
        x == "log_copyITS_persoil" ~ "ITS",
        x == "Carbon" ~ "Total C (%)",
        x == "Nitrogen" ~ "Total N (%)",
        x == "CN_ratio" ~ "C/N ratio",
        x == "Gravimetric.water.content" ~ "Moisture (%)",
        x == "pH" ~ "pH (H\U2082O)",
        TRUE~as.character(x)
      )
    }
  )
# format for publication
tab_meansd <- tab_meansd %>%
  tab_options(
    table.border.top.width = 0, #タイトルの上の線を消す，
    table.border.bottom.width = 0, #注の下の線を消す
    heading.title.font.size = px(18), #タイトルのフォントサイズをいい感じに
    row_group.border.bottom.width = 0, #セクション名の下の線を消す
    table_body.hlines.width = 0, # tableの中の水平線消す
    stub.border.width = 0, # stub列の中の線を消す
    column_labels.border.top.width = 3, # 変数名の行の線を黒く太く
    column_labels.border.top.color = "black", 
    column_labels.border.bottom.width = 2,
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.color = "black", #テーブルの下線を黒く
    table.width = pct(100), # 程よく幅を広げる（数字で調整）
    table.background.color = "white",
    table.font.names = "Times New Roman",
    row_group.border.top.color = "black", #セクション名の上の線を消す
    row_group.border.top.width = 1 #セクション名の上の線を細く
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title")
  ) %>%
  tab_options(
    heading.title.font.weight = "bold"
  )
print(tab_meansd)

# format table for 2way ANOVA result
env2way_summary <- env2way_summary |>
  rownames_to_column(var = "df")
tab_2way <- gt(env2way_summary, rowname_col = "df") |>
  tab_style(
    style = list(cell_text(align = "left")),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = list(cell_text(align = "left")),
    locations = cells_body()
  ) |>
  tab_spanner(
    label = "Shannon Diversity",
    columns = c(shannon_16S, shannon_ITS)
  ) |>
  tab_spanner(
    label = "Quantity (log copies g\U207B\U00B9soil)",
    columns = c(log_copy16S_persoil, log_copyITS_persoil)
  ) |>
  text_transform(
    locations = cells_column_labels(),
    fn=function(x){
      dplyr::case_when(
        x == "shannon_16S" ~ "16S rRNA",
        x == "shannon_ITS" ~ "ITS",
        x == "log_copy16S_persoil" ~ "16S rRNA",
        x == "log_copyITS_persoil" ~ "ITS",
        x == "Carbon" ~ "Total C (%)",
        x == "Nitrogen" ~ "Total N (%)",
        x == "CN_ratio" ~ "C/N ratio",
        x == "Gravimetric.water.content" ~ "Moisture (%)",
        x == "pH" ~ "pH (H\U2082O)",
        TRUE~as.character(x)
      )
    }
  )
tab_2way <- tab_2way %>%
  tab_options(
    table.border.top.width = 0, #タイトルの上の線を消す，
    table.border.bottom.width = 0, #注の下の線を消す
    heading.title.font.size = px(18), #タイトルのフォントサイズをいい感じに
    row_group.border.bottom.width = 0, #セクション名の下の線を消す
    table_body.hlines.width = 0, # tableの中の水平線消す
    stub.border.width = 0, # stub列の中の線を消す
    column_labels.border.top.width = 3, # 変数名の行の線を黒く太く
    column_labels.border.top.color = "black", 
    column_labels.border.bottom.width = 2,
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.color = "black", #テーブルの下線を黒く
    table.width = pct(100), # 程よく幅を広げる（数字で調整）
    table.background.color = "white",
    table.font.names = "Times New Roman",
    row_group.border.top.color = "black", #セクション名の上の線を消す
    row_group.border.top.width = 1 #セクション名の上の線を細く
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title")
  ) %>%
  tab_options(
    heading.title.font.weight = "bold"
  )
print(tab_2way)

#### save files
dir.create("01_summary_envdata_out")
saveRDS(env2way_summary, file = "01_summary_envdata_out/01_2wayANOVA_envdata.obj")
write.csv(env2way_summary, file = "01_summary_envdata_out/01_2wayANOVA_envdata.csv", quote = F)
saveRDS(envTukey, file = "01_summary_envdata_out/01_Tukey_envdata.obj")
saveRDS(env_meansd, file = "01_summary_envdata_out/01_envdata_meansd.obj")
saveRDS(env_meansd_txt, file = "01_summary_envdata_out/01_envdata_meansd_txt.obj")
saveRDS(tab_meansd, file = "01_summary_envdata_out/01_envdata_meansd_table.obj")
saveRDS(tab_2way, file = "01_summary_envdata_out/01_envdata_2way_table.obj")
tab_meansd %>% 
  gtsave(filename = "01_summary_envdata_out/01_envdata_meansd_table.html")

#### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/01_SessionInfo_%s.txt", substr(Sys.time(), 1, 10))) 
