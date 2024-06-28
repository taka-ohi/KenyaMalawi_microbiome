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
library(gtsummary)
library(tibble)
library(flextable)
library(tidyverse)

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


#### make table
# make a dataframe of p-values that can be succeeded to gt and gtsummary
twoway_results <- Do2wayANOVA_table(data = env, factor1 = "Site", factor2 = "Landuse", 
                                    transformation = "sqrt", start_col = 5)
# format the style
env2 <- env |>
  mutate(Site := paste("Site", Site)) |>
  mutate(Landuse = fct_relevel(Landuse, "Natural", "Farm"))
# summarize the data using gtsummary
tbl <-
  # first build a stratified `tbl_summary()` table to get summary stats by two variables
  env2 |>
  dplyr::select(-Sample, -Country) |>
  tbl_strata(
    strata =  Site,
    .tbl_fun =
      ~.x %>%
      tbl_summary(
        by = Landuse,
        missing = "no",
        statistic = all_continuous() ~ "{mean} ± {sd}",
        digits = everything() ~ 1,
        label = list(
          shannon_16S ~ "16S rRNA (Shannon Diversity)",
          shannon_ITS ~ "ITS (Shannon Diversity)",
          log_copy16S_persoil ~ "16S rRNA (log copies g\U207B\U00B9soil)",
          log_copyITS_persoil ~ "ITS (log copies g\U207B\U00B9soil)",
          Carbon ~ "Total C (%)",
          Nitrogen ~ "Total N (%)",
          CN_ratio ~ "C/N ratio",
          Gravimetric.water.content ~ "Moisture (%)",
          pH ~ "pH (H\U2082O)"
        )
      ) %>%
      modify_header(all_stat_cols() ~ "{level}"
                    )
  ) %>%
  # merge the 2way ANOVA results into tbl_summary table
  modify_table_body(
    ~.x %>%
      left_join(
        twoway_results,
        by = c("variable", "row_type")
      )
  ) %>%
  # by default the new columns are hidden, add a header to unhide them
  modify_header(list(
    Site ~ "Site", 
    Landuse ~ "Landuse", 
    interaction ~ "Site × Landuse"
  )) %>%
  # adding spanning header to analysis results
  modify_spanning_header(c(Site, Landuse, interaction) ~ "Two-way ANOVA p-values") %>%
  # format the p-values with a pvalue formatting function
  modify_fmt_fun(c(Site, Landuse, interaction) ~ style_pvalue) %>%
  # update the footnote to be nicer looking
  modify_footnote(all_stat_cols() ~ "Mean ± SD in each location (n=9).")

tbl_formatted <- 
  tbl |> 
  as_gt()|> # success to gt
  gt::tab_header(
    title = "Table. 1 Soil physicochemical properties"
  ) |>
  tab_options(
    table.border.top.width = 0,
    table.border.bottom.width = 0, 
    heading.title.font.size = px(18), 
    row_group.border.bottom.width = 0, 
    table_body.hlines.width = 0, 
    stub.border.width = 0, 
    column_labels.border.top.width = 3, 
    column_labels.border.top.color = "black", 
    column_labels.border.bottom.width = 2,
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.color = "black", 
    table.width = pct(100), 
    table.background.color = "white",
    table.font.names = "Arial",
    row_group.border.top.color = "black", 
    row_group.border.top.width = 1, 
    heading.title.font.weight = "bold"
  ) %>%
  tab_style(
    style = cell_text(align = "left"),
    locations = cells_title("title")
  ) 
print(tbl_formatted)


#### save files
dir.create("03_summary_envdata_out")
saveRDS(env2way_summary, file = "03_summary_envdata_out/03_2wayANOVA_envdata.obj")
write.csv(env2way_summary, file = "03_summary_envdata_out/03_2wayANOVA_envdata.csv", quote = F)
saveRDS(envTukey, file = "03_summary_envdata_out/03_Tukey_envdata.obj")
write.csv(envTukey, file = "03_summary_envdata_out/03_Tukey_envdata.csv", quote = F)
saveRDS(tbl_formatted, file = "03_summary_envdata_out/03_envdata_meansd_table_w_2way.obj")
tbl_formatted |> 
  gtsave(filename = "03_summary_envdata_out/03_envdata_meansd_table_w_2way.tex")
tbl_formatted |> 
  gtsave(filename = "03_summary_envdata_out/03_envdata_meansd_table_w_2way.docx")
tbl_formatted |> 
  gtsave(filename = "03_summary_envdata_out/03_envdata_meansd_table_w_2way.html")


#### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/03_SessionInfo_%s.txt", substr(Sys.time(), 1, 10))) 

#### reference URL
# https://stackoverflow.com/questions/66835663/summary-table-mean-std-error-with-p-values-for-2-way-anova 
