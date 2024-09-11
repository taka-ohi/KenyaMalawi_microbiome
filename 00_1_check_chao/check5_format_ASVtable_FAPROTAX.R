
# ついでに、FAPROTAXもやってみる
### load package
library(dplyr)
library(tibble)

### load data
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)

b_ASV.fapro <- b_ASV.table %>%
  mutate(taxonomy = apply(b_ASV.table, 1, function(row) {
    paste0("D_0__", row["Kingdom"], "; ",
           "D_1__", row["Phylum"], "; ",
           "D_2__", row["Class"], "; ",
           "D_3__", row["Order"], "; ",
           "D_4__", row["Family"], "; ",
           "D_5__", row["Genus"], "; ",
           "D_6__", row["Species"]
           )
  }))

b_ASV.fapro <- b_ASV.fapro |>
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species)

b_ASV.fapro <- b_ASV.fapro |> rownames_to_column("#OTU ID") # just followed the instruction (http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Instructions)

write.table(b_ASV.fapro, "00_1_check_chao/OTU_table.tsv", row.names = F, quote = F, sep = "\t")
