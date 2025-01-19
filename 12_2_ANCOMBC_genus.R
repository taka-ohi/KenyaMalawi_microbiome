####
#### R script for Ohigashi et al (2024)
#### ANCOM for microbial community (class level)
#### 2024.07.20 written by Ohigashi
#### R 4.3.3
####


### load packages
library(ANCOMBC); packageVersion("ANCOMBC")
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(readr); packageVersion("readr")
library(stringr); packageVersion("stringr")
library(DT); packageVersion("DT")


### load data
# environmental data
sample_sheet <- read.table("Data/soil_metadata.txt", header = T, row.names = 1)
# prokaryotes
b_ASV_table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
# fungi
f_ASV_table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)


### format data and create phyloseq objects
## prokaryotes
# get ASV table
b_ASV_sheet <- b_ASV_table[,1:(ncol(b_ASV_table)-7)]

# get taxonomy table
b_tax_sheet <- b_ASV_table[, (ncol(b_ASV_table)-6):ncol(b_ASV_table)]

# check number of variables
dim(sample_sheet); dim(b_ASV_sheet); dim(b_tax_sheet)
all(colnames(b_ASV_sheet) == rownames(sample_sheet))
all(rownames(b_ASV_sheet) == rownames(b_tax_sheet))

# make phyloseq object
b_ASV_ps <- phyloseq(otu_table(b_ASV_sheet, taxa_are_rows = TRUE),
                  sample_data(sample_sheet),
                  tax_table(as.matrix(b_tax_sheet)))

## fungi
# get ASV table
f_ASV_sheet <- f_ASV_table[,1:(ncol(f_ASV_table)-7)]

# get taxonomy table
f_tax_sheet <- f_ASV_table[, (ncol(f_ASV_table)-6):ncol(f_ASV_table)]

# check number of variables
dim(sample_sheet); dim(f_ASV_sheet); dim(f_tax_sheet)
all(colnames(f_ASV_sheet) == rownames(sample_sheet))
all(rownames(f_ASV_sheet) == rownames(f_tax_sheet))

# make phyloseq object
f_ASV_ps <- phyloseq(otu_table(f_ASV_sheet, taxa_are_rows = TRUE),
                     sample_data(sample_sheet),
                     tax_table(as.matrix(f_tax_sheet)))


### ANCOMBC
set.seed(123)
## 1. prokaryotes
# ancombc
b_out.g <- ancombc(data = NULL, assay_name = NULL, 
                tax_level = "Genus", phyloseq = b_ASV_ps, 
                formula = "Landuse", 
                p_adj_method = "BH", prv_cut = 0.1, lib_cut = 1000,
                group = "Landuse", struc_zero = FALSE, neg_lb = TRUE, tol = 1e-5,
                max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                n_cl = 1, verbose = TRUE)

# get result
b_res <- b_out.g$res

# different abundance table
b_diff_tbl <- b_res$diff_abn

# LFC
b_lfc = data.frame(b_res$lfc[, -1] * b_res$diff_abn[, -1], check.names = FALSE) |>
  mutate(taxon = b_res$diff_abn$taxon) |>
  select(taxon, everything())

# SE
b_se = data.frame(b_res$se[, -1] * b_res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon = b_res$diff_abn$taxon) %>%
  select(taxon, everything())
colnames(b_se)[3] = "LanduseSE"

# create a table for visualization
df_b_land <- b_lfc |> 
  left_join(b_se, by = "taxon") |>
  transmute(taxon, LanduseNatural, LanduseSE) |>
  filter(LanduseNatural != 0) |> 
  arrange(desc(LanduseNatural)) |>
  mutate(direct = ifelse((-LanduseNatural) > 0, "Large in Farm", "Large in Natural"))
df_b_land$taxon = factor(df_b_land$taxon, levels = df_b_land$taxon)
df_b_land$direct = factor(df_b_land$direct, 
                            levels = c("Large in Farm", "Large in Natural"))

# filter the ones classified in upper taxonomy
df_b_land_genus <- df_b_land |>
  filter(str_starts(taxon, "Genus"))

# record genus levels
b_sig_genus <- df_b_land_genus$taxon
b_sig_genus <- sub("Genus:", "", b_sig_genus)

# delete character "Genus:" from the taxon but set the level above.
df_b_land_genus$taxon <- sub("Genus:", "", df_b_land_genus$taxon)
df_b_land_genus$taxon <- factor(df_b_land_genus$taxon, levels = b_sig_genus)


## 2. fungi
# ancombc
f_out.g <- ancombc(data = NULL, assay_name = NULL, 
                tax_level = "Genus", phyloseq = f_ASV_ps, 
                formula = "Landuse", 
                p_adj_method = "BH", prv_cut = 0.1, lib_cut = 1000,
                group = "Landuse", struc_zero = FALSE, neg_lb = TRUE, tol = 1e-5,
                max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                n_cl = 1, verbose = TRUE)

# get result
f_res <- f_out.g$res

# different abundance table
f_diff_tbl <- f_res$diff_abn

# LFC
f_lfc <- data.frame(f_res$lfc[, -1] * f_res$diff_abn[, -1], check.names = FALSE) |>
  mutate(taxon = f_res$diff_abn$taxon) |>
  select(taxon, everything())

# SE
f_se <- data.frame(f_res$se[, -1] * f_res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon = f_res$diff_abn$taxon) %>%
  select(taxon, everything())
colnames(f_se)[3] = "LanduseSE"

# create a table for visualization
df_f_land <- f_lfc |> 
  left_join(f_se, by = "taxon") |>
  transmute(taxon, LanduseNatural, LanduseSE) |>
  filter(LanduseNatural != 0) |> 
  arrange(desc(LanduseNatural)) |>
  mutate(direct = ifelse((-LanduseNatural) > 0, "Large in Farm", "Large in Natural"))
df_f_land$taxon = factor(df_f_land$taxon, levels = df_f_land$taxon)
df_f_land$direct = factor(df_f_land$direct, 
                          levels = c("Large in Farm", "Large in Natural"))

# filter the ones classified at upper taxonomy
df_f_land_genus <- df_f_land |>
  filter(str_starts(taxon, "Genus"))

# record genus levels
f_sig_genus <- df_f_land_genus$taxon
f_sig_genus <- sub("Genus:", "", f_sig_genus)

# delete character "genus:" from the taxon but set the level above.
df_f_land_genus$taxon <- sub("Genus:", "", df_f_land_genus$taxon)
df_f_land_genus$taxon <- factor(df_f_land_genus$taxon, levels = f_sig_genus)


### visualization
# prokaryotes
p_b_land.g <- ggplot(data = df_b_land_genus,
                  aes(x = taxon, y = -LanduseNatural, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = -LanduseNatural - LanduseSE, ymax = -LanduseNatural + LanduseSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes in Farm compared to Natural") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank()#,
        # axis.text.x = element_text(angle = 60, hjust = 1)
        ) +
  coord_flip()
print(p_b_land.g)

# fungi
p_f_land.g <- ggplot(data = df_f_land_genus,
                   aes(x = taxon, y = -LanduseNatural, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = -LanduseNatural - LanduseSE, ymax = -LanduseNatural + LanduseSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes in Farm compared to Natural") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank()#,
        # axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  coord_flip()
print(p_f_land.g)


### save data
# dir.create("12_ANCOMBC_out")
saveRDS(b_out.g, "12_ANCOMBC_out/ANCOM_prok_genus.obj")
saveRDS(f_out.g, "12_ANCOMBC_out/ANCOM_fungi_genus.obj")
write.csv(df_b_land_genus, "12_ANCOMBC_out/ANCOM_landuse_prokgenus.csv", quote = F, row.names = F)
write.csv(df_f_land_genus, "12_ANCOMBC_out/ANCOM_landuse_fungigenus.csv", quote = F, row.names = F)
ggsave(filename = "12_ANCOMBC_out/ANCOM_landuse_prokgenus.png", plot = p_b_land.g, width = 6)
ggsave(filename = "12_ANCOMBC_out/ANCOM_landuse_fungigenus.png", plot = p_f_land.g, width = 6)
saveRDS(p_b_land.g, "12_ANCOMBC_out/ANCOM_landuse_prokgenus_plot.obj")
saveRDS(p_f_land.g, "12_ANCOMBC_out/ANCOM_landuse_fungigenus_plot.obj")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/12_SessionInfo_%s.txt", substr(Sys.time(), 1, 10))) 
