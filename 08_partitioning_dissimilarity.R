####
#### R script for Ohigashi et al (2024)
#### Partitioning dissimilarity of perokaryotic and gungal communities
#### 2024.07.10 written by Ohigashi
#### R 4.3.3
####


### load packages
library(adespatial); packageVersion("adespatial")


### load data
# prokaryotes
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
b_ASV <- b_ASV.table[,1:(ncol(b_ASV.table)-7)]
b_ASV.t <- t(b_ASV)

# fungi
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)
f_ASV <- f_ASV.table[,1:(ncol(f_ASV.table)-7)]
f_ASV.t <- t(f_ASV)

# environmental data
env_data <- read.table("Data/soil_metadata.txt", header = T)

# prepare to categorize the pairs of samples
sample_names <- env_data$Sample
combinations <- expand.grid(sample_names, sample_names, stringsAsFactors = FALSE) # list all pairs
combinations <- combinations[combinations[,1] < combinations[,2],] # extract pairs (removed the replaced pairs)
comb <- apply(combinations, 1, function(x) paste(x, collapse = "_"))
comb <- sort(comb)


### calculate Sorensen dissimilarity, replacement, and richness difference
# prokaryotes
# calculation of Sorensen similarity
prok.S <- beta.div.comp(b_ASV.t, coef = "S", quant = FALSE)
prok.S.3 <- cbind((1-prok.S$D),
                  prok.S$repl,
                  prok.S$rich
                  )
colnames(prok.S.3) <- c("Similarity", "Repl", "RichDiff")

# categorize the pairs
prok.S.3_df <- as.data.frame(prok.S.3)
prok.S.3_df$comb <- comb

prok.S.3_df$landuse <- paste(env_data$Landuse[match(substr(prok.S.3_df$comb, 1, 3), env_data$Sample)],
                             env_data$Landuse[match(substr(prok.S.3_df$comb, 5, 7), env_data$Sample)],
                             sep = "_"
                             )
prok.S.3_df$site <- paste(env_data$Site[match(substr(prok.S.3_df$comb, 1, 3), env_data$Sample)],
                          env_data$Site[match(substr(prok.S.3_df$comb, 5, 7), env_data$Sample)],
                          sep = "_"
                          )
prok.S.3_df$country <- paste(env_data$Country[match(substr(prok.S.3_df$comb, 1, 3), env_data$Sample)],
                             env_data$Country[match(substr(prok.S.3_df$comb, 5, 7), env_data$Sample)],
                             sep = "_"
                             )
prok.S.3_df <- prok.S.3_df |>
  rowwise() |>
  mutate(relation = case_when(
    strsplit(site, "_")[[1]][1] == strsplit(site, "_")[[1]][2] ~ "within_site",  # if the 1st and 3rd characters are the same.
    strsplit(site, "_")[[1]][1] != strsplit(site, "_")[[1]][2] && strsplit(country, "_")[[1]][1] == strsplit(country, "_")[[1]][2] ~ "diff_site_within_country",
    TRUE ~ "diff_site_diff_country"  # else
  )) |>
  ungroup()


# fungi
fungi.S <- beta.div.comp(f_ASV.t, coef = "S", quant = FALSE)
fungi.S.3 <- cbind((1-fungi.S$D),
                   fungi.S$repl,
                   fungi.S$rich
                   )
colnames(fungi.S.3) <- c("Similarity", "Repl", "RichDiff")

# categorize the pairs
fungi.S.3_df <- as.data.frame(fungi.S.3)
fungi.S.3_df$comb <- comb

fungi.S.3_df$landuse <- paste(env_data$Landuse[match(substr(fungi.S.3_df$comb, 1, 3), env_data$Sample)],
                              env_data$Landuse[match(substr(fungi.S.3_df$comb, 5, 7), env_data$Sample)],
                              sep = "_"
                              )
fungi.S.3_df$site <- paste(env_data$Site[match(substr(fungi.S.3_df$comb, 1, 3), env_data$Sample)],
                           env_data$Site[match(substr(fungi.S.3_df$comb, 5, 7), env_data$Sample)],
                           sep = "_"
                           )
fungi.S.3_df$country <- paste(env_data$Country[match(substr(fungi.S.3_df$comb, 1, 3), env_data$Sample)],
                              env_data$Country[match(substr(fungi.S.3_df$comb, 5, 7), env_data$Sample)],
                              sep = "_"
                              )
fungi.S.3_df <- fungi.S.3_df |>
  rowwise() |>
  mutate(relation = case_when(
    strsplit(site, "_")[[1]][1] == strsplit(site, "_")[[1]][2] ~ "within_site",  # if the 1st and 3rd characters are the same.
    strsplit(site, "_")[[1]][1] != strsplit(site, "_")[[1]][2] && strsplit(country, "_")[[1]][1] == strsplit(country, "_")[[1]][2] ~ "diff_site_within_country",
    TRUE ~ "diff_site_diff_country"  # else
  )) |>
  ungroup()


### save data
dir.create("08_partitioning_dissimilarity_out")
write.table(prok.S.3_df, "08_partitioning_dissimilarity_out/paritionied_dissimilarity_prok.txt", row.names = F, quote = F, sep = "\t")
write.table(fungi.S.3_df, "08_partitioning_dissimilarity_out/paritionied_dissimilarity_fungi.txt", row.names = F, quote = F, sep = "\t")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/08_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

