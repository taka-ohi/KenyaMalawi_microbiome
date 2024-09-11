####
#### R script for Ohigashi et al (2024)
#### Partitioning dissimilarity of perokaryotic and gungal communities
#### 2024.07.10 written by Ohigashi
#### 2024.09.11 edited by Ohigashi to add statistics
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
# calculation Sorensen dissimilarity
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


### statistics for dissimilarity and its components (replacement and richness differences) added on 2024.09.11
## preparation
# load data
prok.S.3_df <- read.table("08_partitioning_dissimilarity_out/paritionied_dissimilarity_prok.txt", header = T, sep = "\t")
fungi.S.3_df <- read.table("08_partitioning_dissimilarity_out/paritionied_dissimilarity_fungi.txt", header = T, sep = "\t")


### looping for permutation test
scale_cat <- c("within_site", "diff_site_within_country", "diff_site_diff_country")

# prokaryotes
perm_pvals_prok <- list()
for (i in scale_cat){
  # subset by scale
  subdf <- prok.S.3_df |>
    filter(relation == i, !landuse %in% c("Natural_Farm", "Farm_Natural")) |>
    mutate(Dissim = 1-Similarity) |> # put dissimilarity column
    select(Dissim, Repl, RichDiff, landuse)
  
  # get difference of mean values between landuse pairs (observed data)
  subdf_mean <- subdf |>
    group_by(landuse) |>
    summarise(across(c(Dissim, Repl, RichDiff), mean, na.rm = TRUE))
  obs_difmean_dissim <- subdf_mean$Dissim[subdf_mean$landuse == "Farm_Farm"] - subdf_mean$Dissim[subdf_mean$landuse == "Natural_Natural"]
  obs_difmean_repl <- subdf_mean$Repl[subdf_mean$landuse == "Farm_Farm"] - subdf_mean$Repl[subdf_mean$landuse == "Natural_Natural"]
  obs_difmean_rich <- subdf_mean$RichDiff[subdf_mean$landuse == "Farm_Farm"] - subdf_mean$RichDiff[subdf_mean$landuse == "Natural_Natural"]
  
  
  # permute and get difference of mean values for each variable
  varnames <- names(subdf)[1:3]
  perm_diffs <- list()
  set.seed(123)
  for (j in varnames) {
    perm_diff <- replicate(999999, {
      shuffled <- sample(subdf[[j]])
      mean(shuffled[1:(length(shuffled)/2)]) - 
      mean(shuffled[(length(shuffled)/2 + 1):length(shuffled)])
    })
    perm_diffs[[j]] <- perm_diff
  }
  
  # calculate p-values
  p.dissim <- mean(abs(perm_diffs[["Dissim"]]) >= abs(obs_difmean_dissim))
  p.repl <- mean(abs(perm_diffs[["Repl"]]) >= abs(obs_difmean_repl))
  p.rich <- mean(abs(perm_diffs[["RichDiff"]]) >= abs(obs_difmean_rich))
  
  # make a table
  perm_pvals_prok[[i]] <- data.frame(p.dissim = p.dissim,
                                     p.repl = p.repl,
                                     p.rich = p.rich,
                                     row.names = sprintf("%s_NN_vs_FF", i)
                                     )
}



# fungi
perm_pvals_fungi <- list()
for (i in scale_cat){
  # subset by scale
  subdf <- fungi.S.3_df |>
    filter(relation == i, !landuse %in% c("Natural_Farm", "Farm_Natural")) |>
    mutate(Dissim = 1-Similarity) |> # put dissimilarity column
    select(Dissim, Repl, RichDiff, landuse)
  
  # get difference of mean values between landuse pairs (observed data)
  subdf_mean <- subdf |>
    group_by(landuse) |>
    summarise(across(c(Dissim, Repl, RichDiff), mean, na.rm = TRUE))
  obs_difmean_dissim <- subdf_mean$Dissim[subdf_mean$landuse == "Farm_Farm"] - subdf_mean$Dissim[subdf_mean$landuse == "Natural_Natural"]
  obs_difmean_repl <- subdf_mean$Repl[subdf_mean$landuse == "Farm_Farm"] - subdf_mean$Repl[subdf_mean$landuse == "Natural_Natural"]
  obs_difmean_rich <- subdf_mean$RichDiff[subdf_mean$landuse == "Farm_Farm"] - subdf_mean$RichDiff[subdf_mean$landuse == "Natural_Natural"]
  
  
  # permute and get difference of mean values for each variable
  varnames <- names(subdf)[1:3]
  perm_diffs <- list()
  set.seed(123)
  for (j in varnames) {
    perm_diff <- replicate(999999, {
      shuffled <- sample(subdf[[j]])
      mean(shuffled[1:(length(shuffled)/2)]) - 
        mean(shuffled[(length(shuffled)/2 + 1):length(shuffled)])
    })
    perm_diffs[[j]] <- perm_diff
  }
  
  # calculate p-values
  p.dissim <- mean(abs(perm_diffs[["Dissim"]]) >= abs(obs_difmean_dissim))
  p.repl <- mean(abs(perm_diffs[["Repl"]]) >= abs(obs_difmean_repl))
  p.rich <- mean(abs(perm_diffs[["RichDiff"]]) >= abs(obs_difmean_rich))
  
  # make a table
  perm_pvals_fungi[[i]] <- data.frame(p.dissim = p.dissim,
                                     p.repl = p.repl,
                                     p.rich = p.rich,
                                     row.names = sprintf("%s_NN_vs_FF", i)
  )
}

# create tables
perm_pvals_prok.tbl <- rbind(perm_pvals_prok[[1]], perm_pvals_prok[[2]], perm_pvals_prok[[3]])
perm_pvals_fungi.tbl <- rbind(perm_pvals_fungi[[1]], perm_pvals_fungi[[2]], perm_pvals_fungi[[3]])


### save data
dir.create("08_partitioning_dissimilarity_out")

# data
write.table(prok.S.3_df, "08_partitioning_dissimilarity_out/paritionied_dissimilarity_prok.txt", row.names = F, quote = F, sep = "\t")
write.table(fungi.S.3_df, "08_partitioning_dissimilarity_out/paritionied_dissimilarity_fungi.txt", row.names = F, quote = F, sep = "\t")
# p values
write.table(perm_pvals_prok.tbl, "08_partitioning_dissimilarity_out/part_dissim_permpval_prok.txt", row.names = T, quote = F, sep = "\t")
write.table(perm_pvals_fungi.tbl, "08_partitioning_dissimilarity_out/part_dissim_permpval_fungi.txt", row.names = T, quote = F, sep = "\t")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/08_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

