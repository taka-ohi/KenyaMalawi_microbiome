

### load packages and functions
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(ggExtra); packageVersion("ggExtra")
library(cowplot); packageVersion("cowplot")


### load data
# for within-site distance
ws_df <- read.table("05_distance_to_centroid_out/distance_to_centroids_withinsite.txt", header = T, sep = "\t")
ws_df_re <- read.table("00_1_check_chao/distance_to_centroids_withinsite_re.txt", header = T, sep = "\t")
ws_df <- ws_df |>
  mutate(proktaxa_within = ws_df_re$proktaxa_within,
         fungitaxa_within = ws_df_re$fungitaxa_within,
         fungilife_within = ws_df_re$fungilife_within
         )

# for across-site distance
as_df <- read.table("05_distance_to_centroid_out/distance_to_centroids_acrosssite.txt", header = T, sep = "\t")
as_df_re <- read.table("00_1_check_chao/distance_to_centroids_acrosssite_re.txt", header = T, sep = "\t")
as_df <- as_df |>
  mutate(proktaxa_across = as_df_re$proktaxa_across,
         fungitaxa_across = as_df_re$fungitaxa_across,
         fungilife_across = as_df_re$fungilife_across
  )


# environmental data
envdata <- read.csv("Data/soil_metadata.txt", header = T, sep = "\t")
envcategory <- envdata[,1:4]

# combine distance data and envcategory data
# relocate rownames to "Sample" column
ws_df <- ws_df |> rownames_to_column(var = "Sample") 
as_df <- as_df |> rownames_to_column(var = "Sample")
# merge envcategory and distance
ws.cat_df <- merge(envcategory, ws_df, by = "Sample", sort = F)
as.cat_df <- merge(envcategory, as_df, by = "Sample", sort = F)


###  summarize the pairs whose correlations are checked
# make a data frame for x and y variable names (within-site)
ws.pairs <- data.frame(pair = seq(1, 10),
                       x = c(rep("env_within", 6), rep("proktaxa_within", 3), "fungitaxa_within"),
                       y = c(names(ws.cat_df)[5:10], names(ws.cat_df)[7:9], "fungilife_within")
)

# add a column for the path to correlation stats file. use this later.
ws.pairs <- ws.pairs |>
  mutate(
    x_transformed = case_when(
      grepl("env", x) ~ "env",
      grepl("fungilife", x) ~ "fungilife",
      grepl("cfunc", x) ~ "funcc",
      grepl("nfunc", x) ~ "funcn",
      grepl("allfunc", x) ~ "funcall",
      grepl("proktaxa", x) ~ "taxaprok",
      grepl("fungitaxa", x) ~ "taxafungi",
      TRUE ~ x
    ),
    y_transformed = case_when(
      grepl("env", y) ~ "env",
      grepl("fungilife", y) ~ "fungilife",
      grepl("cfunc", y) ~ "funcc",
      grepl("nfunc", y) ~ "funcn",
      grepl("allfunc", y) ~ "funcall",
      grepl("proktaxa", y) ~ "taxaprok",
      grepl("fungitaxa", y) ~ "taxafungi",
      TRUE ~ y
    ),
    filename = paste0(x_transformed, "_", y_transformed)
  ) |>
  select(-x_transformed, -y_transformed)

# make a data frame for across-site heterogeneity
as.pairs <- ws.pairs |>
  mutate(across(everything(), ~ sub("within", "across", .))) # utilize within-site one


### create plots
## 1. correlations between within-site heterogeneity
# set title names
ws.plottitles <- c("Environment × Prokaryotes (within-site)",
                   "Environment × Fungi (within-site)",
                   "Environment × All prokaryotic functions (within-site)",
                   "Environment × C cycle (within-site)",
                   "Environment × N cycle (within-site)",
                   "Environment × Fungal lifestyles (within-site)",
                   "Prokaryotes × All prokaryotic functions (within-site)",
                   "Prokaryotes × C cycle (within-site)",
                   "Prokaryotes × N cycle (within-site)",
                   "Fungi × Fungal lifestyles (within-site)"
)

# set an order for the land use
ws.cat_df$Landuse <- factor(ws.cat_df$Landuse, levels = c("Natural", "Farm"))

# looping for creating plots
ws_plots <- list()
for (i in 1:10){
  # set names for x and y axes
  pair_vars <- ws.pairs$filename[i]
  x_name0 <- strsplit(pair_vars, "_")[[1]][1]
  y_name0 <- strsplit(pair_vars, "_")[[1]][2]
  
  x_name_tr <- case_when(
    x_name0 == "env" ~ "environment",
    grepl("func", x_name0) ~ "functions",
    x_name0 == "taxaprok" ~ "prokaryotic communities",
    x_name0 == "taxafungi" ~ "fungal communities",
    x_name0 == "fungilife" ~ "fungal lifestyles",
    TRUE ~ x_name0
  )
  
  y_name_tr <- case_when(
    y_name0 == "env" ~ "environment",
    grepl("func", y_name0) ~ "functions",
    y_name0 == "taxaprok" ~ "prokaryotic communities",
    y_name0 == "taxafungi" ~ "fungal communities",
    y_name0 == "fungilife" ~ "fungal lifestyles",
    TRUE ~ y_name0
  )
  
  x_name <- sprintf("Distance to centeroids of %s", x_name_tr)
  y_name <- sprintf("Distance to centeroids of %s", y_name_tr)
  
  # obtain the result of correlation test
  cor_res <- read.csv(sprintf("06_correlation_among_heterogeneity_out/cor_%s_within.csv", ws.pairs$filename[i]))
  cor_res <- cor_res |>
    mutate(pval = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01 ~ "p < 0.01",
      p.value < 0.05 ~ "p < 0.05",
      TRUE ~ "n.s."
    ))
  
  cor_anno <- "Site\n"
  for (j in 1:nrow(cor_res)) {
    cor_anno0 <- paste0(cor_res[j, 1], ": ", "r = ", round(cor_res[j, 2], 2), ", ", cor_res[j, 4], "\n")
    cor_anno0 <- sub("-", "\U2212", cor_anno0)
    cor_anno <- paste0(cor_anno, cor_anno0)
  }
  
  # make a list of sites that exhibit a correlation between the pairs
  sig_site <- cor_res |>
    filter(pval != "n.s." & X != "All") |>
    pull(X)
  
  # create plot
  ws_g <- ggplot(ws.cat_df,
                 aes_string(x = ws.pairs$x[i], y = ws.pairs$y[i], color = "Site", shape = "Landuse")) +
    geom_point(size=3.5) + 
    scale_shape_manual(values = c(16, 3))+
    labs(x = x_name,
         y = y_name,
         color = "Site",
         shape = "Land use",
         title = ws.plottitles[i]
    ) + 
    geom_smooth(data = subset(ws.cat_df, Site %in% sig_site),
                method = "lm", se = FALSE, #inherit.aes = FALSE,
                aes(group = Site)) +
    theme_classic()+
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 13, face = "bold", colour = "black"),
          legend.text = element_text(size = 13, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", colour = "black")
    ) +
    annotate("text",
             x = ifelse(ws.pairs$x[i] == "env_within", # just avoiding an overlap with points
                        max(ws.cat_df[[ws.pairs$x[i]]])*0.65,
                        min(ws.cat_df[[ws.pairs$x[i]]])
             ),
             y = max(ws.cat_df[[ws.pairs$y[i]]])*0.83,
             label = cor_anno,
             hjust = "left",
             size = 4
    )
  
  ws_plots[[i]] <- ws_g
  
}

# get a legend for within-site plots
ws_legend <- get_legend(ws_plots[[1]])


## 2. correlations between across-site heterogeneity
# set title names
as.plottitles <- gsub("within", "across", ws.plottitles)

# set an order for the land use
as.cat_df$Landuse <- factor(as.cat_df$Landuse, levels = c("Natural", "Farm"))

# looping for creating plots
as_plots <- list()
for (i in 1:10){
  # set names for x and y axes
  pair_vars <- as.pairs$filename[i]
  x_name0 <- strsplit(pair_vars, "_")[[1]][1]
  y_name0 <- strsplit(pair_vars, "_")[[1]][2]
  
  x_name_tr <- case_when(
    x_name0 == "env" ~ "environment",
    grepl("func", x_name0) ~ "functions",
    x_name0 == "taxaprok" ~ "prokaryotic communities",
    x_name0 == "taxafungi" ~ "fungal communities",
    x_name0 == "fungilife" ~ "fungal lifestyles",
    TRUE ~ x_name0
  )
  
  y_name_tr <- case_when(
    y_name0 == "env" ~ "environment",
    grepl("func", y_name0) ~ "functions",
    y_name0 == "taxaprok" ~ "prokaryotic communities",
    y_name0 == "taxafungi" ~ "fungal communities",
    y_name0 == "fungilife" ~ "fungal lifestyles",
    TRUE ~ y_name0
  )
  
  x_name <- sprintf("Distance to centeroids of %s", x_name_tr)
  y_name <- sprintf("Distance to centeroids of %s", y_name_tr)
  
  # obtain the result of correlation test
  cor_res <- read.csv(sprintf("06_correlation_among_heterogeneity_out/cor_%s_across.csv", as.pairs$filename[i]))
  cor_res <- cor_res |>
    mutate(pval = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01 ~ "p < 0.01",
      p.value < 0.05 ~ "p < 0.05",
      TRUE ~ "n.s."
    ))
  
  cor_anno <- paste0("r = ", round(cor_res[1, 2], 2), ", ", cor_res[1, 4])
  cor_anno<- sub("-", "\U2212", cor_anno)
  
  # create plot
  as_g <- ggplot(as.cat_df,
                 aes_string(x = as.pairs$x[i], y = as.pairs$y[i], color = "Landuse", shape = "Site")) +
    geom_point(size=3.5) + 
    scale_color_manual(values=c("Farm"="tan1", "Natural"="palegreen"))+
    labs(x = x_name,
         y = y_name,
         color = "Land use",
         shape = "Site",
         title = as.plottitles[i]
    ) + 
    theme_linedraw()+
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 13, face = "bold", colour = "black"),
          legend.text = element_text(size = 13, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", colour = "black"),
          legend.position = "none", # to do marginal plot, do not set legend here
          # legend.direction = "horizontal",
          # legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    ) +
    annotate("text",
             x = ifelse(ws.pairs$x[i] == "env_within", # just avoiding an overlap with points
                        max(as.cat_df[[as.pairs$x[i]]])*0.8,
                        min(as.cat_df[[as.pairs$x[i]]])),
             y = max(as.cat_df[[as.pairs$y[i]]]),
             label = cor_anno,
             hjust = "left",
             size = 4
    )
  as_g_fin <- ggMarginal(as_g,
                         groupColour = TRUE,
                         groupFill = TRUE
  )
  
  as_plots[[i]] <- as_g_fin
  
}


# make a plot just to pick the legend up
as_test <- ggplot(as.cat_df,
                  aes_string(x = "proktaxa_across", y = "prok_cfunc_across", color = "Landuse", shape = "Site")) +
  geom_point(size=3.5) + 
  scale_color_manual(values=c("Farm"="tan1", "Natural"="palegreen"))+
  labs(x = "test",
       y = "test",
       color = "Land use",
       shape = "Site",
       title = "test"
  ) + 
  theme_linedraw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 13, face = "bold", colour = "black"),
        legend.text = element_text(size = 13, face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", colour = "black"),
        legend.position = "right", # to do marginal plot, do not set legend here
        # legend.direction = "horizontal",
        # legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )

# get legend for across-site plots
as_legend <- get_legend(as_test) 







### save data
# correlations of within-site distances
for (i in 1:10) {
  # saveRDS(ws_plots[[i]],
  #         file = sprintf("FigCode/FigM_correlation_btw_heterogeneity_out/cor_%s_within.obj",
  #                        ws.pairs$filename[i]))
  ggsave(ws_plots[[i]],
         file = sprintf("00_1_check_chao/cor_%s_within_re.png",
                        ws.pairs$filename[i]),
         bg = "white", width = 8, height = 7)
}

# correlations of across-site distances
for (i in 1:10) {
  # saveRDS(as_plots[[i]],
          # file = sprintf("FigCode/FigM_correlation_btw_heterogeneity_out/cor_%s_across.obj",
          #                as.pairs$filename[i]))
  ggsave(as_plots[[i]],
         file = sprintf("00_1_check_chao/cor_%s_across_re.png",
                        as.pairs$filename[i]),
         bg = "white", width = 8, height = 8)
}

