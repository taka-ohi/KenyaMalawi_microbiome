
### load packages and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")
library(ggsignif); packageVersion("ggsignif")
library(cowplot); packageVersion("cowplot")
source("Function/F1_HelperFunctions.R")
source("Function/F2_HelperFunctions_for_Visualization.R")

### load data (created in 05_distance_to_centroid.R)
# for within-site distance
ws_df <- read.table("00_1_check_chao/distance_to_centroids_withinsite_re.txt", header = T, sep = "\t")
# for across-site distance
as_df <- read.table("00_1_check_chao/distance_to_centroids_acrosssite_re.txt", header = T, sep = "\t")

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


### 1. statistics for each variable
# (also permutation and Tukey HSD tests were conducted for verification. See "XXXX.R")

## 1-1. within-site distances
# two-way ANOVA for within-site distances
ws_anova <- Do2wayANOVA(data = ws.cat_df,
                        factor1 = "Site",
                        factor2 = "Landuse",
                        showingstyle = "asterisk",
                        transformation = "sqrt",
                        start_col = 5)

# tukey test for within-site distances
ws_Tukey <- DoTukeyTest(data = ws.cat_df,
                        factor1 = "Site",
                        factor2 = "Landuse",
                        transformation = "sqrt",
                        start_col = 5
)
# create Site and Landuse columns for later use
ws_Tukey <- ws_Tukey |>
  separate(treat, into = c("Site", "Landuse"), sep = "_", remove = F)

# calculation of Q3 for each treatment
ws.cat_df <- ws.cat_df |> mutate(treat = paste(Site, Landuse, sep = "_"))
ws.Q3 <- ws.cat_df |>
  group_by(treat) |>
  summarise(across(proktaxa_within:fungilife_within, ~ { # changed the range
    Q1 <- quantile(.x, 0.25, na.rm = TRUE)
    Q3 <- quantile(.x, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    non_outliers <- .x[.x >= (Q1 - 1.5 * IQR) & .x <= (Q3 + 1.5 * IQR)]
    max(non_outliers, na.rm = TRUE)
  }, .names = "{col}.nonout_max"))

# create Site and Landuse columns for later use
ws.Q3 <- ws.Q3 |>
  separate(treat, into = c("Site", "Landuse"), sep = "_", remove = F)

# combine Tukey and Q3 data
ws_Tukey.Q3 <- merge(ws_Tukey, ws.Q3, by = c("treat", "Site", "Landuse"), sort = F)


## 1-2. across-site distances
# t test for across-site distances
as_t.test <- DoTTest(data = as.cat_df,
                     factor = "Landuse",
                     showingstyle = "asterisk",
                     transformation = "sqrt",
                     start_col = 5)



### 2. create plots
## 2-1. within-site distances
# get variable names
varnames <- names(ws.cat_df)[5:7] # changed
# set titles for each plot
plottitles <- setNames(c("Prokaryotes (within-site)", 
                         "Fungi (within-site)", 
                         "Fungal lifestyles (within-site)"),
                       varnames)

# set order for factors
ws.cat_df$Landuse <- factor(ws.cat_df$Landuse, levels = c("Natural", "Farm"))
ws_Tukey.Q3$Landuse <- factor(ws_Tukey.Q3$Landuse, levels = c("Natural", "Farm"))

# looping to make plot objects
ws_plots <- list()
for (i in varnames) {
  twoway_res <- paste0("\n", rownames(ws_anova)[1], ": ", ws_anova[[i]][1],
                       "\n", rownames(ws_anova)[2], ": ", ws_anova[[i]][2],
                       "\n", rownames(ws_anova)[3], ": ", ws_anova[[i]][3]
  )
  yRoof <- max(na.omit(ws.cat_df[[i]]))*1.2
  p <- ggplot(ws.cat_df, aes_string(x = "Site", y = i,fill = "Landuse")) +
    geom_boxplot() +
    labs(
      title = plottitles[[i]],
      y="Distance to centroids",
      x="Site",
      fill="Land use")+
    theme_classic()+
    scale_fill_manual(values = c("palegreen", "tan1"))+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.title = element_text(size = 15, face = "bold", colour = "black"),
      legend.text = element_text(size = 15, face = "bold", colour = "black"),
      plot.title = element_text(size=15, face = "bold", colour = "black"),
      axis.text=element_text(size=10, face = "bold", colour = "black"),
      axis.title=element_text(size=15,face="bold", colour = "black"),
      legend.position = "right",
      strip.text.x = element_text(face="bold", size = 12),
      strip.background = element_rect(fill="white", colour="black", size=1)
    )+
    annotate("text", x=Inf, y=yRoof*0.95,
             label=twoway_res,
             hjust="right",
             size = 4)
  
  ifelse(ws_anova[[i]][3] == "n.s.",
         p_fin <- p,
         p_fin <- p +
           geom_text(
             data = ws_Tukey.Q3,
             mapping = aes_string(x = "Site",
                                  y = sprintf("%s*1.1", paste(i, "nonout_max", sep = ".")),
                                  fill = "Landuse",
                                  label = i
             ),
             position = position_dodge(width = 0.75),
             fontface = "bold"
           )
  )
  ws_plots[[i]] <- p_fin
}


## 2-2. across-site distances
# get variable names
varnames_as <- names(as.cat_df)[5:7] # changed
# set titles for each plot
plottitles_as <- setNames(c("Prokaryotes (across-site)", 
                            "Fungi (across-site)",  
                            "Fungal lifestyles (across-site)"),
                          varnames_as)

# set order for factors
as.cat_df$Landuse <- factor(as.cat_df$Landuse, levels = c("Natural", "Farm"))

# make a data frame for annotation (max values and significance for t-test)
as_annotate <- data.frame(
  varnames_as = varnames_as,
  maxs = sapply(varnames_as, function(x) max(as.cat_df[[x]], na.rm = TRUE)),
  sigs = sapply(varnames_as, function(x) as_t.test[[x]][1])
)


# looping to make plot objects
as_plots <- list()

for (i in varnames_as) {
  as_p <- ggplot(as.cat_df, aes_string(x = "Landuse", y = i, fill = "Landuse"
  )) +
    geom_violin(trim = F) +
    geom_boxplot(width=0.2) +
    labs(
      title = plottitles_as[[i]],
      y="Distance to centroids",
      x=NULL)+
    theme_classic()+
    scale_fill_manual(values = c("palegreen", "tan1"))+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      # legend.title = element_text(size = 15, face = "bold", colour = "black"),
      # legend.text = element_text(size = 15, face = "bold", colour = "black"),
      plot.title = element_text(size=15, face = "bold", colour = "black"),
      axis.text=element_text(size=15, face = "bold", colour = "black"),
      axis.title=element_text(size=15,face="bold", colour = "black"),
      legend.position = "none",
      strip.text.x = element_text(face="bold", size = 12),
      strip.background = element_rect(fill="white", colour="black", size=1)
    )
  
  
  anno <- as_annotate |> filter(varnames_as == i)
  
  
  ## add asterisks or "n.s." to the pairs
  as_p_fin <- as_p +
    annotate("text", #dummy
             x=-Inf,
             y = (anno$maxs)*1.25,
             label="",
             hjust="left",
             vjust=1,
             size = 3) +
    geom_signif(
      annotations = anno$sigs,
      y_position=(anno$maxs)*1.2,
      xmin=1.0,
      xmax=2.0,
      tip_length=0.01,
      textsize=7,
      size=1,
      color="black"
    )
  
  as_plots[[i]] <- as_p_fin
}


### save data
# save object and image files
# within-site distances
for (i in varnames) {
  # saveRDS(ws_plots[[i]], file = sprintf("FigCode/FigM_dist.to.cent_out/%s_dist.to.cent.obj", i))
  ggsave(ws_plots[[i]], file = sprintf("00_1_check_chao/%s_dist.to.cent_re.png", i),
         bg = "white", width = 8, height = 7)
}

# across-site distances
for (i in varnames_as) {
  # saveRDS(as_plots[[i]], file = sprintf("FigCode/FigM_dist.to.cent_out/%s_dist.to.cent.obj", i))
  ggsave(as_plots[[i]], file = sprintf("00_1_check_chao/%s_dist.to.cent_re.png", i),
         bg = "white", width = 8, height = 7)
}
