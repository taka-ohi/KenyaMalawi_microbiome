####
#### R script for Ohigashi et al (2025)
#### correlation between soil moisture and RCbray in prokaryotes
#### 2025.02.10 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")


### load data
# bNTI and RCbray results in prokaryotes and fungi
b.RC_prok <- read.csv("08_AssemblyProcess_out/summary_bNTI_RCbray_prok.csv", header = T)
b.RC_fungi <- read.csv("08_AssemblyProcess_out/summary_bNTI_RCbray_fungi.csv", header = T)

# sample data
env <- read.table("Data/soil_metadata.txt", header = T)


### format data
# add a column for soil moisture of sample 1
b.RC_prok$row.moisture <- env$Gravimetric.water.content[match(b.RC_prok$row, env$Sample)]*100
b.RC_fungi$row.moisture <- env$Gravimetric.water.content[match(b.RC_fungi$row, env$Sample)]*100

# add a column for soil moisture of sample 2
b.RC_prok$col.moisture <- env$Gravimetric.water.content[match(b.RC_prok$col, env$Sample)]*100
b.RC_fungi$col.moisture <- env$Gravimetric.water.content[match(b.RC_fungi$col, env$Sample)]*100

# add a column for mean value of soil moisture
b.RC_prok$moisture.pair.mean <- rowMeans(b.RC_prok[, c("row.moisture", "col.moisture")], na.rm = TRUE)
b.RC_fungi$moisture.pair.mean <- rowMeans(b.RC_fungi[, c("row.moisture", "col.moisture")], na.rm = TRUE)

# add a column for difference of soil moisture
b.RC_prok$moisture.pair.diff <- abs(b.RC_prok$row.moisture - b.RC_prok$col.moisture)
b.RC_fungi$moisture.pair.diff <- abs(b.RC_fungi$row.moisture - b.RC_fungi$col.moisture)

# calculate absolute values of bNTI
b.RC_prok$bNTI.abs <- abs(b.RC_prok$bNTI)
b.RC_fungi$bNTI.abs <- abs(b.RC_fungi$bNTI)


### calculate model fitting between bNTI vs moisture in prokaryotes and fungi
# Function to calculate AIC, R², and p-value for linear and quadratic models
model_fit_summary <- function(df, y_var, x_var) {
  # Linear model
  lm_model <- lm(as.formula(paste(y_var, "~", x_var)), data = df)
  lm_AIC <- AIC(lm_model)
  lm_R2 <- summary(lm_model)$r.squared
  lm_p <- summary(lm_model)$coefficients[2, 4]
  
  # Quadratic model
  quad_model <- lm(as.formula(paste(y_var, "~ poly(", x_var, ", 2)")), data = df)
  quad_AIC <- AIC(quad_model)
  quad_R2 <- summary(quad_model)$r.squared
  quad_p <- summary(quad_model)$coefficients[3, 4]  # p-value for the quadratic term
  
  # Combine the results
  result <- data.frame(
    Model = c("Linear", "Quadratic"),
    AIC = c(lm_AIC, quad_AIC),
    R2 = c(lm_R2, quad_R2),
    p.value = c(lm_p, quad_p)
  )
  
  return(result)
}

# prokaryotes
fit_mean_prok <- model_fit_summary(b.RC_prok, "bNTI.abs", "moisture.pair.mean")
fit_diff_prok <- model_fit_summary(b.RC_prok, "bNTI.abs", "moisture.pair.diff")

# fungi
fit_mean_fungi <- model_fit_summary(b.RC_fungi, "bNTI.abs", "moisture.pair.mean")
fit_diff_fungi <- model_fit_summary(b.RC_fungi, "bNTI.abs", "moisture.pair.diff")

# create a list for data frames (values and stats)
fitlist <- list(
  mean_prok = list(b.RC_prok, fit_mean_prok),
  diff_prok = list(b.RC_prok, fit_diff_prok),
  mean_fungi = list(b.RC_fungi, fit_mean_fungi),
  diff_fungi = list(b.RC_fungi, fit_diff_fungi)
)

### plot
# looping to create plots
plots <- list()
for (i in seq_along(fitlist)){
  # get values
  df <- fitlist[[i]][[1]]
  # get stats
  stat.df <- fitlist[[i]][[2]]
  
  # format stat df
  stat.df.ed <- stat.df %>%
    mutate(psig = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01 ~ "p < 0.01",
      p.value < 0.05 ~ "p < 0.05",
      TRUE ~ "n.s."
    )
    )
  # format annotation in plot
  anno <- bquote(atop(
    "Linear: " ~ R^2 == .(round(stat.df.ed[["R2"]][1], 3)) * ", " * .(stat.df.ed[["psig"]][1]), 
    "Quadratic: " ~ R^2 == .(round(stat.df.ed[["R2"]][2], 3)) * ", " * .(stat.df.ed[["psig"]][2])
  ))
  
  # set title and axis title
  title0 <- "βNTI × Soil Moisture"
  if (grepl("prok", names(fitlist)[i])){
    title <- paste("Prokaryotic", title0, sep = " ")
  } else {
    title <- paste("Fungal", title0, sep = " ")
  }
  
  if (grepl("mean", names(fitlist)[i])) {
    xaxis <- "Pair mean of soil moisture (%)"
    xvar <- "moisture.pair.mean"
  } else {
    xaxis <- "Δ Soil moisture (%)"
    xvar <- "moisture.pair.diff"
  }
  
  # plot
  p <- ggplot(df, aes(x = !!sym(xvar), y = bNTI.abs)) +
    geom_point(color = "gray", alpha = 0.5) + 
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "black", se = TRUE, linetype = "dashed") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    labs(x = xaxis, y = "Absolute βNTI", title = title)+
    theme_classic() +
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"),
          legend.text = element_text(size = 12, face = "bold", colour = "black"),
          plot.title = element_text(size = 11, hjust = 0.5, face = "bold", colour = "black"),
          legend.position = "right", # to do marginal plot, do not set legend here
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    ) +
    annotate("text",
             x = max(df %>% pull(!!sym(xvar))),
             y = max(df$bNTI.abs),
             label = anno,
             hjust = "right",
             size = 3)
  
  plots[[names(fitlist)[i]]] <- p
  
}

### save
dir.create("FigCode/FigS_bNTI_cor_out")

for (p in seq_along(plots)) {
  ggsave(sprintf("FigCode/FigS_bNTI_cor_out/%s.png", names(plots)[p]), plot = plots[[names(plots)[p]]],
         width = 8, height = 7, bg = "white")
}

saveRDS(plots, file = "FigCode/FigS_bNTI_cor_out/plot_list.pbj")





