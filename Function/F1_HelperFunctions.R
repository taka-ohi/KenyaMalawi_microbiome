####
#### R script for Ohigashi et al (2024)
#### 
#### Collection of helper functions for descriptive statistics
#### 2024.06.25 Ohigashi
#### 

### required packages
library(dplyr); packageVersion("dplyr")
library(agricolae); packageVersion("agricolae")
library(vegan); packageVersion("vegan")



### function performing 2way ANOVA for multiple variavles
Do2wayANOVA <- function(data, factor1, factor2, showingstyle = "pval", transformation = "sqrt", log_offset = 1, start_col = 1, end_col = ncol(data)) {
  factor1_col <- data[[factor1]]
  factor2_col <- data[[factor2]]
  # Create a 3x1 matrix with initial values of 0
  p_summary <- matrix(rep(0, 3), nrow = 3, ncol = 1)
  # Set row names and column name for the matrix
  row.names(p_summary) <- c(factor1, factor2, paste(factor1, factor2, sep = " × "))  # Replace these with your actual row names
  colnames(p_summary) <- c("dummy")  # Replace this with your actual column name
  
  for (i in start_col:end_col) {
    # Combine the current column with factor columns for splitting by groups
    current_data <- data[, c(i, which(colnames(data) %in% c(factor1, factor2)))]
    group_data <- split(current_data, interaction(factor1_col, factor2_col))
  
    # Determine if any group fails the normality test
    normal <- TRUE
    for (group in group_data) {
      if (shapiro.test(group[[1]])$p.value < 0.05) {
        normal <- FALSE
        break
      }
    }
    
    # Apply transformation if any group is not normal
    if (!normal) {
      if (transformation == "log") {
        if (min(na.omit(data[,i])) + log_offset <= 0) {
          stop("Data contains values that will result in non-positive values after adding the log offset, which are not suitable for log transformation.")
        }
        data[,i] <- log(data[,i] + log_offset)
      } else if (transformation == "sqrt") {
        if (min(na.omit(data[,i])) < 0) {
          data[,i] <- sqrt(data[,i] - min(na.omit(data[,i])))
        } else {
          data[,i] <- sqrt(data[,i])
        }
      } else {
        stop("Invalid transformation type. Use 'sqrt' or 'log'.")
      }
    }
    
    # Model fitting and ANOVA
    formula <- as.formula(paste0("data[,i] ~ ", factor1, "*", factor2))
    model <- lm(formula, data = data)
    anova_summary <- anova(model)
    anova_summary2 <- data.frame(anova_summary[1:3,5])
    colnames(anova_summary2) <- colnames(data[i])
    # row.names(anova_summary2) <- row.names(anova_summary)[1:3]
    # row.names(anova_summary2)[3] <- sub(":", " × ", row.names(anova_summary2)[3])
    row.names(anova_summary2) <- row.names(p_summary)
    p_summary <- cbind(p_summary, anova_summary2)
  }
  
  # Format and return results based on the specified style
  row_names <- row.names(p_summary)
  p_summary_res <- p_summary[,-1]
  if (showingstyle == "pval") {
    p_summary_res_df <- as.data.frame(lapply(p_summary_res, function(x) {
      ifelse(x < 0.001, "p < 0.001",
             ifelse(x < 0.01, "p < 0.01",
                    ifelse(x < 0.05, "p < 0.05", "n.s.")))
    }))
  } else if (showingstyle == "asterisk") {
    p_summary_res_df <- as.data.frame(lapply(p_summary_res, function(x) {
      ifelse(x < 0.001, "***",
             ifelse(x < 0.01, "**",
                    ifelse(x < 0.05, "*", "n.s.")))
    }))
  } else {
    stop("Invalid showing style. Use 'pval' or 'asterisk'.")
  }
  row.names(p_summary_res_df) <- row_names
  return(p_summary_res_df)
}



### function performing 2way ANOVA succeeding to a table (gtsummary)
Do2wayANOVA_table <- function(data, factor1, factor2, transformation = "sqrt", log_offset = 1, start_col = 1, end_col = ncol(data)) {
  factor1_col <- data[[factor1]]
  factor2_col <- data[[factor2]]
  
  twoway_results <- data.frame(factor1=0, factor2=0, interaction=0, variable=NA, row_type=NA)
  colnames(twoway_results)[1] <- factor1
  colnames(twoway_results)[2] <- factor2
  for (i in start_col:end_col) {
    # Combine the current column with factor columns for splitting by groups
    current_data <- data[, c(i, which(colnames(data) %in% c(factor1, factor2)))]
    group_data <- split(current_data, interaction(factor1_col, factor2_col))
    
    # Determine if any group fails the normality test
    normal <- TRUE
    for (group in group_data) {
      if (shapiro.test(group[[1]])$p.value < 0.05) {
        normal <- FALSE
        break
      }
    }
    
    # Apply transformation if any group is not normal
    if (!normal) {
      if (transformation == "log") {
        if (min(na.omit(data[,i])) + log_offset <= 0) {
          stop("Data contains values that will result in non-positive values after adding the log offset, which are not suitable for log transformation.")
        }
        data[,i] <- log(data[,i] + log_offset)
      } else if (transformation == "sqrt") {
        if (min(na.omit(data[,i])) < 0) {
          data[,i] <- sqrt(data[,i] - min(na.omit(data[,i])))
        } else {
          data[,i] <- sqrt(data[,i])
        }
      } else {
        stop("Invalid transformation type. Use 'sqrt' or 'log'.")
      }
    }
    
    # Model fitting and ANOVA
    anovasummary <- 
      as.formula(paste0("data[,i] ~ ", factor1, "*", factor2)) |>
      aov(data = data) %>% 
      broom::tidy() %>%
      select(term, p.value) %>%
      filter(complete.cases(.)) %>%
      pivot_wider(names_from = term, values_from = p.value) %>%
      mutate(
        variable = names(data)[i],
        row_type = "label"
      )
    colnames(anovasummary)[3] <- "interaction"
    twoway_results <- bind_rows(twoway_results, anovasummary)
  }
  
  # Format and return results based on the specified style
  twoway_results <- twoway_results[-1,]
  return(twoway_results)
}



### function performing Tukey's test for multiple variables 
DoTukeyTest <- function(data, factor1, factor2=NULL, transformation = "sqrt", log_offset = 1, start_col = 1, end_col = ncol(data)) {
  # get factor columns
  if (is.null(factor2)) {
    factor1_col <- data[[factor1]]
    treat <- factor1_col
  } else {
    factor1_col <- data[[factor1]]
    factor2_col <- data[[factor2]]
    treat <- paste(factor1_col, factor2_col, sep = "_") # create a column that combines the factors
  }
  results <- data.frame(treat=unique(treat))
  
  for (i in start_col:end_col) {
    # Combine the current column with factor columns for splitting by groups
    current_data <- data[, c(i, which(colnames(data) %in% c(factor1, factor2)))]
    group_data <- split(current_data, interaction(factor1_col, factor2_col))
    
    # Determine if any group fails the normality test
    normal <- TRUE
    for (group in group_data) {
      if (shapiro.test(group[[1]])$p.value < 0.05) {
        normal <- FALSE
        break
      }
    }
    
    # Apply transformation if any group is not normal
    if (!normal) {
      if (transformation == "log") {
        if (min(na.omit(data[,i])) + log_offset <= 0) {
          stop("Data contains values that will result in non-positive values after adding the log offset, which are not suitable for log transformation.")
        }
        data[,i] <- log(data[,i] + log_offset)
      } else if (transformation == "sqrt") {
        if (min(na.omit(data[,i])) < 0) {
          data[,i] <- sqrt(data[,i] - min(na.omit(data[,i])))
        } else {
          data[,i] <- sqrt(data[,i])
        }
      } else {
        stop("Invalid transformation type. Use 'sqrt' or 'log'.")
      }
    }
    
    # Model fitting and Tukey HSD test
    # formula <- as.formula(paste0("data[,i] ~ ", factor1, "*", factor2))
    data <- data |> mutate(treat = treat)
    formula <- as.formula(paste0("data[,i] ~ ", "treat"))
    model <- lm(formula, data = data)
    # tukey_result <- glht(model, linfct = mcp(treat = "Tukey"))
    # tukey_cld <- cld(tukey_result)
    mod_aov <- aov(model)
    tukey <- HSD.test(mod_aov, "treat", console = F)
    tukey_group <- tukey[["groups"]]
    tukey_group <- tukey_group |>
      dplyr::select(groups) |>
      mutate(treat = rownames(tukey_group))
    colnames(tukey_group)[1] <- colnames(data)[i]
    
    # Extract the grouping information
    # group_labels <- tukey_cld$mcletters$Letters
    # variable_name <- colnames(data)[i]
    results <- merge(results, tukey_group, by = "treat")
  }
  
  return(results)
}



### function performing t-test for multiple variables
DoTTest <- function(data, factor, showingstyle = "pval", transformation = "sqrt", log_offset = 1, start_col = 1, end_col = ncol(data)) {
  factor_col <- data[[factor]]
  
  # Create a 1x(ncol(data) - start_col + 1) matrix with initial values of 0
  p_summary <- matrix(rep(0, 1 * (end_col - start_col + 1)), nrow = 1)
  
  # Set row names and column names for the matrix
  row.names(p_summary) <- factor
  colnames(p_summary) <- colnames(data)[start_col:end_col]
  
  for (i in start_col:end_col) {
    current_data <- data[, c(i, which(colnames(data) %in% factor))]
    
    # Split data by the factor
    group_data <- split(current_data, factor_col)
    
    # Determine if any group fails the normality test
    normal <- TRUE
    for (group in group_data) {
      if (shapiro.test(group[[1]])$p.value < 0.05) {
        normal <- FALSE
        break
      }
    }
    
    # Apply transformation if any group is not normal
    if (!normal) {
      if (transformation == "log") {
        if (min(na.omit(data[, i])) + log_offset <= 0) {
          stop("Data contains values that will result in non-positive values after adding the log offset, which are not suitable for log transformation.")
        }
        data[, i] <- log(data[, i] + log_offset)
      } else if (transformation == "sqrt") {
        if (min(na.omit(data[, i])) < 0) {
          data[, i] <- sqrt(data[, i] - min(na.omit(data[, i])))
        } else {
          data[, i] <- sqrt(data[, i])
        }
      } else {
        stop("Invalid transformation type. Use 'sqrt' or 'log'.")
      }
    }
    
    # Perform t-test
    ttest_result <- t.test(data[, i] ~ factor_col)
    
    # Store the p-value in the corresponding column
    p_summary[1, colnames(data)[i]] <- ttest_result$p.value
  }
  
  # Format and return results based on the specified style
  p_summary_res_df <- if (showingstyle == "pval") {
    as.data.frame(lapply(as.data.frame(p_summary), function(x) {
      ifelse(x < 0.001, "p < 0.001",
             ifelse(x < 0.01, "p < 0.01",
                    ifelse(x < 0.05, "p < 0.05", "n.s.")))
    }))
  } else if (showingstyle == "asterisk") {
    as.data.frame(lapply(as.data.frame(p_summary), function(x) {
      ifelse(x < 0.001, "***",
             ifelse(x < 0.01, "**",
                    ifelse(x < 0.05, "*", "n.s.")))
    }))
  } else {
    stop("Invalid showing style. Use 'pval' or 'asterisk'.")
  }
  
  row.names(p_summary_res_df) <- factor
  return(p_summary_res_df)
}




### function calculating distance to centroids by the group
DistToCent <- function(mat, method = "bray", group, name) {
  # Calculate the vegdist
  veg_dist <- vegdist(mat, method = method)
  # Perform betadisper
  beta_disp <- betadisper(veg_dist, group)
  # Extract distances to centroids
  dist_to_cent <- as.data.frame(beta_disp$distances)
  # Set the column name
  colnames(dist_to_cent) <- name
  return(dist_to_cent)
}



### function calculating correlations between variables within a category
CorrInCat <- function(data, var1, var2, method = "pearson", category) {
  treatment <- unique(data[[category]]) # ex. Site
  cor_result <- data.frame(r=NA, p.value=NA) # make a dummy data frame
  for (i in 1:length(treatment)){
    sub_df <- data |> dplyr::filter(data[[category]]==treatment[i]) # subset by category
    var1_col <- sub_df[[var1]]
    var2_col <- sub_df[[var2]]
    cor <- cor.test(var1_col, var2_col, method = method) # correlation test
    cor_summary <- cbind(cor$estimate, cor$p.value)
    colnames(cor_summary) <- c("r", "p.value")
    row.names(cor_summary) <- treatment[i]
    cor_result <- rbind(cor_result, cor_summary)
  }
  cor_result <- na.omit(cor_result) # remove NA row
  return(cor_result)
}



### function converting Lat/Lon data (e.g., N32° -> 32)
dms_to_deg <- function(dms_string) {
  parts <- as.numeric(strsplit(gsub("[^0-9.]", " ", dms_string), " ")[[1]])
  deg <- parts[2] + parts[3]/60 + parts[4]/3600
  if (substring(dms_string, 1, 1) %in% c("W", "S")) {
    deg <- -deg
  }
  return(deg)
}

