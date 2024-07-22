####
#### R script for Ohigashi et al (2024)
#### 
#### Collection of helper functions for visualization
#### 2024.07.02 written by Ohigashi
#### 


### function changing p-values to significance level to show in a figure
changep <- function(pval) {
  p <- function(pval) {
    if (pval < 0.001) {
      return("***")
    } else {
      if (pval < 0.01){
        return("**")
      } else {
        if (pval < 0.05) {
          return("*")
        } else {
          return("n.s.")
        }
      }
    }
  }
  res <- NULL
  for (i in 1:nrow(pval)){
    b <- p(pval[i])
    b2 <- paste(trimws(rownames(pval)[i]), b, sep = " : ")
    res <- paste(res, b2, sep = "\n")
  }
  return(res)
}


### functions to make axis values to scientific ones considering the disits
# function to round the results to show with one decimal place
scaleFUN <- function(x) {
  a <- sprintf("%.1f", x)
  a2 <- sub('^-', '\U2212', format(a))
  a2 <- trimws(a2)
  return(a2)
}
# function to round the results to show with the second decimal place
scaleFUN2 <- function(x) {
  a <- sprintf("%.2f", x)
  a2 <- sub('^-', '\U2212', format(a))
  a2 <- trimws(a2)
  return(a2)
}



### import font for PDF
# library(extrafont); packageVersion("extrafont")
# font_import()
# loadfonts(device = "")


