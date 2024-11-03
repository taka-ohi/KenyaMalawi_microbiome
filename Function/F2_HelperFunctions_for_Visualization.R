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


### function changing p-values to significance (without factor)
changep0 <- function(pval) {
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


colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
            "#A65628", "#F781BF", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
            "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
            "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5",
            "#FFED6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
            "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")


