####
#### R script for Ohigashi et al (2024)
#### 
#### Function for RCbray using parallel computation
#### ref https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r


### parallel calculation for RCbray
library(pbapply)  # For pblapply
raup_crick_abundance_parallel = function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE,
                                         reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE,
                                         cl=NULL){
  
  # initial setting 
  if(plot_names_in_col1){
    row.names(spXsite) <- spXsite[,1]
    spXsite <- spXsite[,-1]
  }
  
  n_sites <- nrow(spXsite)
  gamma <- ncol(spXsite)
  
  results <- matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ceiling(spXsite / max(spXsite)) -> spXsite.inc
  occur <- apply(spXsite.inc, MARGIN=2, FUN=sum)
  abundance <- apply(spXsite, MARGIN=2, FUN=sum)
  
  ## define a function to parallel caluculation for RC
  calculate_rc_for_pair <- function(pair) {
    
    null.one <- pair[1]
    null.two <- pair[2]
    
    null_bray_curtis <- numeric(reps)
    for(i in 1:reps) {
      # empty null communities of size gamma (ASVs)
      com1 <- rep(0, gamma) 
      # add observed number of species to com1 randomly
      com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)] <- 1
      # weighting by spesies occurence frequencies
      com1.samp.sp <- sample(which(com1 > 0), (sum(spXsite[null.one, ]) - sum(com1)), replace = TRUE, prob = abundance[which(com1 > 0)])
      com1.samp.sp <- cbind(com1.samp.sp, 1)
      com1.sp.counts <- as.data.frame(tapply(com1.samp.sp[, 2], com1.samp.sp[, 1], FUN = sum))
      colnames(com1.sp.counts) <- 'counts'
      com1.sp.counts$sp <- as.numeric(rownames(com1.sp.counts))
      com1[com1.sp.counts$sp] <- com1[com1.sp.counts$sp] + com1.sp.counts$counts
      
      # same for com2
      com2 <- rep(0, gamma)
      com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)] <- 1
      com2.samp.sp <- sample(which(com2>0), (sum(spXsite[null.two,]) - sum(com2)), replace=TRUE, prob=abundance[which(com2>0)])
      com2.samp.sp <- cbind(com2.samp.sp, 1)
      com2.sp.counts <- as.data.frame(tapply(com2.samp.sp[,2], com2.samp.sp[,1], FUN=sum))
      colnames(com2.sp.counts) <- 'counts'
      com2.sp.counts$sp <- as.numeric(rownames(com2.sp.counts))
      com2[com2.sp.counts$sp] <- com2[com2.sp.counts$sp] + com2.sp.counts$counts
      
      # calculate Bray-Curtis dissimilarity for the null communities
      null.spXsite <- rbind(com1, com2)
      # null_bray_curtis[i] <- distance(null.spXsite, method='bray-curtis')
      null_bray_curtis[i] <- vegdist(null.spXsite, method='bray') # edited
    }
    
    # calculate empirically observed Bray-Curtis dissimilarity
    # obs.bray <- distance(spXsite[c(null.one, null.two), ], method='bray-curtis')
    obs.bray <- vegdist(spXsite[c(null.one, null.two), ], method='bray') # edited
    # how many null observations is the observed value tied with?
    num_exact_matching_in_null <- sum(null_bray_curtis == obs.bray)
    # how many null values are smaller than the observed dissimilarity?
    num_less_than_in_null <- sum(null_bray_curtis < obs.bray)
    rc <- num_less_than_in_null / reps
    
    if (split_ties) {
      rc <- (num_less_than_in_null + (num_exact_matching_in_null) / 2) / reps
    }
    
    if (!classic_metric) {
      rc <- (rc - 0.5) * 2
    }
    
    return(c(null.one, null.two, round(rc, digits=2)))
  }
  
  # perform parallel calculation
  pairs <- t(combn(1:n_sites, 2))  # create all pairs of the samples
  results_parallel <- pblapply(1:nrow(pairs), function(idx) {
    set.seed(123 + idx) # setting different generation method for random values
    calculate_rc_for_pair(pairs[idx, ])
  }, cl = cl)
  
  # put the result into the matrix
  for (res in results_parallel) {
    results[res[2], res[1]] <- res[3]
  }
  
  if(as.distance.matrix) {
    results <- as.dist(results)
  }
  
  return(results)
}


