### 
# Function to format the output
###

### Prepare raw simulation study results for analysis
FormatOutput <- function(output){
  # Gather each element from al simulation runs and add them together in one data object
  # to make them easier to analyze
  
  results <- list()
  
  # Original & analyzed networks
  results$original_networks_1 <- lapply(output, "[[", "original_networks_1")
  results$original_networks_2 <- lapply(output, "[[", "original_networks_2")
  results$analyzed_networks <- lapply(output, "[[", "analyzed_networks")
  
  # The actual change and analyzed change by the SAOM
  results$actualchange <- lapply(output, "[[", "actualchange")
  results$analyzedchange <- lapply(output, "[[", "original_networks")
  
  # Amount of error in each run
  results$error_neg <- sapply(output, "[[", "error_neg")
  results$error_pos <- sapply(output, "[[", "error_pos")
  results$error_in_change <- sapply(output, "[[", "error_in_change") 
  
  # Full results, as a back-up
  results$results <- lapply(output, "[[", "results")
  
  # All dataframes with results per parameter also include the values of the
  # Error combinations to make analysis per error combination easier.
  
  # Conversion statistics
  results$tconv.max <- sapply(output, "[[", "tconv_max") %>%
    cbind(., sapply(output, "[[", "error_neg"), sapply(output, "[[", "error_pos"), 
          sapply(output, "[[", "error_in_change")) %>% as.data.frame()

  colnames(results$tconv.max) <- c("tconv.max", "Error_Neg", "Error_Pos", "Error_Total")
  
  results$tconv <- as.data.frame(do.call(rbind, lapply(output, "[[", "tconv"))) %>%
    cbind(., sapply(output, "[[", "error_neg"), sapply(output, "[[", "error_pos"), 
          sapply(output, "[[", "error_in_change")) 
  
  colnames(results$tconv) <-  
    c("Rate", "Density", "Reciprocity", "Transitivity-Reciprocity", "3-cycles", "GWESP", 
      "Indegree Popularity", "Indegree Activity", "Homophily-Sex", "Error_Neg", 
      "Error_Pos", "Error_Total")
  
  results$convergence_reruns <- sapply(output, "[[", "rerun_conv")
  
  # Change statistics
  results$hamming <- sapply(output, "[[", "hamming") %>%
    cbind(., sapply(output, "[[", "error_neg"), sapply(output, "[[", "error_pos"), 
          sapply(output, "[[", "error_in_change")) %>% as.data.frame()
  
  results$jaccard <- sapply(output, "[[", "jaccard") %>%
    cbind(., sapply(output, "[[", "error_neg"), sapply(output, "[[", "error_pos"), 
          sapply(output, "[[", "error_in_change")) %>% as.data.frame()
  
  colnames(results$hamming) <- c("hamming", "Error_Neg", "Error_Pos", "Error_Total")
  colnames(results$jaccard) <- c("jaccard", "Error_Neg", "Error_Pos", "Error_Total")
  
  # Parameter values, standard errors and p-values
  results$theta <- as.data.frame(do.call(rbind, lapply(output, "[[", "theta"))) %>%
    cbind(., sapply(output, "[[", "error_neg"), sapply(output, "[[", "error_pos"), 
          sapply(output, "[[", "error_in_change"))
  
  results$se <- as.data.frame(do.call(rbind, lapply(output, "[[", "se"))) %>%
    cbind(., sapply(output, "[[", "error_neg"), sapply(output, "[[", "error_pos"), 
          sapply(output, "[[", "error_in_change"))
  
  colnames(results$theta) <- colnames(results$se) <-  
    c("Rate", "Density", "Reciprocity", "Transitivity-Reciprocity", "3-cycles", "GWESP", 
      "Indegree Popularity", "Indegree Activity", "Homophily-Sex", "Error_Neg", 
      "Error_Pos", "Error_Total")
  
  results$wald <- results$theta[,-c(10:12)]/results$se[,-c(10:12)]
  results$p <- apply(as.matrix(results$wald), c(1,2), function(t) pnorm(abs(t), lower.tail = F)) %>%
    cbind(., sapply(output, "[[", "error_neg"), sapply(output, "[[", "error_pos"), 
          sapply(output, "[[", "error_in_change")) %>%
    as.data.frame()
  
  colnames(results$p)[10:12] <- c("Error_Neg", "Error_Pos", "Error_Total")
  
  return(results)
}

