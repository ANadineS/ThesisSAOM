###
# Function to analyze output
###

### Create 5-number summary of network
NetworkSummary <- function(networklist){
  netwobjects <- lapply(networklist, function(x) 
    network(as.matrix(x), directed = T, matrix.type = "adjacency"))
  
  # Calculate relevant network statistics
  netwstats <- data.frame(
    "Network size" = sapply(netwobjects, network.size),
    "Density" = sapply(netwobjects, gden),
    "Components" = sapply(netwobjects, components),
    "Transitivity" = sapply(netwobjects, function(x) gtrans(x, mode = "graph")),
    "Reciprocity" = sapply(netwobjects, function(x) grecip (x, measure = "edgewise")),
    "Med_Distance" = sapply(netwobjects, Med_Distance)
  )
  
  return(netwstats)
}

### Summarize SAOM parameters & calculate relevant statistics
SummarizeParameters <- function(theta, se, error_combo, realmodel){
  # Select data belonging to a certain error combination
  theta_to_use <- theta[theta$Error_Neg == error_combo[1] & theta$Error_Pos == error_combo[2], 1:9]
  se_to_use <- se[se$Error_Neg == error_combo[1] & se$Error_Pos == error_combo[2], 1:9] 
  
  # True model results
  realtheta <- c(realmodel$rate, realmodel$theta)
  realtheta_se <- c(realmodel$vrate, realmodel$se)
  
  # Calculate the statistics on data belonging to a error combination
  stat <- data.frame(
    Parameter = colnames(theta_to_use)[1:9],
    Error_Neg = rep(error_combo[1], 9),
    Error_Pos = rep(error_combo[2], 9),
    Real_Theta = realtheta,
    Real_SE = realtheta_se,
    Theta = sapply(theta_to_use, mean),
    SE = sapply(se_to_use, mean),
    Bias = colMeans(theta_to_use) - realtheta,
    Rel_Bias = (colMeans(theta_to_use) - realtheta)/realtheta*100,
    Variance = sapply(theta_to_use, var),
    MSE = colMeans((theta_to_use - realtheta)^2),
    Precision_Change = apply(theta_to_use, 2, sd)/realtheta_se,
    Coverage = sapply(1:9, function(i) {
      low <- theta_to_use[,i] - qnorm(0.975)*se_to_use[,i]
      high <- theta_to_use[,i] + qnorm(0.975)*se_to_use[,i]
      mean(realtheta[i] >= low & realtheta[i] <= high)
    }),
    Power = sapply(1:9, function(j){
      wald <- theta_to_use[,j]/se_to_use[,j]
      mean(wald <= qnorm(0.025) | wald >= qnorm(0.975))
    }),
    Percentile_Rank = sapply(1:9, function(i){
      ecdf(theta_to_use[,i])(realtheta[i])
    })
  )
  rownames(stat) <- NULL
  return(stat)
}

### Detect non-convergence 
DetectConvergence <- function(output, error_element, focal_combo, n){
  error_combo <- as.matrix(unlist(unname(error_element[focal_combo,])), byrow=T)
  
  # To save indices of non-converged runs
  notconverged <- numeric()
  
  # To save statistics regarding non convergence
  convergestats <- c(error_combo[1], error_combo[2], rep(0, 30))
  
  # Select data belonging to focal error combination
  tconvmax <- output$tconv.max[output$tconv.max$Error_Neg == error_combo[1] &
                                 output$tconv.max$Error_Pos == error_combo[2], 1] %>% as.numeric()
  
  tconv <- output$tconv[output$tconv$Error_Neg == error_combo[1] &
                          output$tconv$Error_Pos == error_combo[2], 1:9] %>% sapply(., as.numeric)
  theta <- output$theta[output$theta$Error_Neg == error_combo[1] &
                          output$theta$Error_Pos == error_combo[2], 1:9] %>% sapply(., as.numeric)
  se <- output$se[output$se$Error_Neg == error_combo[1] &
                    output$se$Error_Pos == error_combo[2], 1:9] %>% sapply(., as.numeric)
  
  # Check all runs for convergence issues
  for (ind in 1:n) {
    status <- TRUE
    if (is.na(tconvmax[ind])){
      
      # No results due to collinearity issues
      notconverged <- c(notconverged, (focal_combo - 1)*n + ind)
      convergestats[3] <- convergestats[3] + 1
      convergestats[32] <- convergestats[32] + 1
      
    } else{
      
      # No convergence due to max convergence ratio being too high
      if (abs(tconvmax[ind]) > 0.25){
        if (status){
          notconverged <- c(notconverged, (focal_combo - 1)*n + ind)
          convergestats[3] <- convergestats[3] + 1
          status <- FALSE
        }
        
        convergestats[4] <- convergestats[4] + 1
      }
      
      # No convergence due to convergence t-ratio for a specific parameter being too high
      if (any(abs(tconv[ind,]) > 0.1)){
        if (status){
          notconverged <- c(notconverged, (focal_combo - 1)*n + ind)
          convergestats[3] <- convergestats[3] + 1
          status <- FALSE
        }
        
        parameters <- which(abs(tconv[ind,]) > 0.1)
        
        for (i in parameters){
          convergestats[4 + i] <- convergestats[4 + i] + 1
        }
      }
      
      # No convergence due to absolute value of parameter being unreasonably high
      if (any(abs(theta[ind,]) > c(50, rep(10, 8)))){
        if (status){
          notconverged <- c(notconverged, (focal_combo - 1)*n + ind)
          convergestats[3] <- convergestats[3] + 1
          status <- FALSE
        }

        parameters <- which(abs(theta[ind,]) > c(50, rep(10, 8)))
        
        for (i in parameters){
          convergestats[13 + i] <- convergestats[13 + i] + 1
        }
      }
      
      # No convergence due to standard error being unreasonably high
      if (any(abs(se[ind,]) > 10)){
        if (status){
          notconverged <- rbind(notconverged, c(error_combo[1], error_combo[2], ind))
          convergestats[3] <- convergestats[3] + 1
          status <- FALSE
        }

        parameters <- which(abs(se[ind,]) > 10)
        
        for (i in parameters){
          convergestats[22 + i] <- convergestats[22 + i] + 1
        }
      }
    }
  }
  
  output <- list(notconverged, convergestats)
  names(output) <- c("notconverged", "convergestats")
  
  return(output)
}

### Analysis function that executes all analyzing functions above, and combines all results
AnalyzeOutput <- function(raw_output, formatted_output, realmodel, n, error_combo){
  # Convergence results
  convergence_results <- sapply(1:nrow(error_combo), function(x) DetectConvergence(formatted_output, error_combo, x, n),
                                simplify = F)
  noconvergence_indices <- do.call(c, lapply(convergence_results, "[[", "notconverged"))
  convergencestats <- as.data.frame(do.call(rbind, lapply(convergence_results, "[[", "convergestats")))
  
  colnames(convergencestats) <- unname(c("Error_Neg", "Error_Pos", "Tot_NConv", "TconvMax", 
                                         sapply(colnames(formatted_output$theta)[1:9], function(x) paste0("T_", x)),
                                         sapply(colnames(formatted_output$theta)[1:9], function(x) paste0("Par_", x)),
                                         sapply(colnames(formatted_output$theta)[1:9], function(x) paste0("SE_", x)),
                                         "Collin"))
  
  # Create new formatted data without non-converged runs
  filtered_output <- raw_output[-noconvergence_indices]
  output <- FormatOutput(filtered_output)
  
  stats_per_network <- vector("list", nrow(error_combo))
  
  for (i in 1:nrow(error_combo)){
    # Indices belonging to the error-combination - excluding non-converged runs
    indices_all <- ((i-1)*(n)+1):(n*i)
    indices <- indices_all[!indices_all %in% noconvergence_indices]
    
    # Calculate network statistics for each network
    stats_error_combo <- list(
      obs_1 = NetworkSummary(formatted_output$original_networks_1[indices]),
      obs_2 = NetworkSummary(formatted_output$original_networks_2[indices]),
      an_1 = NetworkSummary(lapply(formatted_output$analyzed_networks[indices], function(arr) arr[,,1])),
      an_2 = NetworkSummary(lapply(formatted_output$analyzed_networks[indices], function(arr) arr[,,2]))
    )
    
    stats_per_network[[i]] <- setNames(stats_error_combo, c("obs1", "obs2", "analyzed1", "analyzed2"))
    names(stats_per_network)[i] <- paste0(error_combo[i,1], "_", error_combo[i,2])
  }
  
  # Summarize network statistics per error combination & network type
  summary_networkstats_mean <- lapply(c(1:nrow(error_combo)), function(x) {
    stats <- data.frame(
      Error_Neg = rep(error_combo[x,1], 4),
      Error_Pos = rep(error_combo[x,2], 4),
      Density = unlist(lapply(stats_per_network[[x]], function(y) mean(y$Density))),
      Components = unlist(lapply(stats_per_network[[x]], function(y) mean(y$Components))),
      Transitivity = unlist(lapply(stats_per_network[[x]], function(y) mean(y$Transitivity))),
      Reciprocity = unlist(lapply(stats_per_network[[x]], function(y) mean(y$Reciprocity))),
      Med_Distance = unlist(lapply(stats_per_network[[x]], function(y) mean(y$Med_Distance))),
      Density_sd = unlist(lapply(stats_per_network[[x]], function(y) sd(y$Density))),
      Components_sd = unlist(lapply(stats_per_network[[x]], function(y) sd(y$Components))),
      Transitivity_sd = unlist(lapply(stats_per_network[[x]], function(y) sd(y$Transitivity))),
      Reciprocity_sd = unlist(lapply(stats_per_network[[x]], function(y) sd(y$Reciprocity))),
      Med_Distance_sd = unlist(lapply(stats_per_network[[x]], function(y) sd(y$Med_Distance)))
    )
  }) %>% do.call(rbind, .)
  
  # Summarize network change
  hamming <- output$hamming %>%
    group_by(Error_Neg, Error_Pos) %>%
    summarize(
      ham_mean = mean(hamming),
      ham_sd = sd(hamming))
  
  jaccard <- output$jaccard %>%
    group_by(Error_Neg, Error_Pos) %>%
    summarize(
      jac_mean = mean(jaccard),
      jac_sd = sd(jaccard))
  
  change_statistics <- left_join(hamming, jaccard, by = c("Error_Neg", "Error_Pos"))
  
  # Summarize SAOM results per error combination
  summary_theta_mean <- apply(error_combo, 1, function(x) 
    SummarizeParameters(output$theta, output$se, x, realmodel)) %>%
    do.call(rbind, .)
  
  # Combine results
  results <- list(stats_per_network, summary_networkstats_mean, summary_theta_mean, 
                  change_statistics, noconvergence_indices, convergencestats, output)
  
  names(results) <- c("stats_per_network", "summary_networkstats", "summary_theta", 
                      "change_statistics", "noconvergence_indices", "convergencestats", "filtered_output")
  return(results)
}
