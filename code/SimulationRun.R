###
# Script for one simulation run
###

### One run of the simulation study
SimulationRun <- function(IntroduceError_X, obs1, obs2, sex, controls, simulation_element,
                          extra_args = list()) {
  # Check if (baseline) neg error probability only is supplied, or also an x0 value for sigmoid
  if (length(simulation_element) > 3){
    error_pos <- extra_args$error_pos <- simulation_element[3]
    extra_args$x0 <- simulation_element[4]
    error_neg <- simulation_element[2]
  } else {
    error_pos <- extra_args$error_pos <- simulation_element[3]
    error_neg <- extra_args$error_neg <- simulation_element[2]
  }
  
  # Variable that indicates whether first run converged
  rerun_convergence <- FALSE
  
  # Introduce measurement error based on supplied error-function
  netw_1 <- do.call(IntroduceError_X, c(list(obs1), extra_args))
  netw_2 <- do.call(IntroduceError_X, c(list(obs2), extra_args))
  
  # Detect differences in the changes between observed networks and erroneous network
  actualchange <- ChangeCoding(obs1, obs2)
  analyzedchange <- ChangeCoding(netw_1, netw_2)
  
  # Quantify total amount of error
  error_in_change <- sum(actualchange != analyzedchange) / (nrow(actualchange)^2 - length(diag(actualchange)))
  
  # Stochastic Actor-Oriented Model & preparations
  friendship <- sienaDependent(array(c(netw_1, netw_2), dim = c(133, 133, 2)))
  glasgowdata <- sienaDataCreate(friendship, sex2)

  effects <- getEffects(glasgowdata)
  effects <- includeEffects(effects, name = "friendship", density, recip, cycle3, gwespFF, transRecTrip, inPop, inAct)
  effects <- includeEffects(effects, name = "friendship", sameX, interaction1 = "sex2")

  results <- siena07(controls, data = glasgowdata, effects = effects)
  
  # Rerun of SAOM with `results` as starting point when first run did not converge
  if (abs(results$tconv.max) > 0.25 | any(abs(results$tconv) > 0.1) ){
    results <- siena07(controls, data = glasgowdata, effects = effects, prevAns = results)
    rerun_convergence <- TRUE
  } 

  # Gather output of simulation run
  output <- list(
    original_networks_1 = obs1,
    original_networks_2 = obs2,
    analyzed_networks = array(c(netw_1, netw_2), dim = c(133, 133, 2)),
    actualchange = actualchange,
    analyzedchange = analyzedchange,
    error_neg = error_neg,
    error_pos = error_pos,
    error_in_change = error_in_change,
    results = results,
    theta = c(results$rate, results$theta),
    se = c(results$vrate, results$se),
    el = simulation_element,
    tconv_max = results$tconv.max,
    tconv = results$tconv,
    rerun_conv = rerun_convergence,
    hamming = Hamming(netw_1, netw_2),
    jaccard = Jaccard(netw_1, netw_2)
  )

  return(output)
}

# Safely execute code on supercomputer, so it does not crash with 1 error (/1 singular model)
SafeSimulation <- function(simulation_element){
  tryCatch({
  
    SimulationRun(IntroduceError_Random, netwmatrices[[1]], second_observations[[simulation_element[1]]], 
                  sex2, controls, simulation_element)
    
  }, error = function(e) {
  
  errormessage <- paste("Error in data: ", paste(simulation_element, collapse = ","), " - ", e$message, "\n")
  cat(errormessage, file = "trycatchlog.txt", append = TRUE)
  
  return(NULL)
  
  })
}

