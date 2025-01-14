###
# Script to run simulation study based on input parameters on Habrok cluster
###

setwd("/home3/s3363430/Thesis/Code")

library(snow)
library(parallel)
library(ParallelLogger)

addDefaultFileLogger("ErrorLogs.txt")

# Load relevant functions & data
source("./SimulationRun.R")
source("./HelperFunctions.R")

set.seed(3363430)

second_observations <- readRDS("../Data/Output/SecondObservations.RData")
load("../Data/Input/InputSimulation.RData")
list2env(simulation_input, globalenv())

# Set up cluster & add relevant elements
cl <- getMPIcluster()

clusterExport(cl, c("SafeSimulation", "SimulationRun", "IntroduceError_Random", 
                    "second_observations", "netwmatrices", "sex2", "controls", 
                    "simulation_element", "ChangeCoding", "Hamming", "Jaccard"))

clusterEvalQ(cl, {library(igraph); library(RSiena)})

# Perform simulation study - based on error-combination & second network observation 
# index in simulation_element
random_output <- parApply(cl, as.matrix(simulation_element), 1, SafeSimulation)

# Save results & close cluster
saveRDS(random_output, "../Data/Output/01_Random_Raw_Output.RData")

stopCluster(cl)