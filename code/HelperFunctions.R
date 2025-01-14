###
# Various helper functions, used in different parts of the simulation & data formatting process
###

### Error introduction functions
# Random error introduction
IntroduceError_Random <- function(matr, error_neg, error_pos) {
  matr <- apply(matr, c(1, 2), function(cell) {
    if (cell == 1) {
      # Introduce a false negative with probability error_neg
      rbinom(1, 1, 1 - error_neg)
    } else {
      # Introduce a false positive with probability error_pos
      rbinom(1, 1, error_pos)
    }
  })
  
  # Remove any self-loops that may have been introduced
  diag(matr) <- 0
  
  return(matr)
}

# Systematic error introduction, based on tie embeddedness
IntroduceError_Systematic_Embeddedness <- function(matr, k = 4.5, error_pos, x0){
  # Create undirected graph based on matrix
  networkgraph <- graph_from_adjacency_matrix(matr, mode = "max")
  
  # Determining relationship strength using local edge clustering
  clustering_matrix <- matrix(NA, nrow = 133, ncol = 133)
  actor_names <- V(networkgraph)$name  
  rownames(clustering_matrix) <- actor_names
  colnames(clustering_matrix) <- actor_names
  
  for(edge_id in seq_along(E(networkgraph))){
    actor1 <- ends(networkgraph, edge_id)[1]
    actor2 <- ends(networkgraph, edge_id)[2]
    edge_clust_coef <- LocalEdgeClustering(networkgraph, edge_id)
    
    # Only add tie strengths if there is a tie.
    if (matr[actor1, actor2] == 1) {
      clustering_matrix[actor1, actor2] <- edge_clust_coef
    } 
    
    if (matr[actor2, actor1] == 1) {
      clustering_matrix[actor2, actor1] <- edge_clust_coef
    }
  }
  
  # Incorporating error based on relationship strength using sigmoid function
  matr_final <- matrix(apply(expand.grid(1:nrow(matr), 1:ncol(matr)), 1, function(x) {
    i <- x[1]
    j <- x[2]
    
    if (matr[i, j] == 1) {
      rel_strength = clustering_matrix[i, j]
      retention_probability = SigmoidFunction(rel_strength, x0 = x0, k = k)
      
      # Introduce a false negative with probability error_neg
      rbinom(1, 1, retention_probability)
    } else {
      # Introduce a false positive with probability error_pos
      # Not yet incorporated now
      rbinom(1, 1, error_pos)
    }
  }), nrow = nrow(matr), ncol = ncol(matr))
  
  # Remove any self-loops that may have been introduced
  diag(matr_final) <- 0
  return(matr_final)
}

### Functions related to systematic scenario set-up
# Function to calculate the embeddedness of an edge in a given graph
LocalEdgeClustering <- function(graph, edge){
  actor1 <- ends(graph, edge)[1]
  actor2 <- ends(graph, edge)[2]
  
  common_neighbours <- intersect(neighbors(graph, actor1), neighbors(graph, actor2))
  closed_triangles <- length(common_neighbours)
  
  # potential triangles: any triplet in which the addition of one tie from i to j to another actor k
  # creates a closed triangle (any tie configuration).
  
  # Equal to the sum of the outdegree values of both actors - 2 
  # (tie from - to j itself is included in both degrees) and - closed 
  # triangles (otherwise they are counted twice)
  potential_triangles <- degree(graph, actor1) + degree(graph, actor2) - 2 - closed_triangles
  
  if (potential_triangles == 0){
    return(0)
  } else{
    return(closed_triangles/potential_triangles)
  }
}

# Function to calculate output values of a sigmoid function, given midpoint and steepness
SigmoidFunction <- function(x, x0, k){
  retention <- 1/(1 + exp(-k * (x - x0)))
  return(retention)
}

### Other functions
# Function to find the median distance in the largest component(s)
Med_Distance <- function(net) {
  # Find connected components and get the component membership
  comp_info <- component.dist(net)
  membership <- comp_info$membership  
  
  # Size of the largest component
  largest_size <- max(table(membership)) 
  
  # Calculate geodesic distances for the largest component(s)
  all_distances <- unlist(lapply(unique(membership)[table(membership) == largest_size], function(comp) {
    nodes <- which(membership == comp)
    gdist <- geodist(net[nodes, nodes])$gdist
    gdist[!is.infinite(gdist)] # Ignore infinite distances
  }))
  
  median(all_distances, na.rm = TRUE)
}

#Function to code changes between two observations of the same network to quantify
# total change in change due to error
ChangeCoding <- function(matrix1, matrix2){
  changematrix <- matrix(NA, nrow(matrix1), nrow(matrix1))
  
  # Check for each cell in matrix 1 and 2 what kind of change occurred
  for(i in 1:nrow(matrix1)){
    for (j in 1:ncol(matrix1)){
      if (matrix1[i,j] == 0){
        if (matrix2[i,j] == 0){
          
          # Tie does not exist in both matrices
          changematrix[i,j] <- 0
        } else {
          
          # A new tie was created in matrix 2 that did not exist in matrix 1
          changematrix[i,j] <- 1
        } 
      } else {
        if (matrix2[i,j] == 1){
          
          # Tie exists in both matrices
          changematrix[i,j] <- 2
        } else {
          
          # A tie was removed in matrix 2 that existed in matrix 1
          changematrix[i,j] <- 3
        }
      }
    }
  }
  
  return(changematrix)
}

# Hamming distance function
Hamming <- function(net1,net2) {
  tbl <- table(c(0,0,1,1,net1),c(0,1,0,1,net2))-1
  return(tbl[1,2]+tbl[2,1])
}

# Jaccard index function
Jaccard <- function(net1,net2) {
  tbl <- table(c(0,0,1,1,net1),c(0,1,0,1,net2))-1
  return(tbl[2,2]/(tbl[1,2]+tbl[2,1]+tbl[2,2]))
}

# Function to transform edgelists to adjacency matrices
EdgelistToMatrix <- function(edgelist, n = 133){
  
  # Add self-loops to edgelist to make sure isolated actors are also included
  self_loops <- matrix(c(1:133, 1:133, rep(1, 133)), ncol = 3, nrow = 133, byrow = F)
  final_edgelist <- rbind(self_loops, edgelist)
  
  # Transform edgelist to adjacency matrix
  graph <- graph_from_edgelist(as.matrix(final_edgelist[,-3]), directed = T)
  adj_matrix <- as.matrix(as_adjacency_matrix(graph))
  diag(adj_matrix) <- 0 # Remove self-loops
  
  return(adj_matrix)
}

# Function to calculate the proportion of same-sex ties in an adjacency matrix
PropSameSex <- function(network_matrix, sex_vector) {
  # Create an adjacency matrix indicating same-sex ties
  same_sex_matrix <- outer(sex_vector, sex_vector, FUN = "==") * network_matrix
  
  # Calculate total ties and same-sex ties
  total_ties <- sum(network_matrix, na.rm = TRUE)
  same_sex_ties <- sum(same_sex_matrix, na.rm = TRUE)
  
  # Calculate proportion
  if (total_ties > 0) {
    return(same_sex_ties / total_ties)
  } else {
    return(NA)  # Return NA if there are no ties
  }
}
