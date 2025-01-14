###
# Create graphs to analyse parameter outcomes per error combination
###

### Graph creation
CreateGraph <- function(datasummary, data = NULL, variable, error_neg = 0.1, error_pos = 0,
                        dimensions_x, dimensions_y){
  
  # Select relevant data for selected variable and error combinations
  data <- subset(data, Error_Neg %in% error_neg & Error_Pos == error_pos)
  realvalue <- datasummary[datasummary$Parameter == variable, "Real_Theta"]
  realse <- datasummary[datasummary$Parameter == variable, "Real_SE"]
  
  # Create a density pot for the relevant variable for the relevant error values
  plot <- ggplot(data, aes(x = !!sym(variable), group = Error_Neg, fill = as.factor(Error_Neg))) +
    geom_density(alpha = 0.5, size = 0) +
    scale_fill_manual(
      values = brewer.pal(n = 5, "Set2"), 
      labels = function(x) sprintf("%.3f", as.numeric(x))     
    ) +
    geom_vline(xintercept = realvalue) +
    labs(y = "Density",
         fill = "False negative \n probability") +
    xlim(dimensions_x) +
    ylim(dimensions_y) +
    #scale_fill_discrete(labels = function(x) sprintf("%.3f", as.numeric(x))) +
    theme_minimal() +
    theme(axis.title.y = element_text(margin = margin(r = 15), size = 16),
          axis.title.x = element_text(margin = margin(t = 15), size = 16),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.position = "none") 
  print("plot created")  
  return(plot)
}

### Save created graphs
SaveGraph <- function(variable, error_pos, random_results, directory, error_neg, 
                      dimensions_x, dimensions_y){
  # Create graph and safe to given directory
  directory <- paste0(directory, variable, "_", error_pos, ".jpg")
  CreateGraph(datasummary = random_results$summary_theta, data = random_results$filtered_output$theta, 
              variable = c(variable), error_neg = error_neg, error_pos = error_pos,
              dimensions_x = dimensions_x, dimensions_y = dimensions_y)
  ggsave(directory, width = 2400, height = 1900, units = "px", dpi = 300)
  gc(verbose = FALSE)
}


