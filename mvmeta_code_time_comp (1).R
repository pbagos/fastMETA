library(mvmeta)

# 
# ### REPRODUCE THE RESULTS IN BERKEY ET AL. (1998)
# # INSPECT THE DATA
# berkey98
# # Method of Momments 
# year <- berkey98$pubyear - 1983
# model <- mvmeta(cbind(PD,AL)~year ,S=berkey98[5:7],data=berkey98,method="mm")
# print(summary(model),digits=3)
# 
# 
# # Load the data
# df <- read_delim("fastmeta-main/fastmeta-main/generated_data.tab.txt", 
#                  delim = "\t", escape_double = FALSE, 
#                  trim_ws = TRUE)
# 
# 
# 
# 
# # Convert standard errors to variances
# df$var1 <- df$SE1^2
# df$var2 <- df$SE2^2
# 
# 
# 
# df$covariance <-  0
# 
# 
# # Get unique variable names
# 
# unique_vars <- unique(df$variable)
# start_time <- Sys.time()
# 
# 
# 
# 
# # Iterate over each unique variable name
# for (var_name in unique_vars) {
#   
#   # Subset data for the current variable
#   subset_df <- subset(df, variable == var_name)
#  
#   model <- mvmeta(cbind(BETA1,BETA2) ,S=cbind(subset_df$var1,  subset_df$covariate,subset_df$var2),data=subset_df,method="mm")
#    
#   
# }
# 
# 
# end_time <- Sys.time()
# 
# execution_time <- end_time - start_time
# print(execution_time)


 
# List of file names
files <- list.files(path = "x_meta_data/", pattern = "*.tab", full.names = TRUE)

# Data frame to store results
results <- data.frame(File = character(), ExecutionTime = numeric(), stringsAsFactors = FALSE)

# Iterate over each file
for (file in files) {
  print(paste("Processing file:", file))
  

  
  # Load the data
  df <- read_delim("x_meta_data/", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # Convert standard errors to variances
  df$var1 <- df$SE1^2
  df$var2 <- df$SE2^2
 
  df$covariance <-  0
  
  # Get unique variable names
  unique_vars <- unique(df$variable)
  start_time <- Sys.time()
  # Iterate over each unique variable name
  for (var_name in unique_vars) {
    subset_df <- subset(df, variable == var_name)
 
    # Subset data for the current variable
    subset_df <- subset(df, variable == var_name)
    
    model <- mvmeta(cbind(BETA1,BETA2) ,S=cbind(subset_df$var1,  subset_df$covariate,subset_df$var2),data=subset_df,method="mm")
    
  }
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print(execution_time)
  # Append the result to the data frame
  results <- rbind(results, data.frame(File = basename(file), ExecutionTime = execution_time))
}

# Save the results to a file
write.csv(results, "mvmeta_execution_times.csv", row.names = FALSE)

print("Execution times have been saved to 'mvmeta_execution_times.csv'")
