## R code for analyzing the example data using the proposed MMoM method
## The proposed method is implemented in the ‘mmeta’ function in our 
## R software package ‘xmeta’, and can be installed from CRAN 
## (http://cran.r-project.org/package=xmeta/), the official R 
## package archive.
## Proposed method : Multivariate Method of Moments with working 
## independence assumption
## Input:
## 1. ys: effect sizes (y1 and y2) from n studies in a matrix format 
## with two columns: y1 and y2
## 2. vars: an (nx2) matrix with column 1 being variance of y1 and 
## column 2 being variance of y2
## Note: when y1 or y2 is missing, we put y1 or y2 as 0 and 
## corresponding variances as 10^6.
## Output: beta.hat (the pooled effect sizes for two outcomes) and 
## sigma.hat (covariance matrix of the estimates)
# Start time
# Load necessary libraries
library(readr)

MMoM <- function(ys, vars){
  y1 = ys[,1]
  y2 = ys[,2]
  var1 = vars[,1]
  var2 = vars[,2]
  w1 = 1/(var1)
  w2 = 1/(var2)
  
  y1.weight = sum(w1*y1)/sum(w1)
  y2.weight = sum(w2*y2)/sum(w2)
  n1 = sum(1*(var1 > 10^-4))
  n2 = sum(1*(var2 > 10^-4))
  Q1 = sum((1-w1/sum(w1))^2)
  Q2 = sum((1-w2/sum(w2))^2)
  tau1.2.hat = max(0, (Q1-(n1-1))/((sum(w1)-sum(w1^2)/sum(w1))))
  tau2.2.hat = max(0, (Q2-(n2-1))/((sum(w2)-sum(w2^2)/sum(w2))))
  
  ## variance estimate:
  w1.star = 1/(var1 + tau1.2.hat)
  w2.star = 1/(var2 + tau2.2.hat)
  
  print(sum(w1.star))
  
  
  beta1.hat = sum(w1.star*y1)/sum(w1.star)
  beta2.hat = sum(w2.star*y2)/sum(w2.star)
  var.beta1.hat = 1/sum(w1.star)
  var.beta2.hat = 1/sum(w2.star)
  
  
  mycov.beta = sum((w1.star/sum(w1.star))*(w2.star/sum(w2.star)) * (y1 - beta1.hat)*(y2 - beta2.hat))
  mycov.beta = sum(( var.beta1.hat* var.beta2.hat) * (y1 - beta1.hat)*(y2 - beta2.hat))
  
  beta.hat = c(beta1.hat, beta2.hat)
  sigma.hat = matrix(c(var.beta1.hat, mycov.beta, mycov.beta, 
                       var.beta2.hat), nrow = 2, byrow = T)
  result = list(beta.hat = beta.hat, beta.cov = sigma.hat)
  return(result)
}
 


# Load the data
df <- read_delim("fastmeta-main/fastmeta-main/generated_data.tab.txt", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
 

# Convert standard errors to variances
df$var1 <- df$SE1^2
df$var2 <- df$SE2^2

# Get unique variable names
unique_vars <- unique(df$variable)
start_time <- Sys.time()

# Iterate over each unique variable name
for (var_name in unique_vars) {
  
  # Subset data for the current variable
  subset_df <- subset(df, variable == var_name)
  
  # Construct the ys matrix (effects: beta1, beta2)
  ys <- matrix(c(subset_df$BETA1, subset_df$beta2), ncol = 2, byrow = FALSE)
  
  # Construct the vars matrix (variances: var1, var2)
  vars <- matrix(c(subset_df$BETA2, subset_df$var2), ncol = 2, byrow = FALSE)
  
 
  MMoM.fit <- MMoM(ys, vars)
  
 

}
 
 

end_time <- Sys.time()

execution_time <- end_time - start_time
print(execution_time)