#### setup ####
library(MASS)
library(lme4)
library(psych)
library(MLmetrics)
library(GPArotation)
library(mirt)
library(shiny)
library(stats)
library(mice)

# set seed to check results
set.seed(1234)

# https://psychometroscar.com/simulate-cronbachs-alpha/ <--- VERY helpful


# Since reliability is a measure of the covariance matrix of the items, I decided to work backwards from that.
# Here, we simulate data from a multivariate normal distribution with different covariance matrices.

## Where we're at right now ##
# Can generate data with different number of items and differing levels of reliability

# In general, 
# omega = .3 when sigma = .1
# omega = .5 when sigma = .3
# omega = .8 when sigma = .5
# omega = .95 when sigma = .8

#### Begin simulation ####

# number of items / can also be a vector to be fed into for loop
n_items <- 8
# number of replications per cell
sim_replications <- 10
# housekeeping / row and column names
sigma_values <- c("Sigma = .1", "Sigma = .3", "Sigma = .5", "Sigma = .8")
item_names <- paste(rep("Item", n_items), 1:n_items, sep=" ")
dim_names = list(sigma_values, item_names)
# matrix to store all of the results
all_results <- matrix(rep(0,n_items*n_items),4,n_items,dimnames = dim_names)

# for each item
for (i in 1:n_items){
  # for each case of sigma
  for (j in 1:4){
    # if j = case, set omega to that value
    omega <- switch(j,
                    "0.3", # case 1
                    "0.5", # case 2
                    "0.8", # case 3
                    "0.95", # case 4
                    print("Invalid entry."))
    # if omega = certain value, construct covariance matrix that corresponds to that omega value
    matrix_value <- switch(omega,
                           "0.3"={.1}, 
                           "0.5"={.3}, 
                           "0.8"={.5}, 
                           "0.95"={.8},
                           print("Valid omega value entries are .3, .5, .8, and .95."))
    sigma <- matrix(rep(matrix_value,n_items*n_items),n_items,n_items)
    diag(sigma) <- 1
    
    # vectors for storing results
    results_wave1 <- vector('numeric', sim_replications)
    results_wave2 <- vector('numeric', sim_replications)
    omega_diff <- vector('numeric', sim_replications)
    
    # for each replication, construct data matrix and then calculate and store omega
    for (k in 1:sim_replications){
      data_matrix <- as.data.frame(mvrnorm(sim_replications, c(rep(0,n_items)), sigma))
      
      # create permutatioon matrix for subsetting
      x <- 1:n_items
      permutation_matrix <- do.call(rbind, rep(list(x), sim_replications))
      for (m in 1:sim_replications){
        permutation_matrix[m,] <- sample(n_items)
      }
      
      # subsetting to construct split halves
      permutation_matrix_wave1 <- replace(permutation_matrix, permutation_matrix > (n_items/2), TRUE)
      permutation_matrix_wave1 <- replace(permutation_matrix_wave1, permutation_matrix <= (n_items/2), FALSE)
      permutation_matrix_wave2 <- replace(permutation_matrix, permutation_matrix <= (n_items/2), TRUE)
      permutation_matrix_wave2 <- replace(permutation_matrix_wave2, permutation_matrix > (n_items/2), FALSE)
      
      
      # multiply by 0s and 1s to remove certain values
      wave1 <- data_matrix*permutation_matrix_wave1
      wave2 <- data_matrix*permutation_matrix_wave2
      
      # set removed values to NA
      wave1 <- replace(wave1,wave1 == 0, NA)
      wave2 <- replace(wave2,wave2 == 0, NA)
      
      # impute missing data for omega function
      # possible issue here regarding (non)removal of collinear items for imputation
      imputed_data <- mice(wave1, m=5, maxit = 10, method = 'pmm', seed = 500,remove.collinear=FALSE)
      wave1 <- complete(imputed_data,2)
      imputed_data <- mice(wave2, m=5, maxit = 10, method = 'pmm', seed = 500,remove.collinear=FALSE)
      wave2 <- complete(imputed_data,2)
      
      # calculate reliability for random split halves
      omega_wave1 <- as.numeric(as.character(omega(wave1,1)))[1]
      omega_wave2 <- as.numeric(as.character(omega(wave2,1)))[1]
      
      results_wave1[k] <- omega_wave1
      results_wave2[k] <- omega_wave2
      
      # calculate and store difference in omega between waves
      omega_diff[k] <- results_wave2[k]-results_wave1[k]
    }
    
    # average across all replications for that specific condition and store in matrix
    all_results[j,i] <- mean(omega_diff) #this may be where the error is
  }
}
