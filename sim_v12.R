# setup

library(MASS)
library(lme4)
library(psych)
library(MLmetrics)
library(GPArotation)
library(mirt)
library(shiny)
library(stats)
library(mice)

set.seed(1)

# sample size
n_people <- 100

# number of items
n_items_vec <- c(4,8)

# carryover
carryover_value <- .15

# number of replications per cell
replications_per_cell <- 10

# setup table dimensions + create table to store results
## add back in .8 condition
sigma_values <- c(.1,.3,.5,.8)
col_names <- paste(n_items_vec,rep("Items", length(n_items_vec)), sep=" ")
dim_names = list(sigma_values, col_names)

all_results <- matrix(rep(0,length(n_items_vec)*length(n_items_vec)),length(sigma_values),length(n_items_vec), dimnames = dim_names)


# begin sim
# method: 0 = no split-halves, 1 = split-halves
for (method in 0:1){
# if traditional method
  if (method == 0){
  # for every value of n_items_vec
  for (i in 1:length(n_items_vec)){
    n_items <- n_items_vec[i]
    
    # for every value of sigma_values
    for (j in 1:length(sigma_values)){
      sigma <- sigma_values[j]
      
      # create covariance matrix for latent trait
      covariance_matrix <- matrix(rep(sigma,n_items*n_items),n_items,n_items)
      diag(covariance_matrix) <- 1
      
      # create covariance matrix with carryover effects (if present)
      
      # create mini diagonal chunks
      diagonal_sigma <- sigma/2 + carryover_value
      carryover_matrix_chunk <- matrix(rep((sigma/2),n_items*n_items),n_items,n_items)
      diag(carryover_matrix_chunk) <- diagonal_sigma
      
      # combine all 4 pieces
      top_half <- cbind(covariance_matrix,carryover_matrix_chunk)
      bottom_half <- cbind(carryover_matrix_chunk,covariance_matrix)
      full_covariance_matrix <- rbind(top_half,bottom_half)
      
      # vectors for storing results
      results_wave1 <- vector('numeric', replications_per_cell)
      results_wave2 <- vector('numeric', replications_per_cell)
      omega_diff <- vector('numeric', replications_per_cell)
      
      # for each replication, construct data matrix and then calculate and store omega
      for (k in 1:replications_per_cell){
        # generate data once
        data_matrix <- as.data.frame(mvrnorm(n_people, c(rep(0,n_items*2)), full_covariance_matrix))
        
        # split data into wave 1 and wave 2
        wave1 <- data_matrix[,1:n_items]
        wave2 <- data_matrix[,(n_items+1):(ncol(data_matrix))]
        
        # calculate reliability for random split halves
        omega_wave1 <- as.numeric(as.character(omega(wave1,1)))[1]
        omega_wave2 <- as.numeric(as.character(omega(wave2,1)))[1]
        
        results_wave1[k] <- omega_wave1
        results_wave2[k] <- omega_wave2
        
        # calculate and store difference in omega between waves
        omega_diff[k] <- results_wave2[k]-results_wave1[k]
        }
      all_results[j,i] <- mean(omega_diff)
      }
  }
    all_results_traditional <- all_results
}
  if (method ==1){
  # for every value of n_items_vec
  for (i in 1:length(n_items_vec)){
    n_items <- n_items_vec[i]
    
    # for every value of sigma_values
    for (j in 1:length(sigma_values)){
      sigma <- sigma_values[j]
      
      # create covariance matrix for latent trait
      covariance_matrix <- matrix(rep(sigma,n_items*n_items),n_items,n_items)
      diag(covariance_matrix) <- 1
      
      # create covariance matrix with carryover effects (if present)
      
      # create mini diagonal chunks
      diagonal_sigma <- sigma/2 + carryover_value
      # introduce gaussian white noise
      carryover_matrix_chunk <- matrix(rep((sigma/2),n_items*n_items),n_items,n_items)
      diag(carryover_matrix_chunk) <- diagonal_sigma
      
      # combine all 4 pieces
      top_half <- cbind(covariance_matrix,carryover_matrix_chunk)
      bottom_half <- cbind(carryover_matrix_chunk,covariance_matrix)
      full_covariance_matrix <- rbind(top_half,bottom_half)
      
      # vectors for storing results
      results_wave1 <- vector('numeric', replications_per_cell)
      results_wave2 <- vector('numeric', replications_per_cell)
      omega_diff <- vector('numeric', replications_per_cell)
      
      # for each replication, construct data matrix and then calculate and store omega
      for (k in 1:replications_per_cell){
        # generate data once
        data_matrix <- as.data.frame(mvrnorm(n_people, c(rep(0,n_items*2)), full_covariance_matrix))
        
        # create permutation matrix for subsetting
        x <- 1:(n_items*2)
        permutation_matrix <- do.call(rbind, rep(list(x), n_people))
        for (m in 1:n_people){
          permutation_matrix[m,] <- sample(n_items*2)
        }
        
        # subsetting to construct split halves
        permutation_matrix_wave1 <- replace(permutation_matrix, permutation_matrix > (n_items), TRUE)
        permutation_matrix_wave1 <- replace(permutation_matrix_wave1, permutation_matrix <= (n_items), FALSE)
        permutation_matrix_wave2 <- replace(permutation_matrix, permutation_matrix <= (n_items), TRUE)
        permutation_matrix_wave2 <- replace(permutation_matrix_wave2, permutation_matrix > (n_items), FALSE)
        
        
        # multiply by 0s and 1s to remove certain values
        wave1 <- data_matrix*permutation_matrix_wave1
        wave2 <- data_matrix*permutation_matrix_wave2
        
        # set removed values to NA
        # AKM: try to do NAs before multiplying
        wave1 <- replace(wave1,wave1 == 0, NA)
        wave2 <- replace(wave2,wave2 == 0, NA)
        
        # impute missing data for omega function
        # possible issue here regarding (non)removal of collinear items for imputation
        # AKM: look at doing FIML, or combining across imputed datasets
        # update for above comment: changed parameter for complete to "long" to pool across inputed datasets
        # https://www.rdocumentation.org/packages/mice/versions/2.46.0/topics/complete
        # change m to 20 once code can be run on supercomputer
        imputed_data <- mice(wave1, m=5, maxit = 10, method = 'pmm', seed = 500,remove.collinear=FALSE)
        wave1 <- complete(imputed_data,"long")
        imputed_data <- mice(wave2, m=5, maxit = 10, method = 'pmm', seed = 500,remove.collinear=FALSE)
        wave2 <- complete(imputed_data,"long")
        
        # remove marker variables
        wave1 <- wave1[,-2:-1]
        wave2 <- wave2[,-2:-1]
        
        # calculate reliability for halves
        omega_wave1 <- as.numeric(as.character(omega(wave1,1)))[1]
        omega_wave2 <- as.numeric(as.character(omega(wave2,1)))[1]
        
        results_wave1[k] <- omega_wave1
        results_wave2[k] <- omega_wave2
        
        # calculate and store difference in omega between waves
        omega_diff[k] <- results_wave2[k]-results_wave1[k]
        }
      all_results[j,i] <- mean(omega_diff)
    }
  }
    all_results_split_halves <- all_results
  }
}

# after running both traditional (method == 0) and split-halves (method == 1), compare results
method_diff <- all_results_traditional-all_results_split_halves