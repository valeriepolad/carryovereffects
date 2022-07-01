#### setup and load data and packages ####
library(MASS)
library(lme4)
library(psych)
library(MLmetrics)
library(GPArotation)
library(mirt)
library(mirtCAT)
library(shiny)
library(stats)

# set seed to check results
set.seed(1234)

# import reference data
setwd("~/Documents/research_project")
setwd(paste0(getwd(), "/romero_data"))
reference_data <- read.csv("romero_data.csv")

# split into 4 waves
time_1 <- reference_data[3:12]
time_2 <- reference_data[13:22]
time_3 <- reference_data[23:32]
time_4 <- reference_data[33:42]

# remove completely empty rows in time_1 (row 89 and 93), leaving 114 observations
time_1 <- time_1[-(c(89,93)),]

# calculate omega with 1 factor
omega(time_1,1) 
# omega for wave 1 = .81


#### scratch work for mirt::simdata ####
# In the following code snippet, data is being simulated using the model derived from wave 1 of the Romero dataset.

# fit mirt model to time_1 data with 1 factor
model <- mirt(time_1,1)
summary(model)

# simulate 100 theta values from normal distribution with mu = 0 and sigma = 1
simulated_theta <- matrix(rnorm(100))

# simulate 200 observations using model parameters from time_1
simulated_output <- simdata(model = model, N = 200)
# simulate 200 observations using model parameters from time_1 BUT with simulated theta values
simulated_output2 <- simdata(model = model, N = 200, Theta = simulated_theta)

omega(simulated_output,1)
# omega = .81, with time_1 theta
omega(simulated_output2,1)
# omega = .87, with simulated theta

#### scratch work for mirtCAT ####

# number of items
n_items <- 20
# population mean
mu <- 0
# population sd
sigma <- 1

# develop slope parameters
a1 <- rlnorm(n_items,0,1)

# simulate theta values from normal distribution / matrix of intercepts
d <- rnorm(n_items,mu,1)

# create data frame of parameters to feed into generate.mirt_object()
parameters <- data.frame(a1=a1, d=d)

# generate mirt model
obj <- generate.mirt_object(parameters = parameters, itemtype = '2PL')

# simulate data from mirt model
dat <- simdata(N=200, model=obj)

# plot item trace lines
plot(obj, type = 'trace')

# calculate omega for simulated data with 1 factor
omega(dat,1)
# omega = .82

# graded response model
# look up for dichotomous outcomes

# slope and threshold

# IVs
# keep to 2 levels 
# average loading 0.1 to 0.7 from factor models to IRT

# population mean
mu <- 0
# population standard deviation
sigma <- 1
# latent trait
theta <- matrix(rnorm(100))
# sample size / note: when you have the theta values, sample size is determined by length(theta)
n <- seq(100,1000,100)
# number of items in scale
number_of_items <- seq(5,15,5)
# number of splits of scale (create if/then statement based on number of items)
number_of_splits <- 1:5
# number of items in each split
## number_of_items_per_split <- number_of_items/number_of_splits
# carryover (TBD)
carry_over <- NA
# method: split half implementation OR traditional implementation (no split half)


# DVs

# type I error
# statistical power
# bias (observed - true)
# mean squared error
# reliability (alpha and omega)


# use simdata function in mirt, example below
# https://www.rdocumentation.org/packages/mirt/versions/1.32.1/topics/simdata

# generate.mirt_object(
# parameters,
# itemtype,
# latent_means = NULL,
# latent_covariance = NULL,
# key = NULL,
# min_category = rep(0L, length(itemtype))
# )

# https://psychometroscar.com/simulate-cronbachs-alpha/ <--- VERY helpful


# Since reliability is a measure of the covariance matrix of the items, I decided to work backwards from that.
# Here, we simulate data from a multivariate normal distribution with different covariance matrices.

## Where we're at right now ##
# Can generate data with different number of items and differing levels of reliability

# In general, 
# omega = .3 when value = .1
# omega = .5 when value = .3
# omega = .8 when value = .5
# omega = .95 when value = .8


# simulation

# use switch statement for omega -> value
n_items <- 4
sim_replications <- 100
dim_names = list(c("Sigma = .1", "Sigma = .3", "Sigma = .5", "Sigma = .8"))
all_results <- matrix(rep(0,n_items*n_items),n_items,n_items, dimnames = dim_names)
for (i in 1:n_items){
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
    
    results <- vector('numeric', sim_replications)
    for (k in 1:sim_replications){
      datam <- as.data.frame(mvrnorm(sim_replications, c(rep(0,n_items)), sigma))
      results[k] <- as.numeric(as.character(omega(datam,1)))
    }
    all_results[j,i] <- mean(results)
  }
}


# average size of loading
# variance of loading
# should determine reliability
# covariance matrices where reliability is the same but the variance is different

datam <- as.data.frame(mvrnorm(100, c(rep(0,4)), sigma))

omega(datam)




