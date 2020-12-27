# import requirered packages
library(tidyverse)
library(matlib)
library(FactoMineR)
library(psych)
library(base)
library(sparsepca)
library(emdbook)
library(hrbrthemes)
library(reshape2)
library(ggplot2)
library(MASS)
library(gtools)
library(devtools)
library(tictoc)
library(stats)
library(dplyr)
library(parallel)

# import the source files
file_list <- c('lucas_make_data.R',
               'lucas_run_pca.R',
               'lucas_run_spca.R',
               'lucas_pc_select.R')
for(file in file_list){
  source(file)
}

# open the connection to the clusters
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
clusterExport(cl, 'spca')

# set up the simulation parameters
MC_reps <- 1
num_comp_seq <- c(3)
num_var_seq <- c(5)
variance_of_comps <- rep(0.5, length(num_comp_seq))
error_seq <- c(0.01, 0.05, 0.10)
n_seq <- c(25, 100, 250)
nrFolds <- 3

# define a function to run a single simulation
run_simulation <- function(sim_param_combs, variance_of_comps, nrFolds){
  
  # extract the simulation parameters from sim_param_combs
  num_mc <- sim_param_combs[1]
  p <- sim_param_combs[2]
  c <- sim_param_combs[3]
  e <- sim_param_combs[4]
  n <- sim_param_combs[5]
  
  # set seed for reproducibility 
  set.seed(num_mc)
  
  # generate the simulation data
  res_var <- makeVariance(varianceOfComps = variance_of_comps, p = p, error = e)
  X <- makeDat(n = n, p = p, ncomp = c, variances = res_var)$X
  
  # run the component selection techniques
  final_res <- RUN_PC_SELECT(dat = X, nrFolds = nrFolds)
  
  # extract the results for eigenvector cross-validation for both PCA and SPCA
  res_comp_pca <- final_res$PCA_CV
  res_comp_spca <- final_res$SPCA_CV
  
  # return a list containing the results
  return(c(num_mc = num_mc,
           p = p,
           error = e,
           n = n,
           ncomp = c,
           pca_ncomp = res_comp_pca,
           spca_ncomp = res_comp_spca))
}

# use 'tic' to check runtime
tic('run_simulation')

# create a list of all possible combinations of the simulation parameters
sim_param_combs <- expand.grid(seq(1, MC_reps), 
                               num_var_seq, 
                               num_comp_seq, 
                               error_seq, 
                               n_seq)

# run the simulations
results <- apply(X = sim_param_combs,
                 MARGIN = 1,
                 FUN = run_simulation,
                 variance_of_comps = variance_of_comps,
                 nrFolds = nrFolds)
toc()

# close the connections to the cluster
stopCluster(cl)

# convert the results to a data frame for readability
data.frame(MC = results[1,],
           p = results[2,],
           error = results[3,],
           n = results[4,],
           ncomp = results[5,],
           pca_ncomp = results[6,],
           spca_ncomp = results[7,])


# # next step is to create plots?
# # or to look at if the current CV method is correct
# # and to look at which other methods are useful? e.g. KG-criterion?
# # is parallel analysis really not an option?


