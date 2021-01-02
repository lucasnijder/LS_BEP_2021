

ParallelAnalysis <- function(dat, model_eigenval, type, n_datasets = 20){
  n_real <- nrow(dat)
  p_real <- ncol(dat)
  random_ev_matrix <- matrix(NA, nrow = n_datasets, ncol = min(dim(dat)))
  
  if(type == 'pca'){
    for(i in 1:n_datasets){
      random_dat <- matrix( rnorm(n_real*p_real,mean=mean(dat),sd=sd(dat)), n_real, p_real)
      random_ev <- prcomp(random_dat)$sdev^2
      random_ev_matrix[i,] <- random_ev
    }
  }
  
  if(type == 'spca'){
    for(i in 1:n_datasets){
      random_dat <- matrix( rnorm(n_real*p_real,mean=mean(dat),sd=sd(dat)), n_real, p_real)
      
      # the SPCA version without CV is used because of computational complexity
      random_ev <- spca(random_dat)$eigenvalues
      random_ev_matrix[i,] <- random_ev
    }
  }
  
  # change to 95% CI or some other measure
  random_ev_means <- colMeans(random_ev_matrix)
  
  # [1:length(model_eigenval)] is added because the random matrices have more eigenvalues and therefor the 
  # random_ev_means vector could have more values than the model_eigenval vector. Comparison to a non-existing
  # value leads to TRUE, meaning that the number of components is overestimated.
  res <- model_eigenval > random_ev_means[1:length(model_eigenval)]
  
  num_of_comp <- sum(res)
  return(num_of_comp)
}
