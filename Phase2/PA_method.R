

PA_method <- function(dat, model_eigenval, type, n_datasets = 20){
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
  
  random_ev_CI <- apply(random_ev_matrix, 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
  
  # [1:length(model_eigenval)] is added because the random matrices have more eigenvalues and therefor the 
  # random_ev_means vector could have more values than the model_eigenval vector. Comparison to a non-existing
  # value leads to TRUE, meaning that the number of components is overestimated.
  res <- model_eigenval > random_ev_CI[1:length(model_eigenval)]
  
  num_of_comp <- sum(res)
  return(num_of_comp)
}
