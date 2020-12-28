# extracts the corresponding results. Was created because in the results of EV_run_folds
# --> the pca and spca results are in nested lists
select_model_result <- function(x, type){
  if(type == 'pca'){
    return(x$pca_comp_errors)
  }
  if(type == 'spca'){
    return(x$spca_comp_errors)
  }
}

# a function that extracts the number of components for the specified model type. Was created because 
# --> in the results of EV_run_folds the pca and spca results are in nested lists
extract_num_comps <- function(folds_comp_errors, type, dat){
  folds_comp_errors <- lapply(X = folds_comp_errors, FUN = select_model_result, type = type)
  
  folds_comp_error_matrix <- matrix(unlist(folds_comp_errors), ncol = min(dim(dat))-1, nrow = nrFolds)
  
  PRESS_per_number_of_components <- colSums(folds_comp_error_matrix)
  
  global_min = which.min(PRESS_per_number_of_components)

  return(global_min)
}

# is called for each component and returns the PRESS for that component
EV_PRESS_per_component <- function(c, res_loadings, train_dat, test_dat){
  
  # W represents the loadings of the (S)PCA model
  W <- matrix(res_loadings[,1:c], nrow = ncol(train_dat), ncol = c) # removed 1:
  
  # apply SVD to get P
  Xt_X_W <- t(train_dat) %*% train_dat %*% W
  res_svd <- svd(Xt_X_W)
  P <- res_svd$u %*% t(res_svd$v)
  
  # predict each value x based on the test set minus x as suggested by Bro et al. (2008)
  pred <- matrix(NA, nrow(test_dat), ncol(test_dat))
  
  for(j in 1:ncol(test_dat)){
    TMinusVariableJ <- test_dat[, -j] %*% W[-j, ] # ----------- look at this again! ---------
    pred[, j] <- TMinusVariableJ %*% P[j, ]
  }
  
  # calculate the PRESS for each number of components
  comp_error <- sum((test_dat - pred)^2)
  
  return(comp_error)
}

# is called for each fold and returns a list of PRESS scores for each number of components (for both PCA and SPCA)
EV_run_folds <- function(i, dat, folds){
  fold <- which(folds == i)
  
  # split up the data into train and test sets
  train_dat <- dat[-fold, ]
  test_dat <- dat[fold,]
  
  res_pca <- RUN_PCA(train_dat)
  pca_loadings <- unlist(res_pca$rotation)
  
  # run cross-validation to get the optimal L1 and L2 parameters and return best performing SPCA model
  res_spca <- RUN_SPCA_CV(train_dat)
  spca_loadings <- unlist(res_spca$spca_loadings)
  
  pca_components_vector <- 1:ncol(pca_loadings)
  pca_comp_errors <- lapply(X = pca_components_vector, 
                        FUN = EV_PRESS_per_component, 
                        res_loadings = pca_loadings,
                        train_dat = train_dat,
                        test_dat = test_dat)
  
  spca_components_vector <- 1:ncol(spca_loadings)
  spca_comp_errors <- lapply(X = spca_components_vector, 
                                  FUN = EV_PRESS_per_component, 
                                  res_loadings = spca_loadings,
                                  train_dat = train_dat,
                                  test_dat = test_dat)
  
  return(list(pca_comp_errors = pca_comp_errors, spca_comp_errors = spca_comp_errors))
}

# create a function to run eigenvector cross-validation 
EigenVectorCV <- function(dat, nrFolds=10){
  
  folds <- rep_len(1:nrFolds, nrow(dat))
  
  folds_vector <- 1:nrFolds
  folds_comp_errors <- lapply(X = folds_vector, 
                              FUN = EV_run_folds, 
                              dat = dat, 
                              folds = folds)
  
  pca_global_min = extract_num_comps(folds_comp_errors, 'pca', dat)
  spca_global_min = extract_num_comps(folds_comp_errors, 'spca', dat)

  return(list(pca_global_min = pca_global_min, spca_global_min = spca_global_min))
}

RUN_PC_SELECT <- function(dat, nrFolds){
  
  tic('PC cross-validation')
  res_CV <- EigenVectorCV(dat = dat, nrFolds = nrFolds)
  toc()
  
  return(list(PCA_CV = res_CV[1], SPCA_CV = res_CV[2]))
}
