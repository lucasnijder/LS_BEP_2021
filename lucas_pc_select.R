# create a function which checks which number of components has the lowest PRESS
EV_PRESS_per_component <- function(c, res_loadings, train_dat, test_dat){
  
  # W represents the loadings of the (S)PCA model
  W <- matrix(res_loadings[,1:c], nrow = ncol(train_dat), ncol = c)
  
  # apply SVD to get P
  Xt_X_W <- t(train_dat) %*% train_dat %*% W
  res_svd <- svd(Xt_X_W)
  P <- res_svd$u %*% t(res_svd$v)
  
  # predict each value x based on the test set minus x as suggested by Bro et al. (2008)
  pred <- matrix(NA, nrow(test_dat), ncol(test_dat))
  
  for(j in 1:ncol(train_dat)){
    TMinusVariableJ <- test_dat[, -j] %*% W[-j, ] # ----------- look at this again! ---------
    pred[, j] <- TMinusVariableJ %*% P[j, ]
  }
  
  # calculate the PRESS for each number of components
  comp_error <- sum((test_dat - pred)^2)
  
  return(comp_error)
}

# create a function to run all folds
EV_run_folds <- function(i, dat, folds, type){
  fold <- which(folds == i)
  
  # split up the data into train and test sets
  train_dat <- dat[-fold, ]
  test_dat <- dat[fold,]
  
  if(type == 'pca'){
    res <- RUN_PCA(train_dat)
    # extract PCA loadings
    res_loadings <- unlist(res$rotation)    }
  
  if(type == 'spca'){
    # run cross-validation to get the optimal L1 and L2 parameters and return best performing SPCA model
    res <- RUN_SPCA_CV(train_dat)
    # extract SPCA loadings
    res_loadings <- unlist(res$spca_loadings)
  }
  
  components_vector <- 1:ncol(res_loadings)
  comp_errors <- lapply(X = components_vector, 
                                  FUN = EV_PRESS_per_component, 
                                  res_loadings = res_loadings,
                                  train_dat = train_dat,
                                  test_dat = test_dat)
  return(comp_errors)
}

# create a function to run eigenvector cross-validation 
EigenVectorCV <- function(dat, type, nrFolds=10){
  
  # create a matrix to store the PRESS of a component per fold
  cvError <- matrix(NA, nrow=nrFolds, ncol=min(dim(dat))-1)
  folds <- rep_len(1:nrFolds, nrow(dat))
  
  folds_vector <- 1:nrFolds
  folds_comp_errors <- lapply(X = folds_vector, 
                              FUN = EV_run_folds, 
                              dat = dat, 
                              folds = folds, 
                              type = type)
  
  PRESS_per_number_of_component <- matrix(colSums(cvError), nrow=1, ncol=ncol(cvError))
  
  GlobalMin <- which.min(PRESS_per_number_of_component)
  
  return(GlobalMin)
}

RUN_PC_SELECT <- function(dat, nrFolds){
  
  # Cross-validation
  tic('PCA PC cross-validation')
  pca_res_CV <- EigenVectorCV(dat = dat, type = 'pca', nrFolds = nrFolds)
  toc()
  tic('SPCA PC cross-validation')
  spca_res_CV <- EigenVectorCV(dat = dat, type = 'spca', nrFolds = nrFolds)
  toc()
  
  return(list(PCA_CV = pca_res_CV, SPCA_CV = spca_res_CV))
}
