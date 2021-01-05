run_component <- function(c, res_loadings, train_dat, test_dat){
  # W represents the loadings of the (S)PCA model
  W <- matrix(res_loadings[,1:c], nrow = ncol(train_dat), ncol = c)
  
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
  
  PRESS <- sum((test_dat - pred)^2)
  
  return(PRESS)
}

run_fold <- function(i, dat, folds, type, L1=0.001, L2=0.01){
  fold <- which(folds == i)
  
  train_dat <- dat[-fold,]
  test_dat <- dat[fold,]
  
  if(type == 'pca'){
    res_model <- prcomp(train_dat)
    res_loadings <- res_model$rotation
  }
  if(type == 'spca'){
    res_model <- spca(train_dat, alpha = L1, beta = L2) # no need to restrict the number of components?
    res_loadings <- res_model$loadings
  }
  
  components_vector <- 1:ncol(res_loadings)
  PRESS_list <- lapply(X = components_vector, 
                            FUN = run_component, 
                            res_loadings = res_loadings,
                            train_dat = train_dat,
                            test_dat = test_dat)
  
  return(PRESS_list)
}

CV_method <- function(dat, type, L1=NULL, L2=NULL, nrFolds=10){
  folds <- rep_len(1:nrFolds, nrow(dat))
  
  folds_vector <- 1:nrFolds
  PRESS_matrix <- lapply(X = folds_vector, 
                              FUN = run_fold, 
                              dat = dat, 
                              folds = folds,
                              type = type, L1 = L1, L2 = L2)
  
  PRESS_matrix <- t(matrix(unlist(PRESS_matrix), nrow=nrFolds))
  
  col_sums_PRESS <- colSums(PRESS_matrix)
  
  num_of_comp <- which.min(col_sums_PRESS)
  
  return(num_of_comp)
}