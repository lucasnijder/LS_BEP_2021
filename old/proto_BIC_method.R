
run_component <- function(c, dat, res_loadings, residual_variance0){
  # W represents the loadings of the (S)PCA model
  W <- matrix(res_loadings[,1:c], nrow = ncol(dat), ncol = c)
  T_hat <- dat %*% W
  
  # apply SVD to get P
  Xt_X_W <- t(dat) %*% dat %*% W
  res_svd <- svd(Xt_X_W)
  P <- res_svd$u %*% t(res_svd$v)
  
  residual_variance <- sum((dat - T_hat %*% t(P))^2)
  
  df <- ncol(res_loadings) - 1
  
  I <- nrow(dat)
  
  print(residual_variance)
  print(residual_variance0)
  
  BIC <- (residual_variance / residual_variance0) + (df * (log(I) / I))
  
  return(BIC)
}

BIC_method <- function(dat, type, L1=0.001, L2=0.01){
  
  if(type == 'pca'){
    res_model <- prcomp(dat)
    res_loadings <- res_model$rotation
  }
  if(type == 'spca'){
    res_model <- spca(dat, alpha = L1, beta = L2)
    res_loadings <- res_model$loadings
  }
  
  W0 <- matrix(res_loadings, nrow = ncol(dat))
  T_hat0 <- dat %*% W0
  
  # apply SVD to get P
  Xt_X_W0 <- t(dat) %*% dat %*% W0
  res_svd0 <- svd(Xt_X_W0)
  P0 <- res_svd0$u %*% t(res_svd0$v)
  
  residual_variance0 <- sum((dat - T_hat0 %*% t(P0))^2)/length(dat)
  
  components_vector <- 1:ncol(res_loadings)
  BIC_list <- lapply(X = components_vector, 
                       FUN = run_component, 
                       res_loadings = res_loadings,
                       residual_variance0 = residual_variance0,
                       dat = dat)
  
  print(BIC_list)
  
  num_of_comps <- which.min(unlist(BIC_list))
  
  return(num_of_comps)
}