lseq <- function(from=0.001, to=.01, length.out=3) {
  exp(seq(log(from), log(to), length.out = length.out))
}

spca_cv_1 <- function(param_combs, dat, nfolds, ridge, lasso){
  
  folds <- rep_len(1:nfolds, nrow(dat))
  fold <- which(folds == param_combs[1])
  
  ridge <- param_combs[2]
  lasso <- param_combs[3]
  
  trainDat <- dat[-fold, ]
  testDat <- dat[fold, , drop = FALSE]
  
  res_spca <- spca(X = dat, k = min(dim(dat))-1, alpha = ridge, beta = lasso, max_iter = 10000, verbose = FALSE)
  
  W <- matrix(res_spca$loadings, nrow = ncol(dat), ncol = min(dim(dat))-1)
  
  # apply SVD to get P
  Xt_X_W <- t(trainDat) %*% trainDat %*% W
  res_svd <- svd(Xt_X_W)
  P <- res_svd$u %*% t(res_svd$v)
  
  pred <- matrix(NA, nrow(testDat), ncol(testDat))
  for(i in 1:ncol(dat)){
    TMinusVariableJ <- testDat[, -i] %*% W[-i, ]
    pred[, i] <- TMinusVariableJ %*% P[i, ] 
  }
  
  comb_MSE <- mean((testDat - pred)^2)
  
  return(c(param_combs[1], param_combs[2], param_combs[3], comb_MSE))
}

RUN_SPCA_CV <- function(dat, nfolds=3, ridge_seq=lseq(), lasso_seq=lseq()){

  param_combs <- expand.grid(seq(1, nfolds), ridge_seq, lasso_seq)
  
  res_spca_cv <- parApply(cl = cl, X = param_combs, MARGIN = 1, FUN = spca_cv_1, dat = dat, nfolds=nfolds)
  
  df_res_spca_cv <- data.frame(fold=res_spca_cv[1,],ridge=res_spca_cv[2,],lasso=res_spca_cv[3,],comb_MSE=res_spca_cv[4,])
  
  df_res_spca_cv <- df_res_spca_cv %>%
    group_by(ridge, lasso) %>%
    summarise(comb_MSE = mean(comb_MSE))
  
  row_min <- df_res_spca_cv[which.min(df_res_spca_cv$comb_MSE),]

  best_spca_res <- spca(dat, k =  min(dim(dat))-1, alpha = row_min$ridge, beta = row_min$lasso)
  best_spca_loadings <- best_spca_res$loadings

  return(list(MSE_matrix = df_res_spca_cv, L1 = row_min$ridge, L2 = row_min$lasso, mse = row_min$comb_MSE, spca_loadings = best_spca_loadings))
}
