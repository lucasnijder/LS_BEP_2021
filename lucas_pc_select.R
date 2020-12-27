EigenVectorCV <- function(dat, type, nrFolds=10){
  
  cvError <- matrix(NA, nrow=nrFolds, ncol=min(dim(dat))-1)
  folds <- rep_len(1:nrFolds, nrow(dat))

  for (i in 1:nrFolds) { 
    fold <- which(folds == i)
    
    trainDat <- dat[-fold, ]
    testDat <- dat[fold, , drop = FALSE]
    
    if(type == 'pca'){
      res_loadings <- RUN_PCA(trainDat)    
      best_l1l2 <- list('NO','NO')
    }
    if(type == 'spca'){
      
      res <- RUN_SPCA_CV(trainDat)
      res_loadings <- unlist(res$spca_loadings)
      best_l1l2 <- list(res$L1, res$L2)
    }
    
    for (c in 1:ncol(res_loadings)){
      
      # W represents the loadings of the PCA model
      W <- matrix(res_loadings[,1:c], nrow = ncol(trainDat), ncol = c)

      # apply SVD to get P
      Xt_X_W <- t(trainDat) %*% trainDat %*% W
      res_svd <- svd(Xt_X_W)
      P <- res_svd$u %*% t(res_svd$v)
      
      #Eigenvector crossvalidation Bro Kiers
      pred <- matrix(NA, nrow(testDat), ncol(testDat))
      
      for(j in 1:ncol(trainDat)){
        TMinusVariableJ <- testDat[, -j] %*% W[-j, ] 
        pred[, j] <- TMinusVariableJ %*% P[j, ]
      }
      
      comp_error <- sum((testDat - pred)^2)
      cvError[i,c] <- comp_error
    }
  }
  
  CV_res <- matrix(colSums(cvError), nrow=1, ncol=ncol(cvError))
  
  c_name_l <- c()
  for(a in 1:ncol(CV_res)){
    c_name_l <- append(c_name_l, paste('pc',a,sep=''))
  }
  
  colnames(CV_res) <- c_name_l
  GlobalMin <- which.min(CV_res)
  
  return(list(CV_res, GlobalMin, best_l1l2))
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
