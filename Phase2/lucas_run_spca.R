lseq <- function(from=0.00001, to=0.01, length.out=4) {
  exp(seq(log(from), log(to), length.out = length.out))
}

spca_cv_1 <- function(param_combs, dat, nfolds){
  
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

RUN_SPCA_CV <- function(dat, nfolds=3, ridge_seq=lseq(from=0.000000001,to=0.01,length.out=8), lasso_seq=lseq(from=0.01,to=0.01, length.out=1)){

  param_combs <- expand.grid(seq(1, nfolds), ridge_seq, lasso_seq)
  
  res_spca_cv <- parApply(cl = cl, X = param_combs, MARGIN = 1, FUN = spca_cv_1, dat = dat, nfolds=nfolds)
  
  res_spca_cv <- matrix(res_spca_cv, nrow = 4)
  
  df_res_spca_cv <- data.frame(fold=res_spca_cv[1,],ridge=res_spca_cv[2,],lasso=res_spca_cv[3,],comb_MSE=res_spca_cv[4,])
  
  df_res_spca_cv_with_sd <- df_res_spca_cv %>%
    group_by(ridge, lasso) %>%
    summarise(stdError = sd(comb_MSE))
  df_res_spca_cv_with_mean <- df_res_spca_cv %>%
    group_by(ridge, lasso) %>%
    summarise(mean_comb_MSE = mean(comb_MSE))
  
  df_res_spca_cv <- cbind(df_res_spca_cv_with_mean, df_res_spca_cv_with_sd[,3])
  
  # ------------- HEATMAP OF PARAMETER COMBINATION MSE -----------------------------
  # create heatmap to visually find optimal values for l1 and l2
  # df_heatmap <- df_res_spca_cv %>%
  #   mutate(ridge = factor(formatC(ridge, format = "e", digits = 2)),
  #          lasso = factor(formatC(lasso, format = "e", digits = 2)))
  #
  # p <- ggplot(df_heatmap, aes(x = ridge, lasso)) +
  #   geom_tile(aes(fill = comb_MSE)) +
  #   #geom_text(aes(label = round(comb_MSE, 2))) +
  #   scale_fill_gradient(low = "cyan3", high = "coral1") +
  #   labs(fill = "MSE") +
  #   theme(text = element_text(size=20))

  # ggsave(
  #   filename = "parameter_tuning_"n"_"p"_"error".png",
  #   path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
  #   scale = 1,
  #   dpi = 600,
  #   height = 8,
  #   width = 10,
  #   limitsize = TRUE)
  # ---------------------------------------------------------------------
  
  MSE <- df_res_spca_cv$mean_comb_MSE
  MSEstdError <- df_res_spca_cv$stdError
  
  bestModel <- which.min(MSE)
  MSE1stdErrorRule <- MSE[MSE < (MSE[bestModel] + MSEstdError[bestModel])]
  bestModel1stdErrorRule <- which(MSE1stdErrorRule[length(MSE1stdErrorRule)] == MSE)
  
  row_min <- df_res_spca_cv[bestModel1stdErrorRule,]

  best_spca_res <- spca(dat, k =  min(dim(dat))-1, alpha = row_min$ridge, beta = row_min$lasso)
  best_spca_loadings <- best_spca_res$loadings

  return(list(MSE_matrix = df_res_spca_cv, L1 = row_min$ridge, L2 = row_min$lasso, mse = row_min$mean_comb_MSE, best_spca = best_spca_res))
}
