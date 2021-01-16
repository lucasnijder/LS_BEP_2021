library(tictoc)
library(sparsepca)
library(parallel)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(MASS)
library(ggpubr)
library(grid)
library(sjmisc)
library(nFactors)
library(ggthemes)
library(glue)
library(ggtext)

source('PA_method.R')
source('AF_method.R')
source('OC_method.R')
source('KG_method.R')
source('CV_method.R')
source('Schipper_Deun_make_data.R')
source('SPCA_cv.R')

cl <- makeCluster(7)
clusterExport(cl, 'spca')

n_seq <- c(25, 50, 100)
c_seq <- 3
p_seq <- c(10, 20)
error_seq <- c(0.01,0.1,0.15)

num_mc <- 50
run_type <- 'spca'

combs <- expand.grid(n_seq, c_seq, p_seq, error_seq)

# --------------------------------------- #

tic('start')


run_sim <- function(comb, type){
  
  n <- comb[1]
  c <- comb[2]
  p <- comb[3]
  error <- comb[4]
  
  variances <- makeVariance(varianceOfComps = rep(10,c), p = p, error = error)
  dat <- makeDat(n = n, p = p, ncomp = c, comdis = NULL, variances = variances)$X
   
  if(type == 'pca'){
    pca_res <- prcomp(dat)
    pca_eigenvalues <- pca_res$sdev^2
    
    # eigenvalue pattern methods
    PA_predicted_comps <- PA_method(dat, pca_eigenvalues, 'pca')
    OC_predicted_comps <- OC_method(pca_eigenvalues)
    AF_predicted_comps <- AF_method(pca_eigenvalues)
    KG_predicted_comps <- KG_method(pca_eigenvalues)
    
    # structural model methods
    CV_predicted_comps <- CV_method(dat, 'pca')
  }
  
  if(type == 'spca'){
    res <- RUN_SPCA_CV(dat)
    spca_res <- res[[5]]
    spca_eigenvalues <- spca_res$eigenvalues
    
    # eigenvalue pattern methods
    PA_predicted_comps <- PA_method(dat, spca_eigenvalues, 'spca')
    OC_predicted_comps <- OC_method(spca_eigenvalues)
    AF_predicted_comps <- AF_method(spca_eigenvalues)
    KG_predicted_comps <- KG_method(spca_eigenvalues)
    
    # Because of the problem with the tuned parameters in the CV method we will use the default L1 and L2
    L1 <- res[[2]]
    L2 <- res[[3]]
    
    # structural model methods
    CV_predicted_comps <- CV_method(dat, 'spca', L1=L1, L2=L2)
  }
  
  return(list(n, p, c, error, PA_predicted_comps, OC_predicted_comps, AF_predicted_comps, KG_predicted_comps, CV_predicted_comps))
}

# -------------------------------- #

Process_res <- function(type, pred_comps, res){
  
  if(type == 'PA'){
    s = 1
  }
  
  if(type == 'OC'){
    s = 2
  }
  
  if(type == 'AF'){
    s = 3
  }
  
  if(type == 'KG'){
    s = 4
  }
  
  if(type == 'CV'){
    s = 5
  }
  
  # if more techniques are added, change 5 to the total number of techniques
  pred_comps <- pred_comps[,seq(s, ncol(pred_comps), 5)]
  
  percentage_correct <- rowSums(pred_comps == res[,2])*100/num_mc
  
  res <- cbind(res, percentage_correct)
  # this 4 represents the 4 columns of parameters, so do not change if new techniques are added
  colnames(res)[4+s] <- sprintf('%s_percentage_correct', type)
  
  type_pred_comps <- matrix(0, nrow = 18, ncol = 9)
  
  for(i in 0:8){
    type_pred_comps[,(i+1)] <- rowSums(pred_comps == i)
  }
  
  colnames(type_pred_comps) <- c('PC0','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8')
  
  return(list(type_pred_comps, res))
}

# ----------------------------------- #

Run_MC <- function(num_mc){
  res <- combs
  for(i in 1:num_mc){
    
    set.seed(i+100)
    
    sim_res <- apply(combs, 1, run_sim, type = run_type)
    
    predicted_comps <- t(matrix(unlist(sim_res), nrow=9))[,5:9]
    
    res <- cbind(res, predicted_comps)
  }
  
  pred_comps <- res[,-(1:4)]
  res <- res[,1:4]
  
  run_PA <- Process_res('PA', pred_comps, res)
  PA_pred_comps <- run_PA[[1]]
  res <- run_PA[[2]]
  
  run_OC <- Process_res('OC', pred_comps, res)
  OC_pred_comps <- run_OC[[1]]
  res <- run_OC[[2]]

  run_AF <- Process_res('AF', pred_comps, res)
  AF_pred_comps <- run_AF[[1]]
  res <- run_AF[[2]]
  
  run_KG <- Process_res('KG', pred_comps, res)
  KG_pred_comps <- run_KG[[1]]
  res <- run_KG[[2]]
  
  run_CV <- Process_res('CV', pred_comps, res)
  CV_pred_comps <- run_CV[[1]]
  res <- run_CV[[2]]
  
  return(list(res, PA_pred_comps, OC_pred_comps, AF_pred_comps, KG_pred_comps, CV_pred_comps))
}

# ---------------------------------- #

res_mc <- Run_MC(num_mc)
res <- res_mc[[1]]
PA_pred_comps <- res_mc[[2]]
OC_pred_comps <- res_mc[[3]]
AF_pred_comps <- res_mc[[4]]
KG_pred_comps <- res_mc[[5]]
CV_pred_comps <- res_mc[[6]]

write.table(res, sprintf('%s_res.csv',run_type), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)
write.table(PA_pred_comps, sprintf('%s_PA_pred_comps.csv',run_type), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)
write.table(OC_pred_comps, sprintf('%s_OC_pred_comps.csv',run_type), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)
write.table(AF_pred_comps, sprintf('%s_AF_pred_comps.csv',run_type), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)
write.table(KG_pred_comps, sprintf('%s_KG_pred_comps.csv',run_type), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)
write.table(CV_pred_comps, sprintf('%s_CV_pred_comps.csv',run_type), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)

toc()

stopCluster(cl=cl)

