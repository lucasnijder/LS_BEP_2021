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

cl <- makeCluster(7)
clusterExport(cl, 'spca')

n_seq <- c(25, 50, 100)
c_seq <- 3
p_seq <- c(10, 20)
error_seq <- c(0.01,0.1,0.15)

num_mc <- 100

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
    
    PA_predicted_comps <- ParallelAnalysis(dat, pca_eigenvalues, 'pca')
    OC_predicted_comps <- OC_method(pca_eigenvalues)
    AF_predicted_comps <- AF_method(pca_eigenvalues)
    KG_predicted_comps <- KG_method(pca_eigenvalues)
  }
  
  if(type == 'spca'){
    spca_res <- RUN_SPCA_CV(dat)[[5]]
    spca_eigenvalues <- spca_res$eigenvalues
    
    PA_predicted_comps <- ParallelAnalysis(dat, spca_eigenvalues, 'spca')
    OC_predicted_comps <- OC_method(spca_eigenvalues)
    AF_predicted_comps <- AF_method(spca_eigenvalues)
    KG_predicted_comps <- KG_method(spca_eigenvalues)
  }
  
  return(list(n, p, c, error, PA_predicted_comps, OC_predicted_comps, AF_predicted_comps, KG_predicted_comps))
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
  
  # if more techniques are added, change 4 to the total number of techniques
  pred_comps <- cbind(pred_comps[,seq(s, ncol(pred_comps), 4)])
  
  percentage_correct <- rowSums(pred_comps == res[,2])*100/num_mc
  
  res <- cbind(res, percentage_correct)
  colnames(res)[4+s] <- sprintf('%s_percentage_correct', type)

  for(i in 0:8){
    pred_comps <- row_count(pred_comps, count = i, var = sprintf("PC%s",i))
  }
  
  pred_comps <- data.frame(pred_comps[,-(1:num_mc)])
  
  return(list(pred_comps, res))
}

# ----------------------------------- #

Run_MC <- function(num_mc=10){
  res <- combs
  for(i in 1:num_mc){
    
    set.seed(i)
    
    sim_res <- apply(combs, 1, run_sim, type = 'pca')
    
    predicted_comps <- t(matrix(unlist(sim_res), nrow=8))[,5:8]
    
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
  
  return(list(res, PA_pred_comps, OC_pred_comps, AF_pred_comps, KG_pred_comps))
}

# ---------------------------------- #

res_mc <- Run_MC(num_mc)
res <- res_mc[[1]]
PA_pred_comps <- res_mc[[2]]
OC_pred_comps <- res_mc[[3]]
AF_pred_comps <- res_mc[[4]]
KG_pred_comps <- res_mc[[5]]

Plot_comp_predictions <- function(pred_comps, type, error){
  if(error == 0.01){
    i <- 1
  }
  if(error == 0.1){
    i <- 7
  }
  if(error == 0.15){
    i <- 13
  }
  
  p_list <-  vector('list', 6)
  for(a in 0:5){
    bar_dat <- melt(pred_comps[i+a,])
    bar_dat <- bar_dat %>% mutate(correct_comp = ifelse(value == bar_dat[4,2], T, F))
    
     p <- ggplot(data=bar_dat, aes(x=variable, y=value)) +
      geom_bar(stat="identity", aes(fill = correct_comp))+
      scale_fill_manual(values = c('gray30', 'green')) + 
      xlab("number of PCs") + ylab("% of outcomes") + 
      ylim(0,100) +
      theme(axis.title.x = element_text(size=8),
            axis.title.y = element_text(size=8),
            plot.margin = unit(c(1,1,1,1), "lines"),
            legend.position = 'none')
    
    correct_comps_df <- data.frame(PC3=pred_comps$PC3)
    
    p
    
    p_list[[a+1]] <- p
  }
  
  return(ggarrange(p_list[[1]],
            p_list[[2]],
            p_list[[3]],
            p_list[[4]],
            p_list[[5]],
            p_list[[6]],
            labels = c(sprintf("              A.   n = 25 and p = 10"),
                       sprintf("              B.   n = 50 and p = 10"),
                       sprintf("              C.   n = 100 and p = 10"),
                       sprintf("              D.   n = 25 and p = 20"),
                       sprintf("              E.   n = 50 and p = 20"),
                       sprintf("              F.   n = 100 and p = 20")),
            ncol = 3, nrow = 2,
            font.label = list(size = 10,
                              color = "black"),
            hjust = 0,
            vjust = 1))
}

count = 0
for(type in list(PA_pred_comps, OC_pred_comps, AF_pred_comps, KG_pred_comps)){
  if(count == 0){
    type_str = 'PA'
  }
  if(count == 1){
    type_str = 'OC'
  }
  if(count == 2){
    type_str = 'AF'
  }
  if(count == 3){
    type_str = 'KG'
  }
  
  for(error in list(0.01,0.1,0.15)){
    Plot_comp_predictions(type, 'PCA', error)

    ggsave(
      sprintf(sprintf("%s_%s_%s.png", 'PCA', type_str, error)),
      path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
      scale = 1,
      width = 10, 
      height = 6,
      units = "in",
      dpi = 600,
      limitsize = TRUE)
  }
  count <- count + 1
}

print(sprintf("PA = %s, OC = %s, AF = %s, KG = %s", round(mean(res$PA_percentage_correct),2), 
              round(mean(res$OC_percentage_correct),2),
              round(mean(res$AF_percentage_correct),2),
              round(mean(res$KG_percentage_correct),2)))

stopCluster(cl=cl)

toc()
