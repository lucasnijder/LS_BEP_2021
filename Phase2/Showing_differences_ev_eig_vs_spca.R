set.seed(1)

n_seq <- c(25, 50, 100)
c_seq <- 3
p_seq <- c(10, 20)
error_seq <- c(0.01,0.1,0.15)

combs <- expand.grid(n_seq, c_seq, p_seq, error_seq)

Eigenval_difference <- function(comb, type){
  n <- comb[1]
  c <- comb[2]
  p <- comb[3]
  error <- comb[4]
  
  variances <- makeVariance(varianceOfComps = rep(10,c), p = p, error = error)
  dat <- makeDat(n = n, p = p, ncomp = c, comdis = NULL, variances = variances)$X
  
  eigen_ev <- eigen(cor(dat))$values
  spca_ev <- spca(dat)$eigenvalues
 
  return(list(eigen_ev = eigen_ev[(1:8)], spca_ev = spca_ev[(1:8)]))
}

ev_sim_res <- apply(combs, 1, Eigenval_difference, type = 'pca')

ev_sim_res_matrix <- matrix(unlist(ev_sim_res), nrow=8)
rownames(ev_sim_res_matrix) <- c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8')
spca_ev <- ev_sim_res_matrix[,seq_len(ncol(ev_sim_res_matrix)) %% 2 == 0]
eigen_ev <- ev_sim_res_matrix[,seq_len(ncol(ev_sim_res_matrix)) %% 2 != 0]

ev_df <- data.frame(id = 1:18,type = rep('SPCA eigenvalue',18), t(spca_ev))
ev_df <- rbind(ev_df, data.frame(id = 19:36, type = rep('Eigen eigenvalue',18), t(eigen_ev)))
df_melted <- melt(ev_df, id.vars=c('id','type'))

ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = type, group = id), size=1) +
  xlab('component') + ylab('eigenvalue') + theme_grey()

ggsave(
  "diff_ev_eigen_spca.png",
  path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
  scale = 1,
  dpi = 600,
  limitsize = TRUE)
