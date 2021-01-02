n_seq <- c(25, 50, 100)
c_seq <- 3
p_seq <- c(10, 20)
error_seq <- c(0.01,0.1,0.2)

num_mc <- 1

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
 
  return(list(eigen_ev = eigen_ev, spca_ev = spca_ev))
}

ev_sim_res <- apply(combs, 1, Eigenval_difference, type = 'pca')
