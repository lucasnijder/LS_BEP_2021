copy_of_function_spca <- function (X, k = NULL, alpha = 1e-04, beta = 1e-04, center = TRUE, 
          scale = FALSE, max_iter = 1000, tol = 1e-05, verbose = TRUE) 
{
  X <- as.matrix(X)
  if (any(is.na(X))) {
    warning("Missing values are omitted: na.omit(X).")
    X <- stats::na.omit(X)
  }
  spcaObj = list(u = NULL, v = NULL, d = NULL) # loadings = NULL, transform = NULL, scores = NULL, eigenvalues = NULL, center = center, scale = scale,
  
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(k)) 
    k <- min(n, p)
  if (k > min(n, p)) 
    k <- min(n, p)
  if (k < 1) 
    stop("Target rank is not valid!")
  if (center == TRUE) {
    spcaObj$center <- colMeans(X)
    X <- sweep(X, MARGIN = 2, STATS = spcaObj$center, FUN = "-", 
               check.margin = TRUE)
  }
  else {
    spcaObj$center <- FALSE
  }
  if (scale == TRUE) {
    spcaObj$scale <- sqrt(colSums(X^2)/(n - 1))
    if (is.complex(spcaObj$scale)) {
      spcaObj$scale[Re(spcaObj$scale) < 1e-08] <- 1 + 
        (0+0i)
    }
    else {
      spcaObj$scale[spcaObj$scale < 1e-08] <- 1
    }
    X <- sweep(X, MARGIN = 2, STATS = spcaObj$scale, FUN = "/", 
               check.margin = TRUE)
  }
  else {
    spcaObj$scale <- FALSE
  }
  svd_init <- svd(X)
  Dmax <- svd_init$d[1]
  A <- svd_init$v[, 1:k]
  B <- svd_init$v[, 1:k]
  V <- svd_init$v
  VD = sweep(V, MARGIN = 2, STATS = svd_init$d, FUN = "*", 
             check.margin = TRUE)
  VD2 = sweep(V, MARGIN = 2, STATS = svd_init$d^2, FUN = "*", 
              check.margin = TRUE)
  alpha <- alpha * Dmax^2
  beta <- beta * Dmax^2
  nu <- 1/(Dmax^2 + beta)
  kappa <- nu * alpha
  obj <- c()
  improvement <- Inf
  noi <- 1
  while (noi <= max_iter && improvement > tol) {
    Z <- VD2 %*% (t(V) %*% B)
    svd_update <- svd(Z)
    A <- svd_update$u %*% t(svd_update$v)
    grad <- VD2 %*% (t(V) %*% (A - B)) - beta * B
    B_temp <- B + nu * grad
    idxH <- which(B_temp > kappa)
    idxL <- which(B_temp <= -kappa)
    B = matrix(0, nrow = nrow(B_temp), ncol = ncol(B_temp))
    B[idxH] <- B_temp[idxH] - kappa
    B[idxL] <- B_temp[idxL] + kappa
    R <- t(VD) - (t(VD) %*% B) %*% t(A)
    obj <- c(obj, 0.5 * sum(R^2) + alpha * sum(abs(B)) + 
               0.5 * beta * sum(B^2))
    if (noi > 1) {
      improvement <- (obj[noi - 1] - obj[noi])/obj[noi]
    }
    if (verbose > 0 && (noi - 1)%%10 == 0) {
      print(sprintf("Iteration: %4d, Objective: %1.5e, Relative improvement %1.5e", 
                    noi, obj[noi], improvement))
    }
    noi <- noi + 1
  }

  # problem: find an U that is 41x10
  
  spcaObj$u <- u
  spcaObj$v <- svd_update$v
  spcaObj$d <- svd_update$d
  
  # spcaObj$loadings <- B
  # spcaObj$transform <- A
  # spcaObj$scores <- X %*% B
  # spcaObj$eigenvalues <- svd_update$d/(n - 1)
  # spcaObj$objective <- obj
  # spcaObj$sdev <- sqrt(spcaObj$eigenvalues)
  # spcaObj$var <- sum(apply(Re(X), 2, stats::var))
  # if (is.complex(X))
  #   spcaObj$var <- Re(spcaObj$var + sum(apply(Im(X), 2,
  #                                             stats::var)))
  class(spcaObj) <- "copy_of_function_spca"
  return(spcaObj)
}
