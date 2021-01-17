library(microbenchmark)

n <- 100
p <- 10
c <- 3
error <- 0.1

cl <- makeCluster(7)
clusterExport(cl, 'spca')

variances <- makeVariance(varianceOfComps = rep(10,c), p = p, error = error)
dat <- makeDat(n = n, p = p, ncomp = c, comdis = NULL, variances = variances)$X

mbm <- microbenchmark("PCA" = prcomp(dat),
                      "SPCA" = spca(dat),
                      "SPCA with CV" = RUN_SPCA_CV(dat))

autoplot(mbm, log = T) + theme_grey()

ggsave(
  "runtime_models.png",
  path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
  scale = 1,
  dpi = 600,
  limitsize = TRUE)