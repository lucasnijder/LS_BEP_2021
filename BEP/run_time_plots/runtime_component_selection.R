source('PA_method.R')
source('AF_method.R')
source('OC_method.R')
source('KG_method.R')
source('CV_method.R')
source('Schipper_Deun_make_data.R')
source('SPCA_cv.R')

library(microbenchmark)

n <- 100
p <- 10
c <- 3
error <- 0.1

variances <- makeVariance(varianceOfComps = rep(10,c), p = p, error = error)
dat <- makeDat(n = n, p = p, ncomp = c, comdis = NULL, variances = variances)$X

spca_res <- spca(dat)
spca_ev <- spca_res$eigenvalues

pca_res <- prcomp(dat)
pca_ev <- pca_res$sdev^2

mbm <- microbenchmark("KG" = KG_method(pca_ev),
                      "OC" = OC_method(pca_ev),
                      "AF" = AF_method(pca_ev),
                      "PA" = PA_method(dat, pca_ev, 'pca'),
                      "CV" = CV_method(dat, 'pca'))

autoplot(mbm, log = T) + theme_grey()

ggsave(
  "runtime_component_selection_pca.png",
  path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
  scale = 1,
  dpi = 600,
  limitsize = TRUE)