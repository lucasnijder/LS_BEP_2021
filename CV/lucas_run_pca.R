# call to run PCA, specify which function is to be used (facto/other)
RUN_PCA <- function(dat){

    prcomp_pca_res <- prcomp(dat, rank. = min(dim(dat))-1)
    
    return(prcomp_pca_res)
}