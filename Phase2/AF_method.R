AF_method <- function(model_eigenval){
  num_of_comp <- nScree(model_eigenval)[1]$Components$naf
  
  return(num_of_comp)
}