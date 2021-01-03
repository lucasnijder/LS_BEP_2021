OC_method <- function(model_eigenval){
  num_of_comp <- nScree(model_eigenval)[1]$Components$noc
  
  return(num_of_comp)
}