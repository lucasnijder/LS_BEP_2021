AF_method <- function(model_eigenval){
  af_list <- c()
  for(j in 1:length(model_eigenval)){
    af = model_eigenval[j+1] - 2* model_eigenval[j] - model_eigenval[j-1]
    af_list <- append(af_list, af)
  }
  num_of_comp <- which.max(af_list)
  return(num_of_comp)
}