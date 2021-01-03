KG_method <- function(model_eigenval){
  res <- model_eigenval >= 1
  num_of_comps <- sum(res)

  return(num_of_comps)
}