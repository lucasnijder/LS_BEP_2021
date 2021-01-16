source('copy_of_function_spca.R')
source('copy_of_function_estim_ncp_svd.R')
source('copy_of_function_estim_ncp_spca.R')

res_svd <- copy_of_function_estim_ncp_svd(decathlon)
res_spca <- copy_of_function_estim_ncp_spca(decathlon)