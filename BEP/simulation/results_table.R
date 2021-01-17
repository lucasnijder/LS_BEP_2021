library(data.table)
library(xtable)

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

for(type in c('pca','spca')){
  res <- read.csv(sprintf('%s_res.csv',type))
  
  res <- res[,-2]
  res[,3] <- res[,3] * 100
  setnames(res, old = c('Var1','Var3','Var4',
                        'PA_percentage_correct','OC_percentage_correct',
                        'AF_percentage_correct','KG_percentage_correct',
                        'CV_percentage_correct'), new = c('n','p','error','PA','OC','AF','KG','CV'))
  res <- res[c('error','p','n','KG','OC','AF','PA','CV')]
  res[19,] <- colMeans(res)
  res[19,1:3] <- c('','','')
  print(xtable(res, digits=c(0,0,0,0,0,0,0,0,0)), include.rownames=FALSE)
}

