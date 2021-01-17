highlight = function(x, pat, color="black", family="") {
  x <- substr(x, 3, 3)
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

Plot_per_error <- function(total_res, error){
  
  total_res <- total_res %>% arrange(index)
  total_res <- total_res[,-11]
  
  if(error == 0.01){
    i <- 1
  }
  if(error == 0.1){
    i <- 13
  }
  if(error == 0.15){
    i <- 25
  }
  

  p_list <-  vector('list', 6)
  for(a in seq(0,12,2)){
    bar_dat <- melt(total_res[(i+a):(i+a+1),])
    
    bar_dat <- bar_dat %>% mutate(value = value*100/num_mc)
    
    p <- ggplot(data=bar_dat, aes(x=variable, y=value, fill=type)) +
      geom_bar(stat="identity", position = 'dodge')+
      scale_fill_manual(values = c('coral1', 'cyan3')) +
      scale_x_discrete(labels= function(x) highlight(x, "3", "black")) + 
      xlab("number of PCs") + ylab("% of outcomes") +
      ylim(0,100) +
      theme(axis.title.x = element_text(size=8),
            axis.title.y = element_text(size=8),
            plot.margin = unit(c(1,1,1,1), "lines"),
            legend.position = 'none',
            axis.text.x=element_markdown())
    
    p_list[[(a/2)+1]] <- p
  }
  
  return(ggarrange(p_list[[1]],
                   p_list[[2]],
                   p_list[[3]],
                   p_list[[4]],
                   p_list[[5]],
                   p_list[[6]],
                   labels = c("              A.   n = 25 and p = 10",
                              "              B.   n = 50 and p = 10",
                              "              C.   n = 100 and p = 10",
                              "              D.   n = 25 and p = 20",
                              "              E.   n = 50 and p = 20",
                              "              F.   n = 100 and p = 20"),
                   ncol = 3, nrow = 2,
                   font.label = list(size = 10,
                                     color = "black"),
                   hjust = 0,
                   vjust = 1,
                   common.legend = T))
}

Plot_but_not_per_n <- function(total_res){

  total_res <- total_res[,-11]

  p <- rep(c(10,10,10,20,20,20),3)
  n <- rep(c(25,50,100),6)
  error <- c(rep(0.01,6), rep(0.1,6), rep(0.15,6))
  total_res <- cbind(total_res, p)
  total_res <- cbind(total_res, n)
  total_res <- cbind(total_res, error)
  
  total_res <- total_res %>% group_by(error, p, type) %>% summarise(PC0=sum(PC0),PC1=sum(PC1),PC2=sum(PC2),
                                                         PC3=sum(PC3),PC4=sum(PC4),PC5=sum(PC5),
                                                         PC6=sum(PC6),PC7=sum(PC7),PC8=sum(PC8))
  

  p_list <-  vector('list', 6)
  for(a in seq(1,12,2)){
    bar_dat <- melt(total_res[a:(a+1),], id.vars = c('type','error','p'))
    
    bar_dat <- bar_dat %>% mutate(value = round((value*100)/(num_mc*3),2))

    p <- ggplot(data=bar_dat, aes(x=variable, y=value, fill=type)) +
      geom_bar(stat="identity", position = 'dodge')+
      scale_fill_manual(values = c('coral1', 'cyan3')) +
      scale_x_discrete(labels= function(x) highlight(x, "3", "black")) +
      xlab("number of PCs") + ylab("% of outcomes") +
      ylim(0,100) +
      theme(axis.title.x = element_text(size=8),
            axis.title.y = element_text(size=8),
            plot.margin = unit(c(1,1,1,1), "lines"),
            legend.position = 'none',
            axis.text.x=element_markdown())

    p_list[[(a+1)/2]] <- p
  }

  return(ggarrange(p_list[[1]],
                   p_list[[4]],
                   p_list[[5]],
                   p_list[[2]],
                   p_list[[3]],
                   p_list[[6]],
                   labels = c("              A.   noise = 1% and p = 10",
                              "              B.   noise = 10% and p = 10",
                              "              C.   noise = 15% and p = 10",
                              "              D.   noise = 1% and p = 20",
                              "              E.   noise = 10% and p = 20",
                              "              F.   noise = 15% and p = 20"),
                   ncol = 3, nrow = 2,
                   font.label = list(size = 10,
                                     color = "black"),
                   hjust = 0,
                   vjust = 1,
                   common.legend = T))
}

Plot_PCA_SPCA_res <- function(){
  for(type_str in list('PA', 'OC', 'AF', 'KG', 'CV')){
    
    pca_res <- read.csv(file = sprintf('pca_%s_pred_comps.csv',type_str))
    spca_res <- read.csv(file = sprintf('spca_%s_pred_comps.csv',type_str))
    
    pca_res <- cbind(pca_res, type = rep('pca',nrow(pca_res)))
    pca_res <- cbind(pca_res, index=seq(1:18))

    spca_res <- cbind(spca_res, type = rep('spca',nrow(pca_res)))
    spca_res <- cbind(spca_res, index=seq(1:18))

    total_res <- rbind(pca_res, spca_res)
    
    for(error in list(0.01,0.1,0.15)){
      Plot_per_error(total_res, error)

      ggsave(
        sprintf(sprintf("%s_%s.png", type_str, error)),
        path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
        scale = 1,
        width = 14,
        height = 6,
        units = "in",
        dpi = 600,
        limitsize = TRUE)
    }
    
    Plot_but_not_per_n(total_res)
    ggsave(
      sprintf(sprintf("%s_overall.png", type_str)),
      path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
      scale = 1,
      width = 14,
      height = 6,
      units = "in",
      dpi = 600,
      limitsize = TRUE)
    
  }
}

num_mc <- 50

Plot_PCA_SPCA_res()