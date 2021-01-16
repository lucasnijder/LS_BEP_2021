
Plot_comp_predictions <- function(pred_comps, type, error){
  if(error == 0.01){
    i <- 1
  }
  if(error == 0.1){
    i <- 7
  }
  if(error == 0.15){
    i <- 13
  }
  
  p_list <-  vector('list', 6)
  for(a in 0:5){
    bar_dat <- melt(pred_comps[i+a,])
    bar_dat <- cbind(variable = rownames(bar_dat), bar_dat)
    bar_dat <- bar_dat %>% mutate(correct_comp = ifelse(variable == 'PC3', T, F)) # value == bar_dat[4,2]
    bar_dat <- bar_dat %>% mutate(value = value*100/num_mc)
    
    p <- ggplot(data=bar_dat, aes(x=variable, y=value)) +
      geom_bar(stat="identity", aes(fill = correct_comp))+
      scale_fill_manual(values = c('coral1', 'cyan3')) + 
      xlab("number of PCs") + ylab("% of outcomes") + 
      ylim(0,100) +
      theme(axis.title.x = element_text(size=8),
            axis.title.y = element_text(size=8),
            plot.margin = unit(c(1,1,1,1), "lines"),
            legend.position = 'none')
    
    p_list[[a+1]] <- p
  }
  
  return(ggarrange(p_list[[1]],
                   p_list[[2]],
                   p_list[[3]],
                   p_list[[4]],
                   p_list[[5]],
                   p_list[[6]],
                   labels = c(sprintf("              A.   n = 25 and p = 10"),
                              sprintf("              B.   n = 50 and p = 10"),
                              sprintf("              C.   n = 100 and p = 10"),
                              sprintf("              D.   n = 25 and p = 20"),
                              sprintf("              E.   n = 50 and p = 20"),
                              sprintf("              F.   n = 100 and p = 20")),
                   ncol = 3, nrow = 2,
                   font.label = list(size = 10,
                                     color = "black"),
                   hjust = 0,
                   vjust = 1))
}

print(sprintf("PA = %s, OC = %s, AF = %s, KG = %s, CV = %s", 
              round(mean(res$PA_percentage_correct),2), 
              round(mean(res$OC_percentage_correct),2),
              round(mean(res$AF_percentage_correct),2),
              round(mean(res$KG_percentage_correct),2),
              round(mean(res$CV_percentage_correct),2)))

count = 0
for(type in list(PA_pred_comps, OC_pred_comps, AF_pred_comps, KG_pred_comps, CV_pred_comps)){
  if(count == 0){
    type_str = 'PA'
  }
  if(count == 1){
    type_str = 'OC'
  }
  if(count == 2){
    type_str = 'AF'
  }
  if(count == 3){
    type_str = 'KG'
  }
  if(count == 4){
    type_str = 'CV'
  }
  
  for(error in list(0.01,0.1,0.15)){
    Plot_comp_predictions(type, toupper(run_type), error)
    
    ggsave(
      sprintf(sprintf("%s_%s_%s.png", toupper(run_type), type_str, error)),
      path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
      scale = 1,
      width = 10, 
      height = 6,
      units = "in",
      dpi = 600,
      limitsize = TRUE)
  }
  count <- count + 1
}
