set.seed(2)

variances <- makeVariance(varianceOfComps = rep(10,3), p = 20, error = 0.1)
dat <- makeDat(n = 100, p = 20, ncomp = 3, comdis = NULL, variances = variances)$X

ev <- prcomp(dat)$sdev^2
ev_df <- data.frame(ev=ev, index=1:10)

ggplot(ev_df, aes(x=index, y=ev)) + geom_bar(stat='identity') +
  xlab('component') +
  ylab('eigenvalue') +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(breaks = 0:ceiling(max(ev))) +
  geom_hline(aes(yintercept=1), colour='red') +
  theme(text = element_text(size=12)) 

ggsave(
  filename = "screeplot_example.png",
  path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600,
  limitsize = TRUE)