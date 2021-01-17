set.seed(1)

variances <- makeVariance(varianceOfComps = rep(10,3), p = 10, error = 0.1)
dat <- makeDat(n = 100, p = 10, ncomp = 3, comdis = NULL, variances = variances)$X

pa_res <- fa.parallel(dat, fa = 'pc')

ev <- pa_res$pc.values
ev_sim <- pa_res$pc.sim
ev_df <- data.frame(index=1:10, ev=ev, ev_sim=ev_sim)

ggplot(ev_df) + geom_line(aes(x=index, y=ev), size=0.8) + 
  geom_point(aes(x=index, y=ev), size=2) +
  geom_line(aes(x=index, y=ev_sim), size=0.8, color='red') + 
  xlab('component') +
  ylab('eigenvalue') +
  scale_y_continuous(breaks = 1:ceiling(max(ev))) +
  scale_x_continuous(breaks = 1:10) +
  theme(text = element_text(size=12))

ggsave(
  filename = "PA_example.png",
  path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600,
  limitsize = TRUE)