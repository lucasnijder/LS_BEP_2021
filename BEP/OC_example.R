source("deun_makeData.R")

set.seed(4)

variances <- makeVariance(varianceOfComps = rep(10,3), p = 10, error = 0.01)
dat <- makeDat(n = 50, p = 10, ncomp = 3, comdis = NULL, variances = variances)$X

ev <- prcomp(dat)$sdev^2
x_coord <- c(3, 10)
y_coord <- c(ev[3] , ev[10])

fit <- lm(y_coord ~ x_coord)
intercept <- fit$coefficients[1]
slope <- fit$coefficients[2]

ev_df <- data.frame(ev=ev, index=1:10)

ggplot(ev_df) + geom_line(aes(x=index, y=ev), size=0.8) + 
  geom_point(aes(x=index, y=ev), size=2) +
  geom_point(aes(x=2, y=intercept+2*slope), size=2, color='red') + 
  geom_abline(intercept = intercept, slope=slope, color='red', size=0.8) +
  xlab('component') +
  ylab('eigenvalue') +
  scale_y_continuous(breaks = 1:ceiling(max(ev))) +
  scale_x_continuous(breaks = 1:10) +
  theme(text = element_text(size=12))

ggsave(
  filename = "OC_example.png",
  path = "C:/Users/20175878/Documents/bep_with_version_control/figs",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600,
  limitsize = TRUE)