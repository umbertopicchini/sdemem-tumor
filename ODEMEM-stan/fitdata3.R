

source('indata.R')

mydat3 <- get.my.data('datagroup3')
plot.my.data(dat=mydat3, cols=my.cols1)

library(rstan)
# ENABLE MULTI CORE CALCULATIONS WITH THE FOLLOWING TWO LINES
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())  # detect number of cores on the current processor

y <- mydat3$frame$logvol
x <- mydat3$frame$time-min(mydat3$frame$time)
N <- mydat3$N.sub
ns <- mydat3$n.obs
t <- cumsum(c(0,ns))
M <- sum(ns)
# BRUTAL HARDCODING! THIS IS THE FIRST MEASURED LOVOLUME FOR EACH SUBJECT
logv0 <- c(6.642617,5.109575,4.709530,3.775057,5.655292,5.522660,6.372978,6.622204)

# STARTING PARAMETER VALUES
initpar <- function() {
  list(mu_alpha = 0.5, mu_logbeta = 1.39, mu_logdelta = 0, sigma_alpha =0.3, sigma_beta = 0.5, sigma_delta =0.5, sigma_epsilon =0.1)
} 

my.pars <- c("mu_alpha", "mu_beta", "mu_delta", "sigma_alpha","sigma_beta","sigma_delta","sigma_epsilon")
fit <- stan(file = 'tumorgrowth.stan', init=initpar, iter = 10000, control = list(adapt_delta = 0.95))

# TRACEPLOTS, with and without warm-up:
traceplot(fit, pars = my.pars,inc_warmup = TRUE)  # plot traces with burnin
traceplot(fit, pars = my.pars,inc_warmup = FALSE)  # plot traces without burnin

fit0 <- as.matrix(fit, pars = my.pars)
write.table(fit0, 'newfit3.txt', row.names=FALSE)


pairs(fit, pars = my.pars)
# posterior inference analytics
print(fit, digits=2, pars = my.pars)
# MARGINALS and correlations
pairs(fit,pars = my.pars)

