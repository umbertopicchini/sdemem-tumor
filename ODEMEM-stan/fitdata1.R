

source('indata.R')

mydat1 <- get.my.data('datagroup1')
plot.my.data(dat=mydat1, cols=my.cols1)

library(rstan)
# ENABLE MULTI CORE CALCULATIONS WITH THE FOLLOWING TWO LINES
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())  # detect number of cores on the current processor

y <- mydat1$frame$logvol
x <- mydat1$frame$time-min(mydat3$frame$time)
N <- mydat1$N.sub
ns <- mydat1$n.obs
t <- cumsum(c(0,ns))
M <- sum(ns)
# BRUTAL HARDCODING! THIS IS THE FIRST MEASURED LOVOLUME FOR EACH SUBJECT
logv0 <- c(5.5033, 4.6032, 5.8808, 5.6781, 6.2943)

# STARTING PARAMETER VALUES
initpar <- function() {
  list(mu_alpha = 0.5, mu_logbeta = 1.39, mu_logdelta = 0, sigma_alpha =0.3, sigma_beta = 0.5, sigma_delta =0.5, sigma_epsilon =0.1)
} 

my.pars <- c("mu_alpha", "mu_beta", "mu_delta", "sigma_alpha","sigma_beta","sigma_delta","sigma_epsilon")
fit <- stan(file = 'tumorgrowth.stan', init=initpar, iter = 10000, control = list(adapt_delta = 0.99))

# TRACEPLOTS, with and without warm-up:
traceplot(fit, pars = my.pars,inc_warmup = TRUE)  # plot traces with burnin
traceplot(fit, pars = my.pars,inc_warmup = FALSE)  # plot traces without burnin

fit0 <- as.matrix(fit, pars = my.pars)
write.table(fit0, 'newfit1.txt', row.names=FALSE)


pairs(fit, pars = my.pars)
# posterior inference analytics
print(fit, digits=2, pars = my.pars)
# MARGINALS and correlations
pairs(fit,pars = my.pars)

