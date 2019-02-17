# runs posterior predictive chacks for the fitted data from group 3
#
# Warning! Compared to the paper, the summaries corresponding to the slope beta1 of the first-order autoregression
#          are obtained from a model without intercept, i.e. from E(Y_j)=beta1*Y_{j-1}, instead of E(Y_j)=beta0 + beta1*Y_{j-1}.
#          This is because, as documented in the ar.ols function, using ar.ols(...,demean = FALSE, intercept = TRUE)
#          can sometimes complately fail. Thus we use intercept = FALSE.
#          This of course implies a minor discrepancy with summaries reported in the paper for the BSL case.

library(truncnorm)

source('indata.R')
par(mfrow=c(1,1))
mydat3 <- get.my.data('datagroup3')
plot.my.data(dat=mydat3, cols=my.cols1)


fit3 <- read.table('newfit3.txt', header=TRUE)
# BRUTAL HARDCODING! THIS IS THE FIRST MEASURED LOVOLUME FOR EACH SUBJECT
mylogv0 <- c(6.642617,5.109575,4.709530,3.775057,5.655292,5.522660,6.372978,6.622204)


mymad <- function(x){ # the mean absolute deviation frome the mean (not quite the same as R's MAD function)
  sum(abs(x-mean(x)))/length(x)
}

# Compute summary statistics from data:
sumstat.inter <- function(y, n.obs){
  n.sub <- length(n.obs)
  last.ix <- cumsum(n.obs)
  first.ix <- c(1, last.ix[1:(n.sub-1)]+1)
  second.ix <- first.ix+1
  N<-sum(mydat3$n.obs) # total number of observations
  times_equal_last <- mydat3$frame$time==mydat3$frame$time[N]
  out <- c(mymad(y[first.ix]),mymad(y[second.ix]),mymad(y[times_equal_last]))
  names(out) <- c('MAD FIRST','MAD SECOND','MAD LAST')
  out
}

# Observed data summaries
inter.obs3 <- sumstat.inter(y=mydat3$frame$logvol, n.obs=mydat3$n.obs)

sumstat.intra <- function(y, id, time){
  id.tab <- table(id)
  n.obs <- as.numeric(id.tab)
  ids <- names(id.tab)
  n.sub <- length(n.obs)
  last.ix <- cumsum(n.obs)
  first.ix <- c(1, last.ix[1:(n.sub-1)]+1)
  second.ix <- first.ix+1
  
  mads <- tapply(y, id, mymad)
  y1 <- y[first.ix]
  t1 <- time[first.ix]
  y2 <- y[second.ix]
  yn <- y[last.ix]
  tn <- time[last.ix]
  slope <- (yn-y1)/(tn-t1)
  beta <- rep(NA,n.sub)
  for (i in 1:n.sub){
   # beta[i] <- coef(lm(y[id==ids[i]]~time[id==ids[i]]))[2]
    order1autoreg <- ar.ols(y[id==ids[i]], aic = FALSE, order.max = 1, demean = FALSE, intercept = FALSE)
    beta[i] <- order1autoreg$ar[1]
  }
  out <- as.vector(cbind(mads,y1,y2,slope,beta))
  statnames <- c('MAD','Y1','Y2','SLOPE','BETA')
  names(out) <- as.vector(t(outer(statnames, ids, paste)))
  out
  }

# Observed data summaries
intra.obs3 <- sumstat.intra(y=mydat3$frame$logvol, id=mydat3$frame$id, time=mydat3$frame$time)

# Simulate summary statistics from posterior predictive distribution:
sim.pred.check <- function(n.rep, parms, id, time, initv0, my.seed){
  set.seed(my.seed)
  id.tab <- table(id)
  n.obs <- as.numeric(id.tab)
  n.sub <- length(n.obs)
  n.total <- sum(n.obs)
  t <- c(0,cumsum(n.obs))
  n.parms <- dim(parms)[1]
  x <- time-min(time)

  out.intra <- matrix(NA,n.rep,5*n.sub)
  statnames <- c('MAD','Y1','Y2','SLOPE','BETA')
  colnames(out.intra) <- as.vector(t(outer(statnames, names(id.tab), paste)))
  
  out.inter <- matrix(NA,n.rep,3)
  colnames(out.inter) <- c('MAD.FIRST','MAD.SECOND','MAD.LAST')
    
  for (k in 1:n.rep){
   print(k)
   parm <- as.numeric(parms[sample(1:n.parms, 1),])
   names(parm) <- names(parms)

   logv0 <- initv0
#   logv0 <- rnorm(n.sub, mean=parm['mu_v0'], sd=parm['sigma_v0'])
   b <- rnorm(n.sub, mean=parm['mu_beta'], sd=parm['sigma_beta'])
   d <- rnorm(n.sub, mean=parm['mu_delta'], sd=parm['sigma_delta'])
   a <- rep(NA, n.sub)
   for (i in 1:n.sub){a[i] <- rtruncnorm(1, a=0, b=1, mean=parm['mu_alpha'], sd=parm['sigma_alpha'])}
   sigma <- parm['sigma_epsilon']

   y <- rep(NA, n.total)
   for (n in 1:n.sub){
    for (j in (t[n]+1):t[n+1]){
      y[j] <-rnorm(1,mean=logv0[n]+log((1-a[n])*exp(b[n]*x[j])+a[n]*exp(-1*d[n]*x[j])), sd=sigma)
    }
  }
   out.inter[k,] <- sumstat.inter(y, n.obs)
   out.intra[k,] <-sumstat.intra(y, id, time)
  }

   out <- cbind(out.inter, out.intra)
   as.data.frame(out)
}


# 10000 simulated summaries:
sim.pred3 <- sim.pred.check(n.rep=10000, parms=fit3, id=mydat3$frame$id, time=mydat3$frame$time, initv0=mylogv0, my.seed=345)
# To be compared with
obs.pred3 <- c(inter.obs3, intra.obs3)

write.table(sim.pred3, 'newpostpred3.txt', row.names=FALSE)
