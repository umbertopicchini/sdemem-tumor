#setwd("P:/forskning/tumorgrowth/revision2018")

get.my.data <- function(name){
 dat <- read.table(name)
 names(dat) <- c('time','logvol','id')
 dat <- dat[order(dat$id, dat$time),]
 tab <- table(dat$id)
 ns <- as.numeric(tab)
 is <- as.numeric(names(tab))
 k <- length(tab)

 t <- max(dat$time)
 n <- max(ns)
 imax <- is[ns==n]
 ts <- subset(dat, id==imax[1])$time

 out <- list(name=name, frame=dat, N.sub=k, id=is, n.obs=ns, time=ts)
 out
}

my.cols0 <- rep('black', times=100)
my.cols1=c('blue','red','black','purple','green','orange','turquoise','brown','grey','yellow')
my.cols2 <- rep('grey', times=100)

plot.my.data <- function(dat, cols){
  data <- dat$frame
  is <- dat$id
  k <- dat$N.sub
  ls <- c(min(data$logvol), max(log(1000), max(data$logvol)))
  plot(data$time, data$logvol, type='n', xlim=c(0,1), ylim=ls, xlab='Time', ylab='Log(volume)')
  title(main=dat$name)
  abline(h=log(1000), lty=2, lwd=2, col='grey')
  for (i in 1:k){
    tmp <- subset(data, id==is[i])
    lines(tmp$time, tmp$logvol, lwd=2, col=cols[i])
    points(tmp$time, tmp$logvol, pch=16, col=cols[i])
  }
  0
}
  
#mydat1 <- get.my.data('datagroup1')
#mydat3 <- get.my.data('datagroup3')
#mydat5 <- get.my.data('datagroup5')

#plot.my.data(dat=mydat1, cols=my.cols1)
#plot.my.data(dat=mydat3, cols=my.cols1)
#plot.my.data(dat=mydat5, cols=my.cols1)

