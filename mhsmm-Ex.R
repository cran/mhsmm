pkgname <- "mhsmm"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('mhsmm')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("addStates")
### * addStates

flush(stderr()); flush(stdout())

### Name: addStates
### Title: Adds a bar representing state sequence.
### Aliases: addStates

### ** Examples

  plot(rnorm(100),type='l')
  addStates(rep(c(1,2),each=50))  

  plot(seq(0.01,1,.01),rnorm(100),type='l')
  addStates(rep(c(1,2),each=50),seq(0.01,1,.01))  



cleanEx()
nameEx("dmvnorm.hsmm")
### * dmvnorm.hsmm

flush(stderr()); flush(stdout())

### Name: dmvnorm.hsmm
### Title: Emission ensity function for a multivariate normal emission
###   distribution
### Aliases: dmvnorm.hsmm

### ** Examples

  J<-2
  initial <- rep(1/J,J)
  P <- matrix(c(.3,.5,.7,.5),nrow=J)
  b <- list(mu=list(c(-3,0),c(1,2)),sigma=list(diag(2),matrix(c(4,2,2,3), ncol=2)))
  model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dmvnorm.hsmm)
  model
  train <- simulate(model, nsim=300, seed=1234, rand.emis=rmvnorm.hsmm)
  plot(train,xlim=c(0,100))
  h1 = hmmfit(train,model,mstep=mstep.mvnorm)



cleanEx()
nameEx("dpois.hsmm")
### * dpois.hsmm

flush(stderr()); flush(stdout())

### Name: dpois.hsmm
### Title: Emission density function for Poisson emission distribution
### Aliases: dpois.hsmm

### ** Examples

  J<-3
  initial <- rep(1/J,J)
  P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
  b <- list(lambda=c(1,3,6))
  model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dpois.hsmm)
  model
  train <- simulate(model, nsim=300, seed=1234, rand.emis=rpois.hsmm)
  plot(train,xlim=c(0,100))  
  h1 = hmmfit(train,model,mstep=mstep.pois)



cleanEx()
nameEx("gammafit")
### * gammafit

flush(stderr()); flush(stdout())

### Name: gammafit
### Title: Parameter estimation for the Gamma distribution
### Aliases: gammafit

### ** Examples

  gammafit(rgamma(1000,shape=10,scale=13))



cleanEx()
nameEx("hmmfit")
### * hmmfit

flush(stderr()); flush(stdout())

### Name: hmmfit
### Title: fit a hidden Markov model
### Aliases: hmmfit

### ** Examples

J<-3
initial <- rep(1/J,J)
P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
b <- list(mu=c(-3,0,2),sigma=c(2,1,.5))
model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dnorm.hsmm)
model

train <- simulate(model, nsim=300, seed=1234, rand.emis=rnorm.hsmm)
plot(train,xlim=c(0,100))

init0 <- rep(1/J,J)
P0 <- matrix(1/J,nrow=J,ncol=J)
b0 <- list(mu=c(-3,1,3),sigma=c(1,1,1))
startval <- hmmspec(init=init0, trans=P0,parms.emission=b0,dens.emission=dnorm.hsmm) 
h1 = hmmfit(train,startval,mstep=mstep.norm)

plot(h1$loglik,type='b',ylab='Log-likelihood',xlab='Iteration')
summary(h1)

#proportion of incorrect states
mean(train$s!=predict(h1,train)$s)

#simulate a new test set 
test <- simulate(model, nsim=c(100,200,300), seed=1234,rand.emis=rnorm.hsmm)
mean(test$s!=predict(h1,test)$s)



cleanEx()
nameEx("hsmmfit")
### * hsmmfit

flush(stderr()); flush(stdout())

### Name: hsmmfit
### Title: fit a hidden semi-Markov model
### Aliases: hsmmfit

### ** Examples

J <- 3
init <- c(0,0,1)
P <- matrix(c(0,.1,.4,.5,0,.6,.5,.9,0),nrow=J)
B <- list(mu=c(10,15,20),sigma=c(2,1,1.5))
d <- list(lambda=c(10,30,60),shift=c(10,100,30),type='poisson')
model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)
train <- simulate(model,r=rnorm.hsmm,nsim=100,seed=123456)
plot(train,xlim=c(0,400))
start.poisson <- hsmmspec(init=rep(1/J,J),
  transition=matrix(c(0,.5,.5,.5,0,.5,.5,.5,0),nrow=J),
  parms.emission=list(mu=c(4,12,23),
		sigma=c(1,1,1)),
  sojourn=list(lambda=c(9,25,40),shift=c(5,95,45),type='poisson'),
 dens.emission=dnorm.hsmm)

M=500
# not run (takes some time)
# h.poisson <- hsmmfit(train,start.poisson,mstep=mstep.norm,M=M)
# plot(h.poisson$loglik,type='b',ylab='Log-likelihood',xlab='Iteration') ##has it converged?
# summary(h.poisson)
# predicted <- predict(h.poisson,train)  
# table(train$s,predicted$s) ##classification matrix
# mean(predicted$s!=train$s) ##misclassification rate

d <- cbind(dunif(1:M,0,50),dunif(1:M,100,175),dunif(1:M,50,130))
start.np <- hsmmspec(init=rep(1/J,J),
  transition=matrix(c(0,.5,.5,.5,0,.5,.5,.5,0),nrow=J),
  parms.emission=list(mu=c(4,12,23),
  sigma=c(1,1,1)),
  sojourn=list(d=d,type='nonparametric'),
  dens.emission=dnorm.hsmm)
# not run (takes some time)
# h.np <- hsmmfit(train,start.np,mstep=mstep.norm,M=M,graphical=TRUE)
# mean(predicted$s!=train$s) ##misclassification rate



cleanEx()
nameEx("mstep.mvnorm")
### * mstep.mvnorm

flush(stderr()); flush(stdout())

### Name: mstep.mvnorm
### Title: Performs re-estimation (the M-step) for a multivariate normal
###   emission distribution
### Aliases: mstep.mvnorm

### ** Examples

  J<-2
  initial <- rep(1/J,J)
  P <- matrix(c(.3,.5,.7,.5),nrow=J)
  b <- list(mu=list(c(-3,0),c(1,2)),sigma=list(diag(2),matrix(c(4,2,2,3), ncol=2)))
  model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dmvnorm.hsmm)
  model
  train <- simulate(model, nsim=300, seed=1234, rand.emis=rmvnorm.hsmm)
  plot(train,xlim=c(0,100))
  h1 = hmmfit(train,model,mstep=mstep.mvnorm)



cleanEx()
nameEx("mstep.pois")
### * mstep.pois

flush(stderr()); flush(stdout())

### Name: mstep.pois
### Title: Performs re-estimation (the M-step) for a Poisson emission
###   distribution
### Aliases: mstep.pois

### ** Examples

  J<-3
  initial <- rep(1/J,J)
  P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
  b <- list(lambda=c(1,3,6))
  model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dpois.hsmm)
  model
  train <- simulate(model, nsim=300, seed=1234, rand.emis=rpois.hsmm)
  plot(train,xlim=c(0,100))  
  h1 = hmmfit(train,model,mstep=mstep.pois)



cleanEx()
nameEx("plot.hsmm.data")
### * plot.hsmm.data

flush(stderr()); flush(stdout())

### Name: plot.hsmm.data
### Title: Plot function for hsmm data
### Aliases: plot.hsmm.data

### ** Examples

  J<-3
  initial <- rep(1/J,J)
  P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
  b <- list(mu=c(-3,0,2),sigma=c(2,1,.5))
  model <- hmmspec(init=initial, trans=P, parms.emission=b, dens.emission=dnorm.hsmm)
  
  train <- simulate(model, nsim=300, seed=1234, rand.emis=rnorm.hsmm)
  plot(train,xlim=c(0,100))



cleanEx()
nameEx("predict.hmm")
### * predict.hmm

flush(stderr()); flush(stdout())

### Name: predict.hmm
### Title: Prediction function for hmm
### Aliases: predict.hmm

### ** Examples

##See examples in 'hmmfit'



cleanEx()
nameEx("predict.hmmspec")
### * predict.hmmspec

flush(stderr()); flush(stdout())

### Name: predict.hmmspec
### Title: Prediction function for hmmspec
### Aliases: predict.hmmspec

### ** Examples

     J<-3
     initial <- rep(1/J,J)
     P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
     b <- list(mu=c(-3,0,2),sigma=c(2,1,.5))
     model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dnorm.hsmm)    
     train <- simulate(model, nsim=300, seed=1234, rand.emis=rnorm.hsmm)
     mean(predict(model,train)$s!=train$s) #error rate when true model is known



cleanEx()
nameEx("predict.hsmm")
### * predict.hsmm

flush(stderr()); flush(stdout())

### Name: predict.hsmm
### Title: Prediction for hsmms
### Aliases: predict.hsmm

### ** Examples

##See 'hsmmfit' for examples



cleanEx()
nameEx("predict.hsmmspec")
### * predict.hsmmspec

flush(stderr()); flush(stdout())

### Name: predict.hsmmspec
### Title: Prediction for hsmmspec
### Aliases: predict.hsmmspec

### ** Examples

J <- 3
init <- c(0,0,1)
P <- matrix(c(0,.1,.4,.5,0,.6,.5,.9,0),nrow=J)
B <- list(mu=c(10,15,20),sigma=c(2,1,1.5))
d <- list(lambda=c(10,30,60),shift=c(10,100,30),type='poisson')
model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)
train <- simulate(model,r=rnorm.hsmm,nsim=100,seed=123456)
mean(predict(model,train,M=500)$s!=train$s) #error rate when true model is known



cleanEx()
nameEx("reprocows")
### * reprocows

flush(stderr()); flush(stdout())

### Name: reprocows
### Title: Reproductive data from seven dairy cows
### Aliases: reprocows
### Keywords: datasets

### ** Examples

data(reprocows)
data(reproai)
data(reproppa)
tm = 1600

J <- 3
init <- c(1,0,0)
trans <- matrix(c(0,0,0,1,0,1,0,1,0),nrow=J)
emis <- list(mu=c(0,2.5,0),sigma=c(1,1,1))

N <- as.numeric(table(reprocows$id))
train <- list(x=reprocows$activity,N=N)
class(train) <- "hsmm.data"
tmp <- gammafit(reproppa * 24)
M <- max(N)

d <- cbind(dgamma(1:M,shape=tmp$shape,scale=tmp$scale),
 # ppa sojourn directly estimated from ppa data set
dunif(1:M,4,30),
 # oestrus between 4 and 30 hours
dunif(1:M,15*24,40*24))
 #not-oestrus between 15 and 40 days

startval <- hsmmspec(init,trans,parms.emission=emis,list(d=d,type='gamma'),
  dens.emission=dnorm.hsmm)
#not run (takes some time)
#h.activity <- hsmmfit(train,startval,mstep=mstep.norm,maxit=10,M=M,lock.transition=TRUE)
  



cleanEx()
nameEx("rmvnorm.hsmm")
### * rmvnorm.hsmm

flush(stderr()); flush(stdout())

### Name: rmvnorm.hsmm
### Title: Random number generation from a multivariate normal distributed
###   emission distribution
### Aliases: rmvnorm.hsmm

### ** Examples

  J<-2
  initial <- rep(1/J,J)
  P <- matrix(c(.3,.5,.7,.5),nrow=J)
  b <- list(mu=list(c(-3,0),c(1,2)),sigma=list(diag(2),matrix(c(4,2,2,3), ncol=2)))
  model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dmvnorm.hsmm)
  model
  train <- simulate(model, nsim=300, seed=1234, rand.emis=rmvnorm.hsmm)
  plot(train,xlim=c(0,100))
  h1 = hmmfit(train,model,mstep=mstep.mvnorm)



cleanEx()
nameEx("rpois.hsmm")
### * rpois.hsmm

flush(stderr()); flush(stdout())

### Name: rpois.hsmm
### Title: Random number generation from a Poisson distributed emission
###   distribution
### Aliases: rpois.hsmm

### ** Examples

  J<-3
  initial <- rep(1/J,J)
  P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
  b <- list(lambda=c(1,3,6))
  model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dpois.hsmm)
  model
  train <- simulate(model, nsim=300, seed=1234, rand.emis=rpois.hsmm)
  plot(train,xlim=c(0,100))  
  h1 = hmmfit(train,model,mstep=mstep.pois)



cleanEx()
nameEx("sim.mc")
### * sim.mc

flush(stderr()); flush(stdout())

### Name: sim.mc
### Title: Markov chain simulation
### Aliases: sim.mc

### ** Examples

  p <- matrix(c(.1,.3,.6,rep(1/3,3),0,.5,.5),ncol=3,byrow=TRUE)
  init <- rep(1/3,3)
  sim.mc(init,p,10)  
  
  



cleanEx()
nameEx("simulate.hmmspec")
### * simulate.hmmspec

flush(stderr()); flush(stdout())

### Name: simulate.hmmspec
### Title: Simulation of hidden Markov models
### Aliases: simulate.hmmspec

### ** Examples


J<-3
initial <- rep(1/J,J)
P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
b <- list(mu=c(-3,0,2),sigma=c(2,1,.5))
model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dnorm.hsmm)
train <- simulate(model, nsim=100, seed=1234, rand.emis=rnorm.hsmm)
plot(train)




cleanEx()
nameEx("simulate.hsmmspec")
### * simulate.hsmmspec

flush(stderr()); flush(stdout())

### Name: simulate.hsmmspec
### Title: Simulation for HSMMs
### Aliases: simulate.hsmmspec

### ** Examples

J <- 3
init <- c(0,0,1)
P <- matrix(c(0,.1,.4,.5,0,.6,.5,.9,0),nrow=J)
B <- list(mu=c(10,15,20),sigma=c(2,1,1.5))
d <- list(lambda=c(10,30,60),shift=c(10,100,30),type='poisson')
model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)
train <- simulate(model,rand.emis=rnorm.hsmm,nsim=100,seed=123456)
plot(train,xlim=c(0,400))



cleanEx()
nameEx("smooth.discrete")
### * smooth.discrete

flush(stderr()); flush(stdout())

### Name: smooth.discrete
### Title: Smoothing a discrete time series.
### Aliases: smooth.discrete print.smoothDiscrete summary.smoothDiscrete
###   predict.smoothDiscrete createTransition
### Keywords: models

### ** Examples

## Please see the vignette



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
