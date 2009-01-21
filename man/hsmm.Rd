\name{hsmm}
\alias{hsmm}
\title{hidden semi-Markov model}
\description{
  Estimates parameters of a HSMM using the EM algorithm.
}
\usage{
hsmm(x, model, f, mstep, maxit = 100, sojourn.distribution = c("nonparametric","gamma","poisson")
      ,lock.transition = FALSE, lock.d = FALSE, M=NA,graphical=FALSE)
}
\arguments{
  \item{x}{A hsmm.data object (see Details)}
  \item{model}{ Starting parameters for the model (see Details) }
  \item{f}{ Density function of the emmission distribution }
  \item{mstep}{ Re-estimates the parameters of density function on each iteration }
  \item{maxit}{ Maximum number of iterations }
  \item{sojourn.distribution}{ Waiting or sojourn distribution used (see Details) }
  \item{lock.transition}{ If TRUE will not re-estimate the transition matrix }
  \item{lock.d}{ If TRUE will not re-estimate the sojourn time density }
  \item{M}{Maximum number of time spent in a state (truncates the waiting distribution)}
  \item{graphical}{If TRUE will plot the sojourn densities on each iteration}  
}
\value{
  \item{start}{A vector of the starting probabilities for each state}
  \item{a}{The transition matrix of the embedded Markov chain}
  \item{emission}{A list of the parameters of the emission distribution}
  \item{waiting}{A list of the parameters of the waiting distribution}	
}
\references{ Guedon, Y. (2003),  Estimating hidden semi-Markov chains from discrete sequences, 
	Journal of Computational and Graphical Statistics, Volume 12, Number 3, page 604-639 - 2003 }
\author{ Jared O'Connell }
\seealso{ hsmmspec, simulate.hsmmspec, predict.hsmm}
\examples{
J <- 3
init <- c(0,0,1)
P <- matrix(c(0,.1,.4,.5,0,.6,.5,.9,0),nrow=J)
B <- list(mu=c(10,15,20),sigma=c(2,1,1.5))
d <- list(lambda=c(10,30,60),shift=c(10,100,30),type='poisson')
model <- hsmmspec(init,P,emission=B,sojourn=d,r=rnorm.hsmm)          
train <- simulate(model,nsim=100,seed=123456)
plot(train,xlim=c(0,400))
start.poisson <- hsmmspec(init=rep(1/J,J),transition=matrix(c(0,.5,.5,.5,0,.5,.5,.5,0),nrow=J),emission=list(mu=c(4,12,23),
		sigma=c(1,1,1)),sojourn=list(lambda=c(9,25,40),shift=c(5,95,45)))

M=500
h.poisson <- hsmm(train,start.poisson,f=dnorm.hsmm,mstep=mstep.norm,sojourn.distribution='poisson',M=M)
plot(h.poisson$loglik,type='b',ylab='Log-likelihood',xlab='Iteration') ##has it converged?
summary(h.poisson)
predicted <- predict(h.poisson,train)
table(train$s,predicted$s) ##classification matrix
mean(predicted$s!=train$s) ##misclassification rate

d <- cbind(dunif(1:M,0,50),dunif(1:M,100,175),dunif(1:M,50,130))
start.np <- hsmmspec(init=rep(1/J,J),transition=matrix(c(0,.5,.5,.5,0,.5,.5,.5,0),nrow=J),emission=list(mu=c(4,12,23),
		sigma=c(1,1,1)),sojourn=list(d=d))
h.np <- hsmm(train,start.np,f=dnorm.hsmm,mstep=mstep.norm,sojourn.distribution='nonparametric',M=M,graphical=TRUE)
mean(predicted$s!=train$s) ##misclassification rate
}

