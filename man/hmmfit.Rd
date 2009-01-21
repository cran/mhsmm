\name{hmmfit}
\alias{hmmfit}
\title{fit a hidden Markov model}
\description{
  Estimates parameters of a HMM using the EM algorithm.
}
\usage{
hmmfit(x,start.val,f=dnorm.hsmm,mstep=mstep.norm,tol=1e-08,maxit=1000) 
}
\arguments{
  \item{x}{A hsmm.data object (see Details)}
  \item{start.val}{ Starting parameters for the model (see Details) }
  \item{f}{ Density function of the emmission distribution }
  \item{mstep}{ Re-estimates the parameters of density function on each iteration }
  \item{maxit}{ Maximum number of iterations }
  \item{tol}{ Convergence tolerance }
}
\value{
  \item{start}{A vector of the starting probabilities for each state}
  \item{a}{The transition matrix of the embedded Markov chain}
  \item{emission}{A list of the parameters of the emission distribution}
}
\references{Rabiner, L. (1989), A tutorial on hidden Markov models and selected applications
in speech recognition, Proceedings of the IEEE, 77, 257-286.}
\author{ Jared O'Connell }
\seealso{ sim.hmm, sim.mhmm, predict.hmm}
\examples{
J<-3
initial <- rep(1/J,J)
P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
b <- list(mu=c(-3,0,2),sigma=c(2,1,.5))
model <- hmmspec(init=initial, trans=P, emission=b, r=rnorm.hsmm)
model

train <- simulate(model, nsim=300, seed=1234)
plot(train,xlim=c(0,100))

init0 <- rep(1/J,J)
P0 <- matrix(1/J,nrow=J,ncol=J)
b0 <- list(mu=c(-3,1,3),sigma=c(1,1,1))
startval <- hmmspec(init=init0, trans=P0,emission=b0) 
h1 = hmmfit(train,startval,f=dnorm.hsmm,mstep=mstep.norm)

plot(h1$loglik,type='b',ylab='Log-likelihood',xlab='Iteration')
summary(h1)

#proportion of incorrect states
mean(train$s!=predict(h1,train)$s)

#simulate a new test set 
test <- simulate(model, nsim=c(100,200,300), seed=1234)
mean(test$s!=predict(h1,test)$s)

}

