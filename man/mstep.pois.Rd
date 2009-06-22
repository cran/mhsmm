\name{mstep.pois}
\alias{mstep.pois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs re-estimation (the M-step) for a Poisson emission distribution}
\description{
  Re-estimates the parameters of a Poisson emission distribution as part of the EM algorithm for HMMs and HSMMs.
  This is called by the \code{hmm} and \code{hsmm} functions.  It is a suitable prototype function for users
  wishing to design their own emission distributions.
}
\usage{
mstep.pois(x, wt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A vector of observed values}
  \item{wt}{ A T x J matrix of weights.  Column entries are the weights for respective states. }
}
\details{
  Users may write functions that take the same arguments and return the same values for their
  own custom emission distributions.
}
\value{
  Returns the \code{emission} slot of a \code{hmmspec} or \code{hsmmspec} object
  \item{lambda}{Vector of length J containing the Poisson parameters for each state j}
  ...
}

\seealso{ \code{rpois.hsmm}, \code{dpois.hsmm}}

\author{Jared O'Connell}


\examples{
   J<-3
  initial <- rep(1/J,J)
  P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
  b <- list(lambda=c(1,3,6))
  model <- hmmspec(init=initial, trans=P, emission=b,f=dpois.hsmm)
  model
  train <- simulate(model, nsim=300, seed=1234, r=rpois.hsmm)
  plot(train,xlim=c(0,100))  
  h1 = hmmfit(train,model,mstep=mstep.pois)
}