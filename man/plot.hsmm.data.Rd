\name{plot.hsmm.data}
\alias{plot.hsmm.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot function for hsmm data }
\description{
  Produces a plot of the observed sequences, and displays a coloured bar signifying the hidden states (if available)
}
\usage{
\method{plot}{hsmm.data}(x, ...)
}

\arguments{
  \item{x}{ A \code{hsmm.data} object}
  \item{\dots}{ Arguments passed to \code{plot.ts} }
}

\seealso{ \code{addStates}}
\examples{
  J<-3
  initial <- rep(1/J,J)
  P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
  b <- list(mu=c(-3,0,2),sigma=c(2,1,.5))
  model <- hmmspec(init=initial, trans=P, emission=b,f=dnorm.hsmm)
  
  train <- simulate(model, nsim=300, seed=1234, r=rnorm.hsmm)
  plot(train,xlim=c(0,100))
}
\author{Jared O'Connell}