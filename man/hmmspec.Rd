\name{hmmspec}
\alias{hmmspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Specificatin of HMMs }
\description{
  Creates a model specficiation for a hidden Markov model
}
\usage{
hmmspec(init, trans, emission, r = rnorm.hsmm, f = dnorm.hsmm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{init}{Distribution of states at t=1 ie. P(S=s) at t=1}
  \item{trans}{The transition matrix of the Markov chain}
  \item{emission}{A list containing the parameters of the emission distribution}
  \item{r}{The function used to generate observations from the emission distribution}
  \item{f}{Density function of the emission distribution.}
}
\value{
  A hmmspec object
}
\references{Rabiner, L. (1989), A tutorial on hidden Markov models and selected applications
in speech recognition, Proceedings of the IEEE, 77, 257-286.}
\author{ Jared O'Connell }

\seealso{ simulate.hmmspec, simulate.hmm, hmm, predict.hmm}
