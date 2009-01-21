\name{hsmmspec}
\alias{hsmmspec}
\title{Hidden semi-Markov model specification}
\description{
  Creates a model specification of a hidden semi-Markov model.
}
\usage{
hsmmspec(init, transition, emission, sojourn, r = rnorm.hsmm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{Distribution of states at t=1 ie. P(S=s) at t=1}
  \item{transition}{The transition matrix of the embedded Markov chain (diagonal must be 0)}
  \item{emission}{A list containing the parameters of the emission distribution}
  \item{sojourn}{A list containining the parameters and type of sojourn distribtuion (see Details)}
  \item{r}{The function used to generate observations from the emission distribution}
}
\details{
  This provides a 
}
\value{
  An object of class hsmmspec
}
\references{ Guedon, Y. (2003),  Estimating hidden semi-Markov chains from discrete sequences, 
	Journal of Computational and Graphical Statistics, Volume 12, Number 3, page 604-639 - 2003}
\author{ Jared O'Connell }
\seealso{ hsmm, simulate.hsmmspec, predict.hsmm}
