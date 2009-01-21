\name{simulate.hsmmspec}
\alias{simulate.hsmmspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulation for HSMMs }
\description{
  Simulates values for a specified hidden semi-Markov model
}
\usage{
\method{simulate}{hsmmspec}(object, nsim, seed = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hsmmspec} object }
  \item{nsim}{An integer or vector of integers (for multiple sequences) specifying the length of the sequence(s) }
  \item{seed}{ \code{seed} for the random number generator}
  \item{\dots}{further arguments passed to or from other methods.}

}
\details{
  If \code{nsim} is a single integer then a HSMM of that length is produced.  If \code{nsim} is a vector of integers, then \code{length(nsim)} sequences are generated with respective lengths.
}
\value{
  An object of class hmmdata
  \item{x}{A vector of length \code{sum(N)} - the sequence(s) of observed values}
  \item{s}{A vector of length \code{sum(N)} - the sequence(s) of hidden states}
  \item{N}{A vector of the length of each observation sequence (used to segment x and s)}
}

\references{ Guedon, Y. (2003),  Estimating hidden semi-Markov chains from discrete sequences, 
	Journal of Computational and Graphical Statistics, Volume 12, Number 3, page 604-639 - 2003 }
\author{ Jared O'Connell }


\seealso{ hsmm, simulate.hsmmspec, predict.hsmm}
