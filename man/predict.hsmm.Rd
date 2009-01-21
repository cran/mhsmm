\name{predict.hsmm}
\alias{predict.hsmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for hsmms}
\description{
  Predicts the underlying state sequence for an observed sequence \code{x} given a \code{hsmm} model
}
\usage{
\method{predict}{hsmm}(object, x, method = "viterbi", ...)
}

\arguments{
  \item{object}{ An object of type \code{hsmm}}
  \item{x}{ A vector or dataframe of observations }
  \item{method}{ Only "viterbi" is available at the moment}
  \item{\dots}{further arguments passed to or from other methods.}

}
\details{
    This technique applies the Viterbi algorithm for HSMMs, producing the most likely sequence of states given the observed data.
}
\value{
  Returns a \code{hsmm.data} object, suitable for plotting.
  \item{x }{A vector or data.frame of observations}
  \item{s }{A vector containing the reconstructed state sequence}
  \item{N }{The lengths of each sequence}
}
\references{ Guedon, Y. (2003),  Estimating hidden semi-Markov chains from discrete sequences, 
	Journal of Computational and Graphical Statistics, Volume 12, Number 3, page 604-639 - 2003}
\author{Jared O'Connell}

\seealso{ hsmm }
\examples{
##See 'hsmm' for examples
}