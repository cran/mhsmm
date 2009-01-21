\name{predict.hmm}
\alias{predict.hmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prediction function for hmm }
\description{
  Predicts the underlying state sequence for an observed sequence \code{x} given a \code{hmm} model
}
\usage{
\method{predict}{hmm}(object, x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{hmm}}
  \item{x}{ A vector or data.frame of observations }
  \item{\dots}{further arguments passed to or from other methods.}

}
\details{
  This technique applies the Viterbi algorithm for HMMs, producing the most likely sequence of states given the observed data.
}
\value{
  Returns a \code{hsmm.data} object, suitable for plotting.
  \item{x }{A vector or data.frame of observations}
  \item{s }{A vector containing the reconstructed state sequence}
  \item{N }{The lengths of each sequence}
}
\references{Rabiner, L. (1989), A tutorial on hidden Markov models and selected applications
in speech recognition, Proceedings of the IEEE, 77, 257-286.}

\author{Jared O'Connell}
\seealso{\code{hmm}}
\examples{
##See examples in 'hmm'
}