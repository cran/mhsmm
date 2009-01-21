\name{sim.mc}
\alias{sim.mc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Markov chain simulation }
\description{
  Simulates a Markov chain
}
\usage{
sim.mc(init, transition, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{The distribution of states at the first time step}
  \item{transition}{The transition probability matrix of the Markov chain}
  \item{N}{The number of observations to simulate}
}
\value{
  A vector of integers representing the state sequence.
}
\author{Jared O'Connell}

