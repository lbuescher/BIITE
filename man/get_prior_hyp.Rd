% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/immunogen_predictive_prior.r
\name{get_prior_hyp}
\alias{get_prior_hyp}
\title{P(H)}
\usage{
get_prior_hyp(molecs, h, shape.df = NULL)
}
\arguments{
\item{molecs}{Vector with names of the molecules in the cohort.}

\item{h}{Hypothesis we want the prior of}

\item{shape.df}{Dataframe with shape parameters for the prior Beta distributions. If NULL, a uniform, non-informative prior is used.}
}
\description{
Computes the log of the prior of a given hypothesis. Returns 0 when using a uniform prior.
}

