% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rfuncs.R
\name{TSIRlyap}
\alias{TSIRlyap}
\title{Function to do  Lyapunov exponent calculations from a TSIR simulation}
\usage{
TSIRlyap(I, S, alpha, bt, N)
}
\arguments{
\item{I}{a vector containing the time series of Is}

\item{S}{vector containing the time series of Ss}

\item{alpha}{the exponent on I}

\item{bt}{the seasonal transmission coefficients}

\item{N}{the population size}
}
\value{
An object of class lyap with the lyapunov exponent, values for the Jacobians, parameters and data
}
\description{
Function to do  Lyapunov exponent calculations from a TSIR simulation
}
\examples{
\dontrun{see chapter 10 in book}
}
