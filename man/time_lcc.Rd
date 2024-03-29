% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time_lcc.R
\name{time_lcc}
\alias{time_lcc}
\title{Regular Sequence Generator for Time Variable}
\usage{
time_lcc(time, from, to, n)
}
\arguments{
\item{time}{A numeric vector of unique time values.}

\item{from}{The starting (minimum) value for the time sequence.}

\item{to}{The ending (maximum) value for the time sequence.}

\item{n}{Desired length of the sequence (integer). Typically, a value
between 30 and 50 is adequate.}
}
\value{
A numeric vector containing a regular sequence of time values,
  including the unique values from the input time vector.
}
\description{
Generates a regular sequence for the time variable, including
  the unique values from the input time vector. This function is used
  internally to construct LCC, LPC, and LA curves and their simultaneous
  confidence intervals.
}
\examples{
data(hue)
attach(hue)
time_lcc(time = Time, from = min(Time), to = max(Time), n = 30)
detach(hue)

}
\keyword{internal}
