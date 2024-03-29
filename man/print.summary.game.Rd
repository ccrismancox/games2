% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/games.r
\name{print.summary.game}
\alias{print.summary.game}
\title{Print strategic model summary}
\usage{
\method{print}{summary.game}(x, ...)
}
\arguments{
\item{x}{an object of class \code{summary.game}, typically produced by
running \code{summary} on a fitted model of class \code{game}}

\item{...}{other arguments, currently ignored}
}
\value{
\code{x}, invisibly.
}
\description{
Print output from \code{summary.game}
}
\details{
Prints the standard regression results table from a fitted strategic model,
along with the log-likelihood, AIC, and number of observations.
}
\seealso{
\code{\link{summary.game}}
}
\author{
Brenton Kenkel (\email{brenton.kenkel@gmail.com})
}
