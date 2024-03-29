% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predProbs.r
\name{plot.predProbs}
\alias{plot.predProbs}
\title{Plot predicted probabilities}
\usage{
\method{plot}{predProbs}(x, which = NULL, ask = FALSE, ...)
}
\arguments{
\item{x}{an object of class \code{predProbs} (i.e., a data frame returned by
\code{\link{predProbs}}).}

\item{which}{optional integer specifying which plot (as numbered in the menu
displayed when \code{ask == TRUE}) to make.  If none is given, all available
plots are printed in succession.}

\item{ask}{logical: display interactive menu with options for which plot to
make?}

\item{...}{further arguments to pass to the plotting function.  See
\code{\link{plot.default}} (when the variable on the x-axis is continuous)
or \code{\link{bxp}} (when it is discrete).}
}
\value{
an object of class \code{preplot.predProbs}, invisibly.  This contains
the raw information used by lower-level plotting functions.
}
\description{
Plots predicted probabilities and associated confidence bands, using the
data returned from a call to \code{\link{predProbs}}.
}
\details{
Most \code{predProbs} objects will be associated with multiple plots: one for
each outcome in the estimated model.  These are the three or four terminal
nodes for a \code{\link{egame12}} or \code{\link{egame122}} model
respectively; for an \code{\link{ultimatum}} model, these are the expected
offer and the probability of acceptance.  By default, \code{plot.predProbs}
produces plots for all of them, so only the last will be visible unless the
graphics device is set to have multiple figures (e.g., by setting
\code{par(mfrow = ...)}).  The argument \code{ask} displays a menu to select
among the possible plots for a given object, and \code{which} allows for
this to be done non-interactively.
}
\examples{
data("war1800")
f1 <- esc + war ~ s_wt_re1 + revis1 | 0 | regime1 | balanc + regime2
m1 <- egame12(f1, data = war1800, boot = 10)
pp1 <- predProbs(m1, x = "balanc", n = 5)
pp2 <- predProbs(m1, x = "regime1")

## if "ask" is FALSE and "which" isn't specified, all plots are printed
op <- par(mfrow = c(2, 2))
plot(pp1)
par(op)

\dontrun{
    plot(pp1, ask = TRUE)

    ## Make a plot selection (or 0 to exit):
    ##   1: plot: Pr(~esc)
    ##   2: plot: Pr(esc,~war)
    ##   3: plot: Pr(esc,war)
    ##   4: plot all terms
}

## To change line type for confidence bounds, use argument `lty.ci`
plot(pp1, which = 3, lty.ci = 3)

## All the standard plotting options work too
plot(pp1, which = 3, xlab = "Capabilities", ylab = "Probability", main = "Title")

## Discrete `x` variables are plotted via R's boxplot functionality
plot(pp2, which = 3)
}
\author{
Brenton Kenkel (\email{brenton.kenkel@gmail.com})
}
