% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods.R
\name{ci.Dist}
\alias{ci.Dist}
\title{Confidence Intervals for Bias-corrected Dist Estimators}
\usage{
\method{ci}{Dist}(
  object,
  probability = FALSE,
  alpha = 0.05,
  alternative = c("two.sided", "less", "greater"),
  ...
)
}
\arguments{
\item{object}{An object of class `Dist`, a result of a call to `Dist`}

\item{alpha}{Level of significance to construct confidence interval
(default=0.05)}

\item{alternative}{Indicates the alternative hypothesis to construct
confidence interval and must be one of "two.sided" (default), "less", or
"greater".}

\item{...}{arguments to pass down}
}
\value{
A vector giving lower and upper confidence limits for bias-corrected
  estimator
}
\description{
Computes confidence intervals for bias-corrected estimator
}
\examples{
\dontrun{
## -- Continuing the Dist(.) example:
out <- ci(Est)
out
}
}
\keyword{internal}
