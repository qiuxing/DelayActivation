\name{ModesDetect}
\alias{ModesDetect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Computes the number of modes (local minima/maxima) of each gene.
}

\description{ \code{ModesDetect} ...  }

\usage{
ModesDetect(x.smoothed, delta=0.1)
}
%- maybe also 'usage' for other objects documented here.

\arguments{

  \item{x.smoothed}{An m by T dimensional matrix of \emph{smoothed}
  data, each row represents a gene and each column a time point.}
  \item{delta}{Error tolerance.  Default to 0.1.}
}

\details{ \option{x.smoothed} is the fitted (by a nonparametric
  regression method) values of the original data.  Modes are defined as
  the \dQuote{V} shape or the \dQuote{lambda} shape points. See the
  attached report for more information.
 }

\value{A vector of estimated number of modes for each gene.
}

\references{some references.}

\author{
  Xing Qiu, Hulin Wu, Shuang Wu.
}

\seealso{\code{\link{DelayDetect}}, \code{\link{ActivationDetect}},
  \code{\link{FirstRun}}, \code{\link{PostProc}}. }

\examples{
library(DelayActivation)
data(lung2); Time <- 0:10

## smoothing
mybasis <- create.bspline.basis(range(Time), length(Time)+4-2, 4, Time)
smoothing.par <- fdPar(mybasis, Lfdobj=2, lambda=5.0)
fitted.curves <- smooth.basis(Time, t(lung2), smoothing.par)[["fd"]]
lung2.smooth <- t(eval.fd(Time, fitted.curves))

mymodes <- ModesDetect(lung2.smooth)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{clustering}
\keyword{misc}
