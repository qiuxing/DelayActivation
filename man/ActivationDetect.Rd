\name{ActivationDetect}
\alias{ActivationDetect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A function to determined the up/down-regulation of each gene, and
  their activation time.
}

\description{ \code{ActivationDetect} ...  }

\usage{
ActivationDetect(x.smoothed, Time, Q=0.5)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x.smoothed}{An m by T dimensional matrix of \emph{smoothed}
  data, each row represents a gene and each column a time point.}
  \item{Time}{Time points. It must be a sequence of increasing numbers
    with length equal to the number of columns of x.}
  \item{Q}{Quantile cutoff used in detecting the activation time.}
}

\details{ \option{x.smoothed} is the fitted (by a nonparametric
  regression method) values of the original data.  We first compute the
  global minimum (\option{xmin}) and maximum (\option{xmax}) of
  \option{x.smooth} and determine if this gene is up or down regulated
  by this rule: if \code{xmax - x.smoothed[1] >= x.smoothed[1] - xmin}
  then we call this gene up-regulated; otherwise down-regulated. 

  Activation are defined as the first time point at which
  \code{x.smooth[i] - x.smooth[1] > Q * (xmax - x.smooth[1])} (for
  up-regulated genes) or \code{x.smooth[i] - x.smooth[1] < Q * (xmin -
  x.smooth[1])} (for down-regulated genes).
}

\value{A data frame which contains the following three columns
  \item{UpDown}{Up or down regulation of genes.}
  \item{Activation.Day}{Activation day (time).}
  \item{Activation.Index}{The index of the Activation day time. In other
  words, \code{Time[Activation.Index[i]] == Activation.Day[i]}.}
}

\references{some references.}

\author{
  Xing Qiu, Hulin Wu, Shuang Wu.
}

\seealso{\code{\link{DelayDetect}}, \code{\link{ModesDetect}},
  \code{\link{FirstRun}}, \code{\link{PostProc}}. }

\examples{
library(DelayActivation)
data(lung2); Time <- 0:10

## smoothing
mybasis <- create.bspline.basis(range(Time), length(Time)+4-2, 4, Time)
smoothing.par <- fdPar(mybasis, Lfdobj=2, lambda=5.0)
fitted.curves <- smooth.basis(Time, t(lung2), smoothing.par)[["fd"]]
lung2.smooth <- t(eval.fd(Time, fitted.curves))

act.info <- ActivationDetect(lung2.smooth, Time)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{clustering}
\keyword{misc}
