\name{FirstRun}
\alias{FirstRun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A wrapper to run \code{DelayDetect}, \code{ModesDetect}, and
  \code{ActivationDetect} in one command.  This function is meant to be
  used as the first step in classification analysis, \code{PostProc} is
  the next step.
}

\description{ \code{FirstRun} ...  }

\usage{ FirstRun(x, Time, smoothing.par, delay.min=3, delay.max=5, beta1.cut=0.15, se.adj.cut=0.12, t.cut=3.0, delta=0.1, Q=0.5)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{An m by T dimensional data matrix, each row represents a gene
    and each column a time point.}
  \item{Time}{Time points. It must be a sequence of increasing numbers
    with length equal to the number of columns of x.}
  \item{smoothing.par}{The smoothing parameter used by function
    \code{smooth.basis()} from package \pkg{fda}.}
  \item{delay.min}{
    Minimum number of days used in detecting delay. Default to 3 units.
  }
  \item{delay.max}{
    Maximum number of days used in detecting delay. Default to 5 units.
  }
  \item{beta1.cut}{The cutoff for the estimated slope in the delay
    analysis. Default to 0.15.}
  \item{se.adj.cut}{The cutoff for the adjusted standard error in the delay
    analysis. Default to 0.12.}
  \item{t.cut}{The cutoff for the t-statistic in the delay
    analysis.  Default to 3.0.}
  \item{delta}{Error tolerance.  Default to 0.1.}
  \item{Q}{Quantile cutoff used in detecting the activation time.}
}

\details{ Please see the documentation of \code{DelayDetect},
  \code{ModesDetect}, and \code{ActivationDetect} for more details.}

\value{ A list of these components
  \item{fittded.curves}{Smoothed curves of each gene. It is a functional
    data (fd) object of smoothed gene expression (as functions of time)
    produced by function \code{smooth.basis()} from package \pkg{fda}.}
  \item{meancurs}{Smoothed mean curves of each cluster. It is a
    functional data (fd) object of smoothed gene expression (as
    functions of time) produced by function \code{smooth.basis()} from
    package \pkg{fda}.}
  \item{clust.info}{A data frame containing the results obtained from
    \code{\link{DelayDetect}}, \code{\link{ModesDetect}},
    \code{\link{ActivationDetect}}}.  An additional column,
  \option{Clust}, is provided for convenience. \option{Clust} is defined
  as \code{paste(Delay, UpDown, ActDay, Modes, sep=".")}, the unique
  identifier for clusters defined by these four criteria.
}

\references{some references.}

\author{
  Xing Qiu, Hulin Wu, Shuang Wu.
}

\seealso{\code{\link{DelayDetect}}, \code{\link{ModesDetect}},
  \code{\link{ActivationDetect}}, \code{\link{PostProc}},
  \code{\link{smooth.basis}}, \code{\link{plot.cluster}},
  \code{\link{plotpdf.clusters}}.}

\examples{
library(DelayActivation)
data(lung2); Time <- 0:10

mybasis <- create.bspline.basis(range(Time), length(Time)+4-2, 4, Time)
smoothing.par <- fdPar(mybasis, Lfdobj=2, lambda=5.0)

first.results <- FirstRun(lung2, Time, smoothing.par)
\dontrun{
plotpdf.clusters(lung2, Time, first.results, "firstrun.pdf",
                xlim=c(0,10), ylim=c(-3,3),
                xlab="Time (day)", ylab="Expression")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{clustering}
