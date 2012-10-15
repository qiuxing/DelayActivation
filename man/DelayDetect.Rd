\name{DelayDetect}
\alias{DelayDetect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Detect genes with delayed temporal activations.
}

\description{ \code{DelayDetect} ...  }

\usage{
DelayDetect(x, Time, delay.min=3, delay.max=5, beta1.cut=0.15, se.adj.cut=0.12, t.cut=3.0)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{An m by T dimensional data matrix, each row represents a gene
  and each column a time point.}
  \item{Time}{Time points. It must be a sequence of increasing numbers
  with length equal to the number of columns of x.}
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
}

\details{ For each time point in the range specified by
    \option{delay.min} and \option{delay.max}, each gene is fit to a
    linear regression from the starting time to that point.
    \option{beta} is the slope of the regression line; \option{se} is
    the standard error; \option{t} is the t-statistic.
    \option{se.adj} is defined as \code{se * sqrt{kk-1}}, where
    \code{kk} is the number of time points used in the regression
    analysis. This is to compensate for different sample size in these
    linear regression analysis.

    A gene is labeled as delay iff for at least one time point in the
    range,
    \enumerate{
      \item \code{se.adj < se.adj.cut} and
      \item either \code{abs(beta1)<beta1.cut} or \code{abs(tstat)<t.cut}.
    }
 }

\value{
  \item{delay}{A logic vector of length m, \code{delay[i]} indicates
  whether the ith gene has a delayed temporal pattern or not.}
  \item{delay.stats}{estimated beta, se, se.adj, and t for each gene on
    each day in the time window specified by \option{delay.min} and
    \option{delay.max}.}
  \item{cuts}{The cutoffs of beta, se.adj, and t used in this analysis.}
}

\references{some references.}

\author{
  Xing Qiu, Hulin Wu, Shuang Wu.
}

\seealso{\code{\link{plot.delay}}, \code{\link{ModesDetect}},
  \code{\link{ActivationDetect}}, \code{\link{FirstRun}},
  \code{\link{PostProc}}. }

\examples{
library(DelayActivation)
data(lung2); Time <- 0:10
mydelays <- DelayDetect(lung2, Time)
plot.delay(mydelays, 4)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{clustering}
\keyword{misc}
