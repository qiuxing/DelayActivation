\name{plot.delay}
\alias{plot.delay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Produce some useful plots for detecting genes with delayed temporal patterns.
}

\description{ \code{plot.delay} ...  }

\usage{
plot.delay(delayinfo, day, manual.delay=NULL, manual.regular=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{delayinfo}{
    A matrix produced by function \code{DelayDetect}.
  }
  \item{day}{The day of delay plot. It must be a time point within range
  of minimum and maximum days defined in the delay detection analysis.}
  \item{manual.delay}{A list of the index of genes to be set to "Delay"
    manually. The actual code is: \code{Delay[manual.delay] <- "Delay"}.}
  \item{manual.regular}{A list of the index of genes to be set to "Regular"
    manually. The actual code is: \code{Delay[manual.regular] <- "Regular"}.}

  \item{...}{Other parameters passed to function \code{plot}.}
}

\details{ For each time point in the range specified by
  \option{delay.min} and \option{delay.max}, each gene is fit to a
  linear regression from the starting time to that point.
  \option{beta} is the slope of the regression line; \option{se} is
  the standard error; \option{t} is the t-statistic.  \option{se.adj}
  is defined as \code{se * sqrt{kk-1}}, where \code{kk} is the number
  of time points used in the regression analysis. This is to
  compensate for different sample size in these linear regression
  analysis.

  A gene is labeled as delay iff for at least one time point in the
  range,
  \enumerate{
    \item \code{se.adj < se.adj.cut} and
    \item either \code{abs(beta1)<beta1.cut} or \code{abs(tstat)<t.cut}.
  }

  \code{plot.delay} plots \option{se.adj} versus \option{beta} and
  illustrates the three criteria employed to detect delays by red lines.
}

\references{some references.}

\author{
  Xing Qiu, Hulin Wu, Shuang Wu.
}

\seealso{\code{\link{DelayDetect}}}

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
