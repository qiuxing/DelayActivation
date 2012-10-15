\name{plotpdf.clusters}
\alias{plotpdf.clusters}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ A wrapper of \code{plot.cluster} to produce one PDF file of
  spaghetti plots for all clusters defined in \option{info.mat}. }

\description{ \code{plotpdf.clusters} ...  }

\usage{ plotpdf.clusters(x, Time, info.mat, filename="clusters.pdf",
fig.size=3.0, ncol=4, ...)  }
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{An m by T dimensional data matrix, each row represents a gene
    and each column a time point.}
  \item{Time}{Time points. It must be a sequence of increasing numbers
    with length equal to the number of columns of x.}
  \item{info.mat}{The information matrix produced by either
  \code{FirstRun} or preferably, \code{PostProc}.}
  \item{filename}{Name of the output PDF file.}
  \item{fig.size}{size (in inches) of each individual spaghetti
  plot. Default to 3.0.}
  \item{ncol}{Number of subfigures per row. Default to 4. Number of rows
  is computed automatically.}
  \item{...}{Other parameters passed to function \code{plot.cluster}.}
}

\details{ ...
}

\references{some references.}

\author{
  Xing Qiu, Hulin Wu, Shuang Wu.
}

\seealso{\code{\link{plot.cluster}}}

\examples{
library(DelayActivation)
data(lung2); Time <- 0:10

mybasis <- create.bspline.basis(range(Time), length(Time)+4-2, 4, Time)
smoothing.par <- fdPar(mybasis, Lfdobj=2, lambda=5.0)

\dontrun{
first.results <- FirstRun(lung2, Time, smoothing.par)
rr2 <- PostProc(first.results, cluster.min=60)

plotpdf.clusters(lung2, Time, rr2, "firstrun.pdf",
                xlim=c(0,10), ylim=c(-3,3), 
                xlab="Time (day)", ylab="Expression")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{clustering}
\keyword{misc}
