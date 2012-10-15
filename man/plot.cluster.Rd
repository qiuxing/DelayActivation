\name{plot.cluster}
\alias{plot.cluster}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Produce a spaghetti plot for a given cluster defined in \option{info.mat}. }

\description{ \code{plot.cluster} ...  }

\usage{ plot.cluster(x, Time, info.mat, clustname, ...)  }
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{An m by T dimensional data matrix, each row represents a gene
    and each column a time point.}
  \item{Time}{Time points. It must be a sequence of increasing numbers
    with length equal to the number of columns of x.}
  \item{info.mat}{The information matrix produced by either
  \code{FirstRun} or preferably, \code{PostProc}.}
  \item{clustname}{The unique name of the cluster.}
  \item{...}{Other parameters passed to \code{plot}.}
}

\details{ ...
}

\references{some references.}

\author{
  Xing Qiu, Hulin Wu, Shuang Wu.
}

\seealso{\code{\link{plot.delay}}, \code{\link{plotpdf.clusters}},
  \code{\link{FirstRun}}, \code{\link{PostProc}}.}

\examples{
library(DelayActivation)
data(lung2); Time <- 0:10

mybasis <- create.bspline.basis(range(Time), length(Time)+4-2, 4, Time)
smoothing.par <- fdPar(mybasis, Lfdobj=2, lambda=5.0)

\dontrun{
first.results <- FirstRun(lung2, Time, smoothing.par)
rr2 <- PostProc(first.results, cluster.min=60)

clust1 <- rr2$clust.info$Clust[1]
plot.cluster(lung2, Time, rr2, clust1,
                xlim=c(0,10), ylim=c(-3,3), 
                xlab="Time (day)", ylab="Expression")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{clustering}
\keyword{misc}
