## Hulin's ad hoc method.  See /Consulting/Hulin_Wu/Micro_Array-Time_Course/stepwise-clustering/delay_and_activation for more details.

DelayDetect <- function(x, Time, delay.min=3, delay.max=5, beta1.cut=0.15, se.adj.cut=0.12, t.cut=3.0){
### we check the "flatness" of time course microarray data (x)
### from day "delay.min" up to day "delay.max".  "Time" is assumed to
### be monotonic increasing.
  ngenes <- dim(x)[1]; genenames <- rownames(x)
  kk <- which(Time >= delay.min & Time <= delay.max)
  timenames <- paste("Day", Time[kk], sep="")
### First, let's re-center data so every expression on first day is
### always zero
  x0 <- sweep(x, MARGIN=1, x[,1], "-")

###  Delay detection is done by a series of linear regressions between
###  the expressions (centered on Day 0) and Time, from Day delay.min
###  to Day delay.max. Note these regressions do not have an intercept
###  term.
  beta1 <- matrix(0, ngenes, length(kk))
  colnames(beta1) <- timenames
  se <- matrix(0, ngenes, length(kk))
  colnames(se) <- timenames
  for (k in kk){
    t.k <- Time[1:k]; y.k <- t(x0[,1:k])
    rr.k <- coef(summary(lm(y.k~0+t.k)))
    beta1[, paste("Day",Time[k],sep="")] <- sapply(1:ngenes, function(x) rr.k[[x]][,"Estimate"])
    se[, paste("Day",Time[k],sep="")] <- sapply(1:ngenes, function(x) rr.k[[x]][,"Std. Error"])
  }

### t statistic is defined as beta1/se.
  tstat <- beta1/se
### we need to adjust the standard error by sqrt(k-1) so to make sure
### SE computed on different days are comparable.
  se.adj <- sweep(se, MARGIN=2, sqrt(kk-1), "*")

### A gene is declared a delayed gene if for at least one day in the
### range of (delay.min, delay.max) the following two conditions hold:
### 1) the adjusted standard error of beta1 (se.adj) is less than
### se.adj.cut; 2) either the absolute value of t-statistic is less than
### t.cut or the absolute value of beta1 is less than beta1.cut.
  DELAY.kk <- (se.adj<se.adj.cut) & ( (abs(beta1)<beta1.cut) | (abs(tstat)<t.cut) )
  ## use character strings instead of T/F for readability
  DELAY <- ifelse(rowSums(DELAY.kk)>0, "Delay", "Regular")
  delay.stats <- abind(beta1, se, se.adj, tstat, along=3)
  dimnames(delay.stats) <- list(genenames, timenames, c("beta1", "se", "se.adj", "tstat"))

  return(list("delay"=DELAY, "delay.stats"=delay.stats, "Time"=Time,
              "cuts"=c("beta1.cut"=beta1.cut, "se.adj.cut"=se.adj.cut, "t.cut"=t.cut)))
}

plot.delay <- function(delayinfo, day, manual.delay=NULL, manual.regular=NULL, ...){
  ## A convenient function to plot delay detections on a particular day.
  DELAY <- delayinfo$delay=="Delay"
  ## take the setdiff between manual.delay and delay
  manual.delay2 <- intersect(which(!DELAY), manual.delay)
  manual.regular2 <- intersect(which(DELAY), manual.regular)
  cuts <- delayinfo$cuts
  beta1.cut <- cuts["beta1.cut"]; se.adj.cut <- cuts["se.adj.cut"]; t.cut <- cuts["t.cut"]
  mycoefs <- delayinfo$delay.stats[, paste("Day", day, sep=""), c("se.adj", "beta1")]
  plot(mycoefs[!DELAY,], ...)
  if (!is.null(DELAY)) points(mycoefs[DELAY,], col=2)
  if (!is.null(manual.delay2)) points(mycoefs[manual.delay2,], col=3, pch=2)
  if (!is.null(manual.regular2)) points(mycoefs[manual.regular2,], col=4, pch=3)
  abline(0, t.cut/sqrt(day-1), col=2); abline(0, -t.cut/sqrt(day-1), col=2)
  abline(v=se.adj.cut, col=2, lty=2)
  abline(h=beta1.cut, col=2)
  abline(h=-beta1.cut, col=2)
}


## Detect number of modes
.nmodes <- function(x, delta=0.1) {
  ## determine "roughly" how many modes (local minima/maxima) a given gene has
  diffs <- diff(x)
  signs <- sign(sapply(diffs, function(xi) ifelse(abs(xi)>delta, xi, 0)))
  ## if sign == 0 we move to the next t.
  signs <- signs[signs !=0]
  return(sum(abs(diff(signs)) == 2))
}
ModesDetect <- function(x.smoothed, delta=0.1){
  ## this function should be applied to the smoothed data
  genenames <- rownames(x.smoothed)
  mymodes <- apply(x.smoothed, 1, .nmodes, delta=delta)
  names(mymodes) <- genenames
  return (mymodes)
}


## Detect the activation day and up/down based on the first day
## expression level reaches Q (default: 50\%) of its global optima.
.activation <- function(v, Time, Q=0.5){
  vmin <- min(v); vmax <- max(v)
  if (vmax>=abs(vmin)) {                #global maximum
    vsign <- 1
    for (i in 2:length(Time)){
      if (v[i] > Q*vmax) {
        act.ind <- i; break
      }}
  } else {                              #looking for minimum instead
    vsign <- -1
    for (i in 2:length(Time)){
      if (v[i] < Q*vmin) {
        act.ind <- i; break
      }}}
  act.days <- Time[act.ind]
  return(c("vsign"=vsign, "Activation.Day"=act.days, "Activation.Index"=act.ind))
}

ActivationDetect <- function(x.smoothed, Time, Q=0.5){
  ## this function should be applied to the smoothed data
  genenames <- rownames(x.smoothed)
  x.smoothed0 <- sweep(x.smoothed, MARGIN=1, x.smoothed[,1], "-")
  rr <- t(apply(x.smoothed0, 1, .activation, Time=Time, Q=Q))
  ## use character strings instead of -1/1 for readability
  UpDown <- ifelse(rr[, "vsign"]==1, "Up", "Down")

  activation.info <- data.frame("UpDown"=UpDown, "Activation.Day"=rr[, "Activation.Day"], "Activation.Index"=rr[, "Activation.Index"])
  rownames(activation.info) <- genenames
  return(activation.info)
}


.actperiod <- function(v, Time, Q=0.5){
  vmin <- min(v); vmax <- max(v)
  if (vmax>=abs(vmin)) {                #global maximum
    vsign <- 1
    act.period.ind <- v > Q*vmax
  } else {                              #looking for minimum instead
    vsign <- -1
    act.period.ind <- v < Q*vmin
  }
  return(list("vsign"=vsign, "ActPeriod.Index"=act.period.ind))
}

ActivationPeriodDetect <- function(x.smoothed, Time, Q=0.5){
  ## reports a period of activation instead of just the first day
  genenames <- rownames(x.smoothed)
  x.smoothed0 <- sweep(x.smoothed, MARGIN=1, x.smoothed[,1], "-")
  UpDown <- c()
  actperiod.days <- matrix(0, nrow=length(genenames), ncol=length(Time))
  rownames(actperiod.days) <- genenames
  colnames(actperiod.days) <- paste("Day",Time,sep="")
  for (gg in genenames){
    rr <- .actperiod(x.smoothed0[gg,], Time=Time, Q=Q)
    UpDown[gg] <- ifelse(rr[["vsign"]]==1, "Up", "Down")
    actperiod.days[gg,] <- rr[["ActPeriod.Index"]]
  }
  return(data.frame(UpDown, actperiod.days))
}

FirstRun <- function(x, Time, smoothing.par, delay.min=3, delay.max=5, beta1.cut=0.15, se.adj.cut=0.12, t.cut=3.0, delta=0.1, Q=0.5){
  ## A convenient wrapper to do the first run in one step
  fitted.curves <- smooth.basis(Time, t(x), smoothing.par)[["fd"]]
  x.smoothed <- t(eval.fd(Time, fitted.curves))
  Delay <- DelayDetect(x, Time, delay.min=delay.min, delay.max=delay.max, beta1.cut=beta1.cut, se.adj.cut=se.adj.cut, t.cut=t.cut)[["delay"]]
  num.modes <- ModesDetect(x.smoothed, delta=delta)
  Modes <- paste("Modes", num.modes, sep="")
  act.info <- ActivationDetect(x.smoothed, Time, Q=Q)
  act.day <- act.info[,"Activation.Day"]
  UpDown <- act.info[,"UpDown"]
  ActDay <- paste("ActDay",act.day,sep="")
  Clust <- paste(Delay, UpDown, ActDay, Modes, sep=".")

  clust.info <- data.frame(Delay=Delay,
                           Activation.Day=act.day,
                           UpDown=UpDown,
                           Modes=num.modes,
                           Clust=Clust)

  params <- list(delay.min=delay.min, delay.max=delay.max,
                 beta1.cut=beta1.cut, se.adj.cut=se.adj.cut,
                 t.cut=t.cut, delta=delta, Q=Q)

  ## Compute the mean curves
  meancurs <- NULL
  cnames <- sort(unique(Clust))
  for (clust in cnames){
    cl.k <- which(Clust==clust)
    meancurs[[clust]] <- mean.fd(fitted.curves[cl.k])
  }

  return(list(fitted.curves=fitted.curves,
              meancurs=meancurs,
              clust.info=clust.info,
              params=params))
}


.merge.clust <- function(fvec, min.size){
  ## given a factor fvec which indicates the cluster membership (such
  ## as modes and/or activation days) and min.size, re-label this
  ## vector so that the left-most and right-most clusters are merged
  ## to the middle, and the merged clusters all have size >= min.size
  fvec2 <- fvec
  clust.names <- sort(unique(fvec)); nc <- length(clust.names)
  clust.sizes <- table(fvec)[as.character(clust.names)]
  ## starting from the left (small number of modes/activation day).
  ## left.merge.ind is the smallest index of cluster with a
  ## cumulative size greater equal to min.size
  left.merge.ind <- min((1:nc)[cumsum(clust.sizes)>=min.size])
  if (left.merge.ind>1){
    left.merge.clust <- clust.names[left.merge.ind]
    fvec2 <- replace(fvec, fvec<left.merge.clust, left.merge.clust)
  }
  ## starting from the right (large number of modes).
  right.merge.ind <- max((1:nc)[rev(cumsum(rev(clust.sizes)))>=min.size])
  if (right.merge.ind<nc){
    right.merge.clust <- clust.names[right.merge.ind]
    fvec2 <- replace(fvec2, fvec2>right.merge.clust, right.merge.clust)
  }
  return(fvec2)
}

## A utility function needed by PostProc
.splitclustname <- function(clusts) {
  cl2 <- as.character(clusts)
  tokens <- strsplit(cl2, '.', fixed=TRUE)
  tokens2 <- do.call(rbind, tokens)
  Activation.Day <- as.integer(substring(tokens2[,3], 7))
  Modes <- as.integer(substring(tokens2[,4], 6))
  return (data.frame(Delay=tokens2[,1], Activation.Day=Activation.Day,
                     UpDown=tokens2[,2], Modes=Modes, Clust=cl2))
}

PostProc <- function(firstrun.info, manual.delay=NULL,
                     manual.regular=NULL, act.min=30, modes.min=50,
                     cluster.min=10){
  ## post-process results generated by DelayActivation.  You can
  ## specify manual overrides here.
  fitted.curves <- firstrun.info[["fitted.curves"]]
  clust.info <- firstrun.info[["clust.info"]]
  genenames <- rownames(clust.info)
  params <- firstrun.info[["params"]]
  Delay <- clust.info[, "Delay"]
  if (!is.null(manual.delay)) Delay[manual.delay] <- "Delay"
  if (!is.null(manual.regular)) Delay[manual.regular] <- "Regular"
  ## act.min and modes.min specifies the minimum size of activation
  ## cluster and modes cluster.  The actual algorithm is a bit
  ## complicated so please see .merge.clust() for details.
  num.modes2 <- .merge.clust(clust.info[, "Modes"], modes.min)
  ## Apply the same algorithm to activation days
  act.day2 <- .merge.clust(clust.info[, "Activation.Day"], act.min)
  ## After determining the new act/mode-clusters, all genes will be
  ## given a unique cluster based on merged c(delay, updown, act.day,
  ## modes).
  Modes <- paste("Modes", num.modes2, sep="")
  ActDay <- paste("ActDay", act.day2, sep="")
  Clust <- paste(Delay,
                 clust.info[,"UpDown"],
                 ActDay, Modes, sep=".")
  cnames <- sort(unique(Clust))
  csize <- sapply(cnames, function(cn) sum(Clust==cn))
  ## since we may have manually selected delay/regular genes and merged some modes groups, we have to re-compute the mean curves
  meancurs <- NULL
  for (clust in cnames){
    cl.k <- which(Clust==clust)
    meancurs[[clust]] <- mean.fd(fitted.curves[cl.k])
  }

  ## Now if such a cluster has less than cluster.min number of genes,
  ## it will be merged to the "nearest" cluster determined by L2
  ## distance
  small.clusters <- cnames[csize<cluster.min]
  large.clusters <- cnames[csize>=cluster.min]
  Clust2 <- Clust
  if (!is.null(small.clusters)){
    for (sc in small.clusters){
      L2vec <- sapply(large.clusters, function(cn) inprod(meancurs[[sc]]-meancurs[[cn]], meancurs[[sc]]-meancurs[[cn]]))
      closest.large.clust <- large.clusters[which.min(L2vec)]
      Clust2 <- replace(Clust2, Clust2==sc, closest.large.clust)
    }
  }

  ## re-compute all the information based on the merged clusters
  clust.info2 <- .splitclustname(Clust2)
  ## ## group all information together and return
  ## clust.info2 <- data.frame(Delay=Delay,
  ##                           Activation.Day=act.day2,
  ##                           UpDown=clust.info[,"UpDown"],
  ##                           Modes=num.modes2,
  ##                           Clust=Clust2)
  rownames(clust.info2) <- genenames

  ## Re-compute the new mean curves to reflect merges
  meancurs <- NULL
  cnames <- sort(unique(Clust2))
  for (clust in cnames){
    cl.k <- which(Clust2==clust)
    meancurs[[clust]] <- mean.fd(fitted.curves[cl.k])
  }

  return(list(fitted.curves=fitted.curves,
              meancurs=meancurs,
              clust.info=clust.info2,
              params=params))
}

plot.cluster <- function(x, Time, info.mat, clustname, ...){
  ## plot an individual cluster.
  fitted.curves <- info.mat$fitted.curves
  Clust <- info.mat$clust.info$Clust
  cl <- which(Clust==clustname); csize <- length(cl)
  T.grid <- seq(min(Time), max(Time), 0.2)

  meancur <- eval.fd(T.grid, mean.fd(fitted.curves[cl]))
  plot(0, type="n", main=paste(clustname, "(n=", csize, ")", sep=""), ...)
  for (i in cl) lines(Time, x[i,], col="grey", lwd=0.2)
  lines(T.grid, meancur, lwd=2.0, col="red")
}


plotpdf.clusters <- function(x, Time, info.mat, filename="clusters.pdf", fig.size=3.0, ncol=4, ...){
  ## where info.mat is one returned by either FirstRun() or PostProc()
  Clust <- info.mat$clust.info$Clust
  cnames <- sort(unique(Clust))

  fig.rows <- ceiling(length(table(Clust))/ncol)
  pdf(filename, width=ncol*fig.size, height=fig.rows*fig.size)
  par(mfrow=c(fig.rows, ncol))
  for (clustname in cnames){
    plot.cluster(x, Time, info.mat, clustname, ...)
  }
  dev.off()
}

## plot.clusters(lung, Time, rr2, "clusters-lung.pdf", xlim=c(0,10), ylim=c(-3,3), xlab="Time (day)", ylab="Expression")
