#' Find local maxima and local minima of data sequence
#'
#' @param x numerical vector contains local maxima (minima)
#' @param partial logical value indicating if the two endpoints will be considered
#' @param decreasing logical value indicating whether to find local minima
#'
#' @return a vector of locations of local maxima or minima
#' @export
#'
#' @examples
#' a <- 100:1
#' which.peaks(a*sin(a/3))
#'
which.peaks = function(x, partial=FALSE, decreasing=FALSE){
  if (decreasing){
    if (partial){
      which(diff(c(FALSE,diff(x)>0,TRUE))>0)
    }else {
      which(diff(diff(x)>0)>0)+1
    }
  }else {
    if (partial){
      which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
    }else {
      which(diff(diff(x)>=0)<0)+1
    }
  }
}

#' Gaussiann process peak height
#'
#' Simulate the peak height and its density of the smoothed Gaussian process
#' @param n number of data points from a Gaussian distribution.
#' @param gam bandwidth of gausian kernel.
#' @param add.height logical value indicates if peak height should be plotted
#' @param add.density logical value indicates if density function should be plotted
#' @param ... arguments in \code{rnorm}, such as mean, sd
#' @returns location of peaks, plots of peak height and its density
#' @export
#' @examples
#' simu_peak_height(n = 100, add.height = TRUE, add.density = TRUE)
#'
simu_peak_height = function(n = 100, gam = 4, add.height = TRUE, add.density = TRUE, ...) {
  # generate smoothed gaussian process
  noise <- rnorm(n, ...)
  sdata <- smth.gau(noise, gam=gam)
  std <- sqrt(1/(2*gam*sqrt(pi))) # standard deviation
  sdata <- sdata/std
  # detect positive peaks
  peaks0 <- which.peaks(sdata)
  # plots
  if(add.height) {
    plot(sdata, type="l", lwd=2, xlab="", ylab="", main = "Peak Height")
    abline(h=0, lty=2)
    points(peaks0, sdata[peaks0], pch=21, bg="red", col="red")
    shape::Arrows(peaks0, rep(0,length(peaks0)), peaks0, sdata[peaks0]-0.01, code=3, col="blue",
                  arr.type="triangle",arr.adj=1,lwd=2, arr.width=0.2,arr.length=0.15)
  }
  if(add.density) {
    sdata <- smth.gau(rnorm(1000000,...), gam)
    sdata <- sdata/std # standardize
    peaks <- which.peaks(sdata)
    height <- sdata[peaks]
    # theoretical peak height density
    f.peak = function(x)
      sqrt(2/3)*dnorm(sqrt(3/2)*x) + sqrt(2*pi/3)*x*dnorm(x)*pnorm(x/sqrt(2))
    x <- seq(-3, 4, length.out=1000)
    value.th <- f.peak(x)
    ylim <- c(0, max(value.th, hist(height, freq=FALSE)$density)) * 1.2
    hist(height, breaks=60, freq=FALSE, xlab="", ylab="", main="", xlim=c(-3,4), ylim=c(0,0.5))
    lines(x, value.th, lwd=2, col="red", lty=1)
    legend("topleft",legend = c("simulated dist.", "theoretical dist."),
           lwd = c(2,2), col = c("black","red"), bty = "n")
  }
  return(peaks0)
}

#' Find the cross points of data sequence and a horizontal line
#'
#' @param data a vector
#' @param u a user-specified thereshold value
#' @param n number to data points to fit the cureve of \code{data}
#' @export
#' @returns locations of cross points
#'
cross_point = function(data, u, n=10000) {
  fit <- approx(1:length(data), data, n=n)
  above <- fit$y > u
  # TRUE if increasing, FALSE if decreasing
  intersect.points <- which(diff(above) != 0)
  # start with an increasing point and end with a decreasing point
  if (length(intersect.points) >= 2) {
    start <- which(!above[intersect.points])[1]
    end <- tail(which(above[intersect.points]),1)
    locations <- fit$x[intersect.points[start:end]]
  }
  else locations <- NULL
  return(locations)
}

#' compute tail probability based on empirical cdf
#'
#' @param x a numerical vector
#' @returns a matrix with two columns of x and tail probability \eqn{F(x)}
#'
ecdf.tail = function(x) {
  fecdf <- ecdf(x)
  x <- sort(x)
  cbind(x, 1 - fecdf(x))
}

#' Search endpoints of a peak extent
#'
#' @param x location of a peak
#' @inheritParams cross_point
#' @param data a numerical vector
#' @returns distance of a peak to its left and right endpoints
#' @export
#'
search_endpoint = function(x, u, data){
  # left-direction search
  lcount <- 1
  while(data[x-lcount]>u & lcount < x-1) {lcount <- lcount + 1}
  lcount <- ifelse(lcount==x-1, NA, lcount)
  lcount <- lcount - 1 + (data[x-lcount+1]-u)/(data[x-lcount+1] - data[x-lcount])
  # right-direction search
  rcount <- 1
  while(data[x+rcount]>u & rcount < length(data)-x){rcount <- rcount + 1}
  rcount <- ifelse(rcount==length(data)-x,NA,rcount)
  rcount <- rcount -1 + (data[x+rcount-1]-u)/(data[x+rcount-1] - data[x+rcount])
  return(c(lcount, rcount))
}

#' Theoretical pdf of peak extent
#'
#' @param x vector of quantities
#' @inheritParams cross_point
#' @param gam bandwidth of Gaussian kernel
#' @returns density of peak extent
#' @export
#'
f_x = function(x, u, gam) {
  lam1 <- 1/(2 * gam^2)
  lam2 <- 3/(4 * gam^4)
  if(lam1<=0 | lam2<=0) stop("lambda1 and lambda2 must be postive numbers")
  det <- lam2 - lam1^2
  den <- lam2*x^4 - 16*lam1*x^2+64
  if(det<=0) stop("lambda2 is smaller than lambda1 square")
  num <- lam2*x^2 - 8*lam1
  p1 <- sqrt(lam2/(2*pi))*(1-pnorm(sqrt(lam2/det)*u)) +
    lam1*dnorm(u)*(1-pnorm(-lam1*u/sqrt(det)))
  p2 <- 128*x*det/den^(3/2)*dnorm(8*u/sqrt(den))
  p3 <- (1+u^2*num^2/(den*det))*(1-pnorm(num*u/(sqrt(det*den))))
  p4 <- u*num/(sqrt(det*den))*dnorm(u*num/sqrt(det*den))
  return(1/p1*p2*(p3-p4))
}

#' Theoretical tail cdf of peak extent
#'
#' @inheritParams f_x
#' @inheritParams f_x
#' @inheritParams f_x
#' @seealso \code{f_x}
#' @returns right-tail probability of peak extent
#' @export
#'
Fx = function(x, u, gam) {
  lam1 <- 1/(2 * gam^2)
  lam2 <- 3/(4 * gam^4)
  if(lam1<=0 | lam2<=0) stop("lambda1 and lambda2 must be postive numbers")
  det <- lam2 - lam1^2
  den <- lam2*x^4 - 16*lam1*x^2 + 64
  if(det<=0) stop("lambda2 is smaller than lambda1 square")
  num <- lam2*x^2-8*lam1
  p1 <- sqrt(lam2/(2*pi))*(1-pnorm(sqrt(lam2/det)*u))
  p2 <- num/sqrt(den)*dnorm(8*u/sqrt(den))*(1-pnorm(num*u/sqrt(det*den)))
  p3 <- lam1*dnorm(u)*(1-pnorm(-lam1*u/sqrt(det)))
  return((p1-p2)/(p1+p3))
}

#' Theoretical pdf of peak mass
#'
#' @param v vector of quantities
#' @inheritParams f_x
#' @inheritParams f_x
#' @returns density of peak mass
#' @export
#'
g_v = function(v, u, gam) {
  lam1 <- 1/(2 * gam^2)
  lam2 <- 3/(4 * gam^4)
  if(lam1<=0 | lam2<=0) stop("lambda1 and lambda2 must be postive numbers")
  det <- lam2 - lam1^2
  p1 <- sqrt(lam2/(2*pi))*(1-pnorm(sqrt(lam2/det)*u))
  p2 <- lam1*dnorm(u)*(1-pnorm(-lam1*u/sqrt(det)))
  p3 <- 2048/(81*sqrt(2*pi*det)*v^5)
  f <- function(x) x^6*exp(-(lam1*(x+u)-32*x^3/(9*v^2))^2/(2*det))*dnorm(x+u)
  p4 <- integrate(f=f, lower=0, upper=Inf)$value
  # integrate vs. trapz
  return(1/(p1+p2)*p3*p4)
}

#' Theoretical tail cdf of peak mass
#'
#' @param V vector of quantiles
#' @inheritParams f_x
#' @inheritParams f_x
#' @seealso \code{g_v}
#' @returns right-tail probability of peak height
#' @export
#'
Gv = function(V, u, gam) {
  lam1 <- 1/(2 * gam^2)
  lam2 <- 3/(4 * gam^4)
  Z <- rnorm(5000)
  f1 = function(x) {
    mZ <- abs(sqrt(lam2-lam1^2)*Z - lam1*x)
    ind1 <- which(Z >= x*lam1/sqrt(lam2-lam1^2))
    mZ[ind1] <- 0
    return(mean(mZ*dnorm(x)))
  }
  f2 = function(x,V) {
    mZ <- abs(sqrt(lam2-lam1^2)*Z-lam1*x)
    ind1 <- which(Z >= x*lam1/sqrt(lam2-lam1^2))
    ind2 <- which(Z <= (x*lam1 - 32*(x-u)^3/(9*V^2))/sqrt(lam2-lam1^2))
    mZ[c(ind1,ind2)] <- 0
    return(mean(mZ)*dnorm(x))
  }
  vx = seq(u,15,0.01)
  temp1 = pracma::trapz(x=vx, y=sapply(vx,f1))
  temp2 = pracma::trapz(x=vx, y=sapply(vx,f2,V=V))
  return(temp2/temp1)
}

#' Gaussiann process peak extent
#'
#' Simulate the peak extent and its density of the smoothed Gaussian process
#' @inheritParams  simu_peak_height
#' @inheritParams  simu_peak_height
#' @param u threshold value
#' @param u.th threshold value for varifying the theoretical distribution of peak extent, supposed to be larger than 3
#' @param add.extent logical value indicates if plot the peak extent
#' @inheritParams simu_peak_height
#' @param add.cdf logical value indicates if plot the right tail probability distribution
#' @param ... arguments in \code{rnorm}, such as mean, sd
#' @returns locations of cross points, plots of peak extent and its right-tail probability
#' @export
#' @seealso [cross_point], [search_endpoint]
#' @examples
#' simu_peak_extent(n = 200, u = 1)
#'
simu_peak_extent = function(n = 100, gam = 4, u, u.th = 3.5, add.extent = TRUE,
                            add.density = TRUE, add.cdf = TRUE, ...) {
  noise <- rnorm(n,...)
  sdata <- smth.gau(noise, gam=gam)
  std <- sqrt(1/(2*gam*sqrt(pi))) # standard deviation
  sdata <- sdata/std
  crossx <- cross_point(sdata, u, n=10000)
  if (is.null(crossx)) stop("u is too large, see the plot of (smoothed) data")
  # plots
  peaks <- which.peaks(sdata)
  peaks <- peaks[sdata[peaks]>=u]
  if (add.extent) {
    plot(sdata, type="l", lwd=2, xlab="", ylab="")
    abline(h=0)
    abline(h = u, lty=2, col = "darkgray")
    points(peaks, sdata[peaks], pch=21, bg="red", col="red")
    cross1 = crossx[seq(1,length(crossx),2)]
    cross2 = crossx[seq(2,length(crossx),2)]
    shape::Arrows(cross1, u, cross2, u, code=3, col="blue",
           arr.type="triangle", arr.adj=1, lwd=2,
           arr.width=0.15, arr.length=0.1)
  }
  if (add.density | add.cdf) {
    u <- u.th
    noise <- rnorm(10000000,...)
    sdata <- smth.gau(noise, gam)
    sdata <- sdata/sd(sdata) # standardize
    peaks <- which.peaks(sdata)
    upeaks <- peaks[sdata[peaks]>u]
    width <- sapply(lapply(upeaks, search_endpoint, u=u, data=sdata), sum)
    tail <- ecdf.tail(width)
    if(add.cdf) {
      vx <- seq(0, 15, 0.1)
      plot(vx, sapply(vx,Fx, u=u, gam=gam), ylim=c(0,1), col="red", lwd=2, type="l",
           xlab="Excursion Extent", ylab="Tail Probability", main=paste("u = ", u, sep=""))
      lines(smooth.spline(tail[,1], tail[,2], df=10), lwd=2, col="black")
      legend("topright", legend=c("theoretical dist.", "simulated dist."),
             col=c("red","black"), lwd=rep(2,2), bty="n")
    }
    if(add.extent) {
      vx <- seq(0, 15, 0.1)
      value.th <- fx(vx, u=u,gam=gam)
      ylim <- c(0, max(value.th, hist(width, freq=FALSE)$density)) * 1.2
      hist(width, ylim=ylim, freq=FALSE, xlab="", ylab="", main="")
      lines(vx, value.th, lwd=2, col="red", lty=1)
      legend("topright",legend = c("simulated dist.", "theoretical dist."),
             lwd = c(2,2), col = c("black","red"), bty = "n")
    }
  }
  return(crossx)
}

#' Calculate the area under a curve
#'
#' @inheritParams search_endpoint
#' @inheritParams search_endpoint
#' @inheritParams search_endpoint
#' @returns a vector of area under the curve
#' @export
#' @seealso \code{search_endpoint}
#'
auc = function(x, u, data) {
  endpoints <- search_endpoint(x, u, data)
  temp <- (x-round(endpoints[1])) : (x+round(endpoints[2]))
  traparea <- pracma::trapz(x=temp, y=data[temp])- u *(tail(temp,1)-temp[1])
  traparea <- ifelse(is.na(traparea), 0, traparea)
  return(traparea)
}

#' Peak detection based on peak height, extent and mass
#'
#' @param data a data sequence, if not NULL, other parameters generating signal will not be used
#' @param l length of simulated data sequence
#' @param loc a vector of locations of peaks
#' @param sigma a vecctor of standard deviation for each peak, same length as \code{loc}
#' @param scale a vector of scale for each peak, same length as \code{loc}
#' @inheritParams simu_peak_height
#' @inheritParams simu_peak_extent
#' @param alpha significant level for peak detection
#' @param by detect peaks by which, a chracter of \code{height}, \code{extent} or \code{mass}
#' @param ... arguments in \code{rnorm}, such as mean, sd
#' @returns locations of significant peaks by \code{height}, \code{extent} and \code{mass}
#' @export
#' @seealso [cross_point], [search_endpoint]
#' @examples
#' loc=c(100,250,400,650,800,900)
#' sigma=c(4,10,15,20,10,2)
#' scale=c(0.4,0.5,1.2,1.8,1.2,1.2)
#' gam = 4
#' simu_peak_detect(l=1000, loc=loc, sigma=sigma, scale=scale, gam = gam, u=0.02, by="height")
#'
simu_peak_detect = function(data = NULL, l, loc, sigma, scale, gam=4, u, alpha=0.05,
                            by=c("height","extent","mass"), ...) {
  if (is.null(data)) {
    gen.signal = function(l, loc, sigma, scale){
      out <- rep(0, l)
      n <- length(loc)
      for (i in 1:n) {
        range <- (loc[i]-3*sigma[i]) : (loc[i] + 3*sigma[i])
        out[range] <- scale[i] * dnorm(range, loc[i], sigma[i])
      }
      return(out)
    }
    signal <- gen.signal(l=l, loc=loc, sigma=sigma, scale=scale)
    noise <- rnorm(l, sd=0.02)
    data <- signal + noise
    sdata <- smth.gau(data, gam)
    plot(data, type="l", col="darkseagreen3", xlab="", ylab="")
    abline(h=0)
    lines(signal,col="red",lwd=2)
    lines(sdata,col="black",lwd=2)
    legend("topleft", legend=c("original data", "signal", "smoothed data"), lwd = rep(2,3),
           col=c("darkseagreen3", "red", "black"), bty="n")
  }
  else {
    plot(data, type="l", col="darkseagreen3", xlab="", ylab="")
    abline(h=0)
    lines(sdata,col="black",lwd=2)
    legend("topleft", legend=c("original data","smoothed data"), lwd = rep(2,2),
           col=c("darkseagreen3", "black"), bty="n")
  }
  sdata <- smth.gau(data, gam)
  peaks <- which.peaks(sdata)
  upeaks <- peaks[sdata[peaks]>u]
  # plot data
  col1 <- "cyan"
  col2 <- "blue"
  # plot peak detection
  if (by == "height") {
    plot(sdata, type="l", col="black", lwd=2, xlab="", ylab="")
    abline(h=0, col="grey")
    abline(h=u, col="deeppink1", lty=2)
    peaks_dif <- setdiff(peaks, upeaks)
    points(peaks_dif, sdata[peaks_dif], pch=24, bg=col1, col="black")
    points(upeaks, sdata[upeaks], pch=22, bg=col2, col="black")
    legend("topleft", legend=c("candidate peaks","significant peaks"), pch=c(24,22), bg=c(col1,col2), bty="n")
    return(upeaks)
  }
  # plot peak extent
  if (by == "extent") {
    crossx <- cross_point(sdata, u)
    cross1 <- crossx[seq(1,length(crossx),2)]
    cross2 <- crossx[seq(2,length(crossx),2)]
    width <- sapply(lapply(upeaks, search_endpoint, u=u, data=sdata), sum)
    pvalue <- sapply(width, Fx, u=u, gam=gam)
    sig.peaks <- upeaks[pvalue <= alpha]
    plot(sdata, type="l", col="black", lwd=2, xlab="", ylab="")
    abline(h=0)
    abline(h=u, col="deeppink1", lty=2)
    peaks_dif <- setdiff(upeaks, sig.peaks)
    cdt.width1 <- cross1[which(upeaks %in% sig.peaks)]
    true.width1 <- cross1[-which(upeaks %in% sig.peaks)]
    cdt.width2 <- cross2[which(upeaks %in% sig.peaks)]
    true.width2 <- cross2[-which(upeaks %in% sig.peaks)]
    shape::Arrows(cdt.width1, u, cdt.width2, u, code=3, col=col2,
                  arr.type="triangle", arr.adj=1, lwd=3,
                  arr.width=0, arr.length=0.0)
    shape::Arrows(true.width1, u, true.width2, u, code=3, col=col1,
                  arr.type="triangle", arr.adj=1, lwd=3,
                  arr.width=0, arr.length=0.0)
    legend("topleft", legend=c("candidate peaks","significant peaks"), lwd=c(2,2),
           col=c(col1,col2), bty="n")
    return(sig.peaks)
  }
  # plot peak mass
  if (by == "mass") {
    plot(sdata, type="l", col="black", lwd=2, xlab="", ylab="")
    abline(h=0)
    abline(h=u, col="deeppink1", lty=2)
    area <- sapply(upeaks, auc, u=u, data=sdata)
    upeaks <- upeaks[!duplicated(area)]
    area <- area[!duplicated(area)]
    pvalue <- sapply(area, Gv, u=u, gam=gam)
    sig.peaks <- upeaks[pvalue <= alpha]
    peaks_dif <- setdiff(upeaks, sig.peaks)
    shade.area = function(data, xmn, xmx, n=10000, col="grey"){
      fit <- approx(1:length(data), data, n=n)
      xmn.ind <- which(fit$x==xmn)
      xmx.ind <- which(fit$x==xmx)
      if (fit$y[xmn.ind] < u) xmn.ind <- xmn.ind + 1
      if (fit$y[xmx.ind] < u) xmx.ind <- xmx.ind - 1
      shade.x <- fit$x[xmn.ind:xmx.ind]
      shade.y <- fit$y[xmn.ind:xmx.ind]
      polygon(shade.x, shade.y, col=col)
    }
    for (i in 1:length(peaks_dif)){
      lx <- 2 * which(upeaks %in% peaks_dif) -1
      rx <- 2 * which(upeaks %in% peaks_dif)
      shade.area(sdata, xmn=crossx[lx[i]], xmx=crossx[rx[i]], col=col2)
    }
    if (length(sig.peaks) !=0) {
      for (i in 1:length(sig.peaks)){
        lx <- 2 * which(upeaks %in% sig.peaks) -1
        rx <- 2 * which(upeaks %in% sig.peaks)
        shade.area(sdata, xmn=crossx[lx[i]], xmx=crossx[rx[i]], col=col1)
      }
      legend("topleft", legend=c("candidate peaks","significant peaks"), lwd=c(2,2),
             col=c(col1,col2), bty="n")
      return(sig.peaks)
    }
    else return(NULL)
  }
}







