#' Theoretical excursion extent (cross-surface) distribution for 2D Gaussian process
#'
#' @param l a value of area
#' @param u a user-specified thereshold value
#' @param gam bandwidth of Gaussian kernel
#' @returns right-tail probability
#' @export
#'
F2x = function(l, u, gam=4) {
  n <- 10000
  # denominator
  F1 = function(x){
    Lam <- matrix(rnorm(n*2),ncol=2)
    ind  <- (Lam[,1] < Lam[,2] & Lam[,2] < x/sqrt(2)) + 0
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2)) * ind) * dnorm(x)
  }
  # numerator
  F2 = function(x,l){
    Lam <- matrix(rnorm(n*2),ncol=2)
    ind1 <- (Lam[,1] < Lam[,2] & Lam[,2] < x/sqrt(2)) + 0
    ind2 <- (4*pi^2*(x-u)^2 >= l^2*(Lam[,2]-x/sqrt(2))*(Lam[,1]-x/sqrt(2))/(2*gam^4)) + 0
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2))*ind1*ind2)*dnorm(x)
  }
  vx <- seq(u,15,0.05)
  temp2 <- pracma::trapz(vx,sapply(vx,F2,l=l))
  temp1 <- pracma::trapz(vx,sapply(vx,F1))
  return(temp2/temp1)
}

#' Theoretical excursion volume (volume above cross-surface) distribution for 2D Gaussian process
#'
#' @param l a value of volume
#' @param u a user-specified thereshold value
#' @param gam bandwidth of Gaussian kernel
#' @returns right-tail probability
#' @export
#'
G2v = function(l, u, gam=4) {
  n <- 10000
  # denominator
  F1 = function(x){
    Lam <- matrix(rnorm(n*2),ncol=2)
    ind  <- (Lam[,1] < Lam[,2] & Lam[,2] < x/sqrt(2)) + 0
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2)) * ind) * dnorm(x)
  }
  # numerator
  F2 = function(x, l) {
    Lam <- matrix(rnorm(n*2),ncol=2)
    ind1 <- (Lam[,1] < Lam[,2] & Lam[,2] < x/sqrt(2)) + 0
    ind2 <- (pi^2*(x-u)^4 >= l^2*(Lam[,2]-x/sqrt(2))*(Lam[,1]-x/sqrt(2))/(2*gam^4)) + 0
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2))*ind1*ind2)*dnorm(x)
  }
  vx <- seq(u, 15, 0.05)
  temp2 <- pracma::trapz(vx, sapply(vx, F2, l=l))
  temp1 <- pracma::trapz(vx, sapply(vx, F1))
  return(temp2/temp1)
}

#' Calulate cross sectional area of \eqn{f(x,y) = u}
#'
#' @param data a dataframe/matrix of smoothed 3D data
#' @param u a user-specified thereshold value
#' @returns a list of all cross sectional area
#' @export
#' @seealso \code{F2x}
#' @examples
#'
#' gam <- 4
#' u <- 3.5
#' nr <- nc <- 1000   # 10000 * 10000
#' noise <- matrix(rnorm(nr * nc), nr, nc)
#' sdata <- smoothie::kernel2dsmooth(noise, "gauss", sigma=gam, nx=nr, ny=nc)
#' sdata <- sdata/sqrt(1/(4 * pi * gam^2))
#' xx <- cross_area(sdata, u=u)
#' xlim <- seq(0,30,0.5)
#' theoval <- sapply(xlim, F2x, u, gam)
#' plot(smooth.spline(xlim, theoval,df=5),lwd=2, type="l", ylab="Tail Probability", col="red",
#' ylim=c(0,1), xlab="Excursion Cross-sectional Area")
#' area <- ecdf.tail(xx$area$area)
#' lines(smooth.spline(area[,1], area[,2], df=5), lwd=2, col="black")
#' legend("topright", legend=c("Theoretical dist.", "Numerical dist."), col=c("red","black"), lwd=rep(2,2), bty="n")
#' ## detect peaks by cross-sectional surface area
#' pvalue <- sapply(xx$area$area, F2x, u=u, gam=gam)
#' cluster.id <- xx$area$group[which(pvalue <= 0.05)]
#' lapply(cluster.id, function(i) xx$cluster[xx$cluster$group == i,])
#'
cross_area = function(data, u) {
  x <- data >= u
  nc <- ncol(data)
  nr <- nrow(data)
  find.line = function(x) {
    # x is the data which contains only TRUE or FALSE
    columns <- which(colSums(x) > 0) # non-False columns
    f = function(column) {
      runs <- rle(x[,column]) # consecutive repeat
      runs_true <- which(runs$value==TRUE)
      end <- cumsum(runs$length)[runs_true]
      start <- cumsum(runs$length)[runs_true-1] + 1
      # bottom element is TRUE
      if (0 %in% (runs_true-1)) start <- c(1,start)
      cbind(column, start, end)
    }
    as.data.frame(do.call("rbind", lapply(columns, f)))
  }
  ##
  find.cluster = function(x) {
    # x is the output of function find.line
    match.line <- function(x1,x2) {
      for(i in 1:nrow(x2)) {
        ind <- ifelse(!(x1$start > x2[i,]$end | x1$end < x2[i,]$start), i, 0)
        if(ind!=0) break
      }
      return(x2[ind,])
    }
    find.line.neighbor = function(y){
      # given a line, find all its neighboring lines
      cluster <- y
      column <- y$column + 1
      while(column %in% x$column) {
        px <- x[x$column==column,]
        y <- match.line(y,px)
        if(nrow(y)==0) break # no match
        cluster <- rbind(cluster,y)
        column <- column + 1
      }
      return(cluster)
    }
    out <- NULL
    group <- 0
    columns <- unique(x$column)
    for(j in columns) {
      pat <- x[x$column==j,]
      if(nrow(pat) == 0) next
      for(k in 1:nrow(pat)) {
        group <- group + 1
        out <- rbind(out,cbind(find.line.neighbor(pat[k,]), group))
      }
      x <- dplyr::anti_join(x,out[,1:3],by=colnames(x)) # delete matched lines
    }
    return(out)
  }
  ##
  line2point = function(x) {
    # x is line-data
    temp <- rbind(cbind(x$column,x$end),cbind(rev(x$column),rev(x$start)))
    colnames(temp) <- c("column","row")
    temp <- as.data.frame(temp)
    return(cbind(temp[!duplicated(temp),], group=x$group[1]))
  }
  ##
  cross.area = function(x) {
    # x is line-data
    if ("start" %in% colnames(x)) x <- line2point(x)
    area <- ifelse(nrow(x)<3, 0, abs(pracma::polyarea(x[,1], x[,2])))
    return(data.frame(group=x$group[1], area=area))
  }
  ##
  outer = function(x) {
    # x is line-data
    t1 <- x[1,]
    t1$column <- t1$column - 1
    t2 <- x[nrow(x),]
    t2$column <- t2$column + 1
    x$start <- x$start -1
    x$end <- x$end + 1
    temp <- rbind(t1, x, t2)
    ## boundary issues
    temp <- temp[!(temp$column %in% c(0,nc+1)),] #column boundary
    temp$start[temp$start==0] <- 1
    temp$start[temp$start==(nr+1)] <- nr
    temp$end[temp$end==(nr+1)] <- nr # row boundary
    return(temp)
  }
  ##
  fit.point = function(x) {
    inner_point <- line2point(x)
    outer_point <- line2point(outer(x))
    euc.dist = function(x1, x2) sqrt(sum((x1 - x2)^2))
    interpolation = function(x){
      f1 = function(x) inner_point[apply(inner_point,1,euc.dist,x2=x) <= 1,]
      f2 = function(x1,y1,x2,y2) ((x2-x1)*u + x1*y2-x2*y1)/(y2-y1) # linear interpolation
      inner_match = f1(x)
      inpl_value = NULL
      x = as.numeric(x)
      for(i in 1:nrow(inner_match)) {
        t <- inner_match[i,]
        t <- as.numeric(t)
        iden_col <- (x[1] == t[1])
        if(euc.dist(t,x)==0) inpl_value <- rbind(inpl_value,x[1:2])
        else{
          temp <- ifelse(iden_col, f2(x[2], data[x[2], x[1]], t[2], data[t[2], t[1]]),
                         f2(x[1], data[x[2], x[1]], t[1], data[t[2], t[1]]))
          if(iden_col) inpl_value <- rbind(inpl_value, cbind(x[1], temp))
          else inpl_value <- rbind(inpl_value, cbind(temp, x[2]))
        }
      }
      colnames(inpl_value) <- c("column","row")
      return(as.data.frame(inpl_value, stringsAsFactors = FALSE))
    }
    if(nrow(x)==1) inner_point # single line cluster
    else as.data.frame(cbind(do.call(rbind, apply(outer_point, 1, interpolation)),
                             group=x$group[1]))
  }
  ##
  clusters <- find.cluster(find.line(x))
  outer_point <- do.call("rbind",lapply(split(clusters,clusters$group),outer))
  fit_point <- do.call("rbind",lapply(split(clusters,clusters$group),fit.point))
  area_lower <- do.call("rbind",lapply(split(clusters,clusters$group),cross.area))
  area_upper <- do.call("rbind",lapply(split(clusters,clusters$group),function(x) cross.area(outer(x))))
  area_fit <- do.call("rbind",lapply(split(clusters,clusters$group),function(x) cross.area(fit.point(x))))
  # return(list(cluster = clusters, outer_point = outer_point, fit_point = fit_point,
  #             area_lower = area_lower, area_upper = area_upper, area_fit = area_fit))
  return(list(cluster = clusters, fit_point = fit_point, area = area_fit))
}

#' Calculate volume above cross sectional of \eqn{f(x,y) >= u}
#'
#' @param output.cross_area output of the function \code{cross_area}
#' @param data Gaussian kernel smoothed data
#' @returns a data frame of all cross sectional volumes
#' @export
#' @seealso \code{cross_area}
#' @examples
#'
#' gam <- 4
#' u <- 3
#' nr <- nc <- 1000   # 10000 * 10000
#' noise <- matrix(rnorm(nr * nc), nr, nc)
#' sdata <- smoothie::kernel2dsmooth(noise, "gauss", sigma=gam, nx=nr, ny=nc)
#' sdata <- sdata/sqrt(1/(4 * pi * gam^2))
#' xx <- cross_area(sdata, u=u)
#' xlim <- seq(0, 30, 0.2)
#' theoval <- sapply(xlim, G2v, u=u, gam)
#' plot(smooth.spline(xlim, theoval, df=20),lwd=2, type="l", ylab="Tail Probability", col="red",
#' ylim=c(0,1), xlab="Excursion Volume", main=paste("u = ", u, sep=""))
#' yy <- cross_volume(xx, sdata)
#' volume <- ecdf.tail(yy$volume)
#' lines(smooth.spline(volume[,1], volume[,2], df=3), lwd=2, col="black")
#' legend("topright", legend=c("Theoretical dist.", "Numerical dist."),
#' col=c("red","black"), lwd=rep(2,2), bty="n")
#' plotly::plot_ly(z = ~ sdata, type="surface")
#' #' ## detect peaks by volume above surface
#' pvalue <- sapply(yy$volume, G2v, u=u, gam=gam)
#' cluster.id <- yy$group[which(pvalue <= 0.05)]
#' lapply(cluster.id, function(i) xx$cluster[xx$cluster$group == i,])
#'
cross_volume = function(output.cross_area, data) {
  ## fill all the points that are inside the boundary
  fill_point = function(x, data, clusters) {
    # x is boundary
    # isfit: logic value indicates if x is fitted points
    x0 <- x
    x <- clusters[clusters$group==x0$group[1],]
    t <- nrow(x)
    points = data.frame()
    for(i in 1:t) {
      temp <- x[i,]$start:x[i,]$end
      points <- rbind(points, cbind(rep(x[i,]$column, length(temp)), temp))
    }
    colnames(points) <- c("column", "row")
    out = cbind(points, z = data[cbind(points$row,points$column)]-u, group=x$group[1])
    out <- rbind(out,cbind(column=x0$column, row=x0$row, z=0, group=x0$group[1]))
    out[order(out$column),]
    return(out)
  }
  fit.fill <- do.call("rbind",lapply(split(output.cross_area$fit_point, output.cross_area$fit_point$group),
                                     fill_point, data = data, clusters=output.cross_area$cluster))
  ##
  getVolume = function(df) {
    # df is a dataframe from funtion fill_point
    if(nrow(df) < 3 | length(unique(df$column))==1 | length(unique(df$row))==1) out = NA
    else {
      # find triangular tesselation of (x,y) grid
      res <- geometry::delaunayn(as.matrix(df[,c("column","row")]), full=TRUE, options="Qz")
      # calulates sum of truncated prism volumes
      out <- sum(mapply(function(triPoints, A) A/3*sum(df[triPoints,"z"]),
                        split.data.frame(res$tri,seq_along(res$areas)),
                        res$areas))
    }
    return(as.data.frame(cbind(group = df$group[1], volume = out, peak.height = max(df$z))))
  }
  temp <- do.call("rbind", lapply(split(fit.fill, fit.fill$group), getVolume))
  temp[which(is.na(temp$volume)), "volume"] <- 0
  return(temp)
}
