#' Compute convolution function using FFT
#' @description Compute convolution function using FFT, similar to \code{'conv'} in matlab
#' @references Matlab document on \code{'conv'}: \url{https://www.mathworks.com/help/matlab/ref/conv.html}
#' @param u numerical vector
#' @param v numerical vector, don't need to have the same length as \code{u}
#' @param shape if 'same', return central part of the convolution and has the same size as \code{u};
#'              otherwise return the whole sequence of size \eqn{lenth(u)+length(v)-1}.
#'
#' @return a vector of convolution, as specified by shape.
#' @export
#' @examples
#' u = c(-1,2,3,-2,0,1,2)
#' v = c(2,4,-1,1)
#' w = conv(u,v,'same')
conv = function(u, v, shape = c("same","full")){
  shape <- match.arg(shape)
  lx <- length(u)
  ly <- length(v)
  n <- lx + ly - 1
  w <- fft(fft(c(u,rep(0,n-lx))) * fft(c(v,rep(0,n-ly))),inverse=TRUE)/n
  w <- Re(w)
  if(shape=="same") w <- w[floor(ly/2+1):(floor(ly/2)+lx)]
  return(w)
}

#' Smoothing data using Gaussian kernel
#'
#' @param x numeric vector of values to smooth
#' @param gamma bandwidth of Gaussian kernel
#'
#' @return vector of smoothed values
#' @export
#' @examples
#' smth.gau(x=rnorm(1000), gamma=20)
#'
smth.gau = function(x, gamma){
  # Gaussian kernel
  .kern = function(x,v=1){
    temp = exp(-x^2/(2*v^2))
    return(temp/sum(temp))}
  k = ifelse(2*6*gamma <= 0.9*length(x),6,floor(0.9*length(x)/(2*gamma))) #6sigma
  Lwindow = (-k*gamma):(k*gamma)
  w = .kern(Lwindow,v=gamma)
  sx = conv(x,w,"same")
  # adjusted weights
  adj.w = function(w){
    hw = floor(length(w)/2)
    a = cumsum(w)[(hw+1):length(w)]
    return(c(a,rep(1,length(x)-2*length(a)),rev(a)))
  }
  return(sx/adj.w(w))
}

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
#' a = 100:1
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
