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
#' @param gam bandwidth of Gaussian kernel
#'
#' @return vector of smoothed values
#' @export
#' @examples
#' smth.gau(x=rnorm(1000), gam=20)
#'
smth.gau = function(x, gam){
  # Gaussian kernel
  .kern = function(x,v=1){
    temp = exp(-x^2/(2*v^2))
    return(temp/sum(temp))}
  k = ifelse(2*6*gam <= 0.9*length(x),6,floor(0.9*length(x)/(2*gam))) #6sigma
  Lwindow = (-k*gam):(k*gam)
  w = .kern(Lwindow,v=gam)
  sx = conv(x,w,"same")
  # adjusted weights
  adj.w = function(w){
    hw = floor(length(w)/2)
    a = cumsum(w)[(hw+1):length(w)]
    return(c(a,rep(1,length(x)-2*length(a)),rev(a)))
  }
  return(sx/adj.w(w))
}
