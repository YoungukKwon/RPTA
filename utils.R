library(ggplot2)

### square norm (L2) of vectors
normVector <- function(x) {
  sqrt(sum(x^2))
}

### soft-thresholding function 
soft_thres <- function(x, lambda) {
  tmp = abs(x) - lambda/2
  res = sign(x) * ifelse(tmp > 0, tmp, 0)
  res
}


### generate gaussian process with given covariance 
gaussprocess <- function(from = 0, to = 1, K = function(s, t) {min(s, t)},
                         start = 0, nt = 101) {
  
  t <- seq(0,nt-1)/(nt-1)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  path <- MASS::mvrnorm(mu = rep(0, times = nt), Sigma = Sigma)
  path <- path - path[1] + start  # Must always start at "start"
  
  return(data.frame("t" = t, "xt" = path))
}

### Huber's weight function (From RobRSVD package) ###
huberWeightLS <- function(data, k) {
  w = k/abs(data)
  w[w > 1] = 1
  return(w)
}

### Tukey's bisquare (biweight) weight function ###
bisqWeightLS <- function(data, k) {
  w = (1 - (data/k)^2)^2
  w[abs(data) > k] <- 0 
  return(w)
}

### smoothing spline function (From RobRSVD package) ###
ssmatls <- function(n) {
  h = 1/(n + 1)
  Qmat = mat.or.vec(n, n - 1)
  for (j in 2:(n - 1)) {
    Qmat[j - 1, j] = 1/h
    Qmat[j, j] = -2/h
    Qmat[j + 1, j] = 1/h
  }
  Qmat = Qmat[, 2:(n - 1)]
  Rmat = mat.or.vec(n - 1, n - 1)
  for (i in 1:(n - 2)) {
    Rmat[i, i] = (2/3) * h
    Rmat[i, i + 1] = h/6
    Rmat[i + 1, i] = h/6
  }
  Rmat[n - 1, n - 1] = (2/3) * h
  Rmat = Rmat[2:(n - 1), 2:(n - 1)]
  y = Qmat %*% solve(Rmat) %*% t(Qmat)
  return(list(y = y, Qmat = Qmat, Rmat = Rmat, h = h))
}

### draw 2d heatmap using ggplot2 
draw_2dmap <- function(mat, colscale, colimit=NULL) {
  nr = dim(mat)[1]
  nc = dim(mat)[2]
  
  df = expand.grid(X=1:nr, Y=1:nc)
  df$Z = as.numeric(mat)
  if (is.null(colimit)) {
    ggobj = ggplot(df) +
      geom_tile(aes(X, Y, fill=Z)) +
      scale_fill_gradientn(colors=colscale)
  } else {
    ggobj = ggplot(df) +
      geom_tile(aes(X, Y, fill=Z)) +
      scale_fill_gradientn(colors=colscale, limits=colimit)
  }
  
  ggobj
}

# Recursion for difference operator matrix
# taken from fpca.ssvd function in refund
makeDiffOper <- function(degree, dim){
  if(degree==0){
    return(diag(dim))
  } else {
    return(diff(makeDiffOper(degree-1, dim)))
  }
}

# Adjusting sign 
sign_adjust <- function(truefun, estimfun) {
  res = estimfun
  if (sum(truefun*estimfun) < 0) {
    res = (-1)*estimfun
  }
  return(res)
}

compare_plot <- function(t, truefun, estimfun, title="") {
  
  truefun = truefun/normVector(truefun)
  
  mean.h1 = mean(truefun) - min(truefun)
  mean.h2 = max(truefun) - mean(truefun) 
  
  plot(t, truefun, type="l", xlab="", ylab="", main=title,
       ylim=c(min(truefun)-mean.h1, max(truefun)+mean.h2), col="gray")
  
  if (sum(truefun * estimfun) < 0) {
    lines(t, -estimfun)
  } else {
    lines(t, estimfun)
  }
}

## (simulation) evaluation functions for obtaining mean and standard deviation
mean_sd <- function(x) {
  mu_x = signif(mean(x), 3)
  sd_x = signif(sd(x), 3)
  
  res = paste0(mu_x , " (", sd_x, ")")
  
  return(res)
}


