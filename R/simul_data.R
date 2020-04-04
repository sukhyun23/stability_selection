rmultinorm <- function(n, p, rho) {
  sigma <- matrix(0, p, p)
  diag(sigma) <- 1
  
  # sigma[upper.tri(sigma)] <- runif(sum(upper.tri(sigma)), rho, 1)
  # rhos <- rep(rho, sum(upper.tri(sigma)))
  # rhos <- rhos*sample(c(1,-1), length(rhos), replace = T)
  
  sigma[upper.tri(sigma)] <- rho
  
  ind <- lower.tri(sigma)
  sigma[ind] <- t(sigma)[ind] 
  
  suppressWarnings(xx <- mvtnorm::rmvnorm(n, rep(0,p), sigma))
  
  n_minus <- round(p/2)
  idx_minus <- sort(sample(p, n_minus))
  
  xx[, idx_minus] <- -xx[, idx_minus]
  xx
}

x_generator <- function(
  n = 100, p = 1000, block = 20, seed = 1, cor_within_block = 0.9
) {
  n_block <- p/block  
  
  set.seed(seed)
  x_block <- replicate(
    block, rmultinorm(n, n_block, cor_within_block), simplify = F
  )
  
  x <- do.call(cbind, x_block)
  x <- apply(x, 2, function(x) (x-mean(x))/sd(x))
  x
}
y_generator <- function(
  x, noise_signal = 0.5, b, n_true = 50, idx_true, seed = 1
) {
  if (missing(b)) {
    b <- round(runif(n_true, min = 1.5, max = 3), 1)  
  }
  
  y <- b %*% t(x[, idx_true])
  set.seed(1)
  y <- c((y - mean(y))/sd(y) + rnorm(n, 0, noise_signal))
  
  y
}
