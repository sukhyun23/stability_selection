asd

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

n <- 100
p <- 1000
block <- 20

x <- x_generator(
  n = n, p = p, block = block, seed = 1, cor_within_block = 0.75
)


idx_true <- seq(1, p, by = p/block) # block당 하나
b <- c(
  1.6, 2.9, 1.9, 2.9, 2.8, 2, 2.4, 2.3, 1.5, 2.9,
  2.4, 2.2, 2, 2.4, 2.3, 2.4, 2.8, 2.2, 2.5, 2.2
)
y <- y_generator(
  x = x, n_true = block, b = b, noise_signal = 0.5, 
  idx_true = idx_true, seed = 1
)

m <- lm(y ~ x[, idx_true])
summary(m)

lasso <- glmnet(x, y)
plot(lasso)

coef_lasso <- as.matrix(coef(lasso))[-1,]
plot(
  -lasso$lambda, 
  rep(NA, length(lasso$lambda)), 
  type='n', ylim = c(min(coef_lasso),max(coef_lasso)), ylab = 'coefficient',
  xlab = '-lambda'
)
for (i in 1:p) {
  line_col <- 'grey50'
  if (i %in% idx_true) line_col <- 'red'
  
  line_size <- 1
  if (i %in% idx_true) line_size <- 2
  
  line_y <- coef_lasso[i,]
  lines(-lasso$lambda, line_y, col=line_col, lwd=line_size)
}



stability_selection <- function(x, y, alpha, n_iter = 100, idx_true) {
  n <- nrow(x)
  p <- ncol(x)
  
  # lambda grid
  lambdas <- glmnet(x, y, alpha = 1, nlambda = 100)$lambda
  # lambdas <- sort(lambdas)
  selected <- matrix(0, nrow = p, ncol = length(lambdas))
  
  for (i in 1:n_iter) {
    # add perturbation(rescaling)
    x_rescaled <- apply(x, 2, function(x) x + runif(n, alpha, 1))
    
    # resampling index
    idx <- sample(n)
    idx1 <- idx[1:(n/2)]
    idx2 <- idx[-(1:(n/2))]
    
    # model fitting
    model1 <- glmnet(x_rescaled[idx1,], y[idx1], alpha = 1, lambda = lambdas)
    coef1 <- as.matrix(coef(model1))
    nonzero1 <- abs(sign(coef1))[-1,]
    
    model2 <- glmnet(x_rescaled[idx2,], y[idx2], alpha = 1, lambda = lambdas)
    coef2 <- as.matrix(coef(model2))
    nonzero2 <- abs(sign(coef2))[-1,]
    
    selected <- selected + nonzero1 + nonzero2
  }
  
  ratio_lambda <- apply(selected, 2, function(x) x/(n_iter*2))
  
  result <- list(
    ratio_lambda = ratio_lambda,
    lambdas = lambdas, 
    n = n, p = p, idx_true = idx_true
  )
  structure(result, class = 'stab')
}
plot.stab <- function(object) {
  plot(
    x = -object$lambdas, 
    y = rep(NA, length(object$lambdas)), 
    type='n', ylab = 'ratio', xlab = '-lambda',
    ylim = c(min(object$ratio_lambda), max(object$ratio_lambda))
  )
  
  for (i in 1:object$p) {
    line_col <- 'grey75'
    if (i %in% object$idx_true) line_col <- 'red'
    
    line_size <- 1
    if (i %in% object$idx_true) line_size <- 2
    
    line_y <- object$ratio_lambda[i,]
    line_y[which.min(line_y)] <- 0
    
    lines(-object$lambdas, line_y, col = line_col, lwd = line_size)
  }
}
selection <- stability_selection(x, y, 0.2, 1000, idx_true)
plot(selection)


knitr::knit('./README/README.Rmd')

