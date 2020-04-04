asd



rmultinorm <- function(n, p, rho) {
  sigma <- matrix(0, p, p)
  diag(sigma) <- 1
  
  # sigma[upper.tri(sigma)] <- runif(sum(upper.tri(sigma)), rho, 1)
  sigma[upper.tri(sigma)] <- rep(rho, sum(upper.tri(sigma)))
  
  
  ind <- lower.tri(sigma)
  sigma[ind] <- t(sigma)[ind] 
  
  suppressWarnings(xx <- mvtnorm::rmvnorm(n, rep(0,p), sigma))
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
y_generator <- function(x, noise_signal = 0.5, idx_true, seed = 1) {
  b <- c(
    1.6, 2.9, 1.9, 2.9, 2.8, 2, 2.4, 2.3, 1.5, 2.9, 
    2.4, 2.2, 2, 2.4, 2.3, 2.4, 2.8, 2.2, 2.5, 2.2, 
    2.5, 2, 2.4, 2.4, 3, 2.6, 2.2, 2.6, 2, 2.5, 2.9,
    1.7, 1.5, 1.9, 2.1, 1.9, 1.9, 2.4, 2.9, 1.6, 1.8,
    2.4, 2.7, 2.4, 1.9, 3, 2.8, 1.5, 2, 1.6
  )  
  
  y <- b %*% t(x[, idx_true])
  set.seed(1)
  y <- c((y - mean(y))/sd(y) + rnorm(n, 0, noise_signal))
  
  y
}

n <- 100
p <- 1000
block <- 20
x <- x_generator(n = n, p = p, block = block)

idx_true <- seq(1, p, by = block)
y <- y_generator(x, 1.5, idx_true)

m <- lm(y ~ x[, idx_true])
summary(m)

lasso <- glmnet(x, y)
plot(lasso)

coef_lasso <- as.matrix(coef(lasso))[-1,]
plot(
  -lasso$lambda, 
  rep(NA, length(lasso$lambda)), 
  type='n', ylim = c(min(coef_lasso),max(coef_lasso)), ylab = 'ratio', xlab = '-lambda'
)
for (i in 1:p) {
  line_col <- 'grey50'
  if (i %in% idx_true) line_col <- 'red'
  
  line_size <- 1
  if (i %in% idx_true) line_size <- 2
  
  line_y <- coef_lasso[i,]
  lines(-lasso$lambda, line_y, col=line_col, lwd=line_size)
}


selection <- stability_selection(x, y, 0.2, 1000, idx_true)

plot(selection)


list.dirs()

knitr::knit('README.Rmd')

