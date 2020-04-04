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
