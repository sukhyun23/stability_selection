asd



library(dplyr)
library(data.table)
library(corrplot)
library(glmnet)
library(mvtnorm)
library(caret)
library(stabs)



p <- 10
x <- sy(100, 1000) 






# generate simultation data
synthetic <- function(x, n, p) {
  replicate(p, (sample(c(2,-2), 1)*x + rnorm(n, 0, runif(1, 0.1, 0.4))))
}

n <- 100
p <- 1000

p_true <- 5

set.seed(1)
x_true <- replicate(p_true, rnorm(n))
x_list <- list()
for (i in 1:p_true) {
  result <- synthetic(x_true[ ,i], n, 20)
  x_list[[i]] <- result
}
x_corr <- do.call(cbind, x_list)
x_list[[2]] %>% data.frame() %>% plot()

x <- cbind(x_true, x_corr)
# x_noise <- replicate(p-ncol(x), rnorm(n))
x_noise <- sy(n, p-ncol(x))

x <- cbind(x, x_noise)
x <- as.matrix(scale(x))


# corrplot(cor(result), method = 'shade', type ='lower', addCoef.col = 'black')
b <- c(
  -1.7, 1.9, 1.8, 2.1, 2.5, 2.1, -1.7, 1.9, 1.8, 2.1,
  -1.7, 1.9, 1.8, 2.1, 2.5, 2.1, -1.7, 1.9, 1.8, 2.1
)
b <- c(
  -1.7, 1.9, 1.8, 2.1, 2.5
)

idx_true <- 1:p_true

# set.seed(1)
# x <- replicate(p, rnorm(n))
# x <- as.matrix(scale(x))

x <- sy(n, p)
x %>% dim()

y <- b %*% t(x[, idx_true])
y <- c((y - mean(y))/sd(y) + rnorm(n, 0, 0.4))

lm(y ~ x[,1:5]) %>% summary()

# viz
# plot(data.frame(cbind(x[, idx_true], y)), pch = 19)
# mat_cor <- cor(cbind(x[, idx_true], y))
# corrplot(mat_cor, method = 'shade', type ='lower', addCoef.col = 'black')


# corrplot(cor(x[c()]), method = 'shade', type ='lower', addCoef.col = 'black')

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

object <- selection
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


model <- stability_selection(x, y, 0.2, 100)
model %>% plot()


ratio <- apply(selected, 1, sum)/(length(lambdas)*n_iter*2)

ratio10 <- tail(sort(ratio), 10)
idx_col <- (as.numeric(gsub('V', '', names(ratio10))) %in% idx_true) + 1
barplot(ratio10, col = idx_col, ylim = c(0,1))



lasso <- glmnet(x, y)
coef_lasso <- as.matrix(coef(lasso))[-1,]
plot(
  -lasso$lambda, 
  rep(NA, length(lasso$lambda)), 
  type='n', ylim = c(-0.4,0.4), ylab = 'ratio', xlab = '-lambda'
)
for (i in 1:p) {
  line_col <- 'grey50'
  if (i %in% idx_true) line_col <- 'red'
  
  line_size <- 1
  if (i %in% idx_true) line_size <- 2
  
  line_y <- coef_lasso[i,]
  lines(-lasso$lambda, line_y, col=line_col, lwd=line_size)
}





apply(selected, 1, sum)[idx_true]
apply(selected, 1, length)[idx_true]


apply(nonzero_lasso, 1, sum)[-idx_true]
apply(nonzero_lasso, 1, sum) %>% sort(decreasing = T)

nonzero_lasso


idx_lambda <- which(model_lasso$lambda == model_lasso$lambda.1se)

nbootstrap = 200
nsteps = 20
alpha = 0.2

dimx <- dim(x)
n <- dimx[1]
p <- dimx[2]
halfsize <- as.integer(n/2)
freq <- matrix(0,1,p)

for (i in seq(nbootstrap)) {
  
  # Randomly reweight each variable
  xs <- t(t(x)*runif(p,alpha,1))
  
  # Ramdomly split the sample in two sets
  perm <- sample(dimx[1])
  i1 <- perm[1:halfsize]
  i2 <- perm[(halfsize+1):n]
  
  # run the randomized lasso on each sample and check which variables are selected
  cv_lasso <- lars::cv.lars(xs[i1,],y[i1],plot.it=FALSE, mode = 'step')
  idx <- which.max(cv_lasso$cv - cv_lasso$cv.error <= min(cv_lasso$cv))
  coef.lasso <- coef(lars::lars(xs[i1,],y[i1]))[idx,]
  freq <- freq + abs(sign(coef.lasso))
  
  cv_lasso <- lars::cv.lars(xs[i2,],y[i2],plot.it=FALSE, mode = 'step')
  idx <- which.max(cv_lasso$cv - cv_lasso$cv.error <= min(cv_lasso$cv))
  coef.lasso <- coef(lars::lars(xs[i1,],y[i1]))[idx,]
  freq <- freq + abs(sign(coef.lasso))
  print(freq)
}

glmnet
x
xs <- t(t(x)*runif(p,0.2,1))


