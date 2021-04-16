# KMLMM Project: RVM Study
# Jake Watson, 16/01/2021
# UPC

# Return sigmoid of input
sigmoid <- function(X) {
  return(1.0/(1.0 + exp(-X)))
}

# Return euclidean distances of rows between 2 matrices
distsq <- function(X, Y) {
  library(pracma)
  return(distmat(X, Y)^2)
}

# Set of kernel functions to return inner product between 2 matrices
# Available kernels:
# linear distance
# gaussian
# radial basis function (rbf)
Kernel <- function(X, Y, kernel, sigma = 0.1) {
  library(kernlab)
  library(pracma)
  
  N <- nrow(X)
  d <- ncol(Y)
  
  if (strcmp(kernel, "distance")) {
    return(sqrt(distsq(X,Y)))
  } 
  
  if (strcmp(kernel, "gauss")) {
    return(exp(-distsq(X, Y)))
  }
  
  if (strcmp(kernel, "radial")) {
    return(exp(-distsq(X, Y))*sigma)
  }
  
  else {
    stop("Error: invalid kernel function!")
  }
}

# Wrapper for main RVM estimator function
# converts inputs to matrices
# adds bias columns, and strips bias at the end
#
# OUTPUTS
#   vectors      indices of the relevance vectors in the training set
#   weights      weights of the relevance vectors
#   bias         scalar bias value
#
# INPUTS
#   X            data representation in original form
#   t            target values
#   mode         string, either (REGRESSION | CLASSIFICATION)
#   kernel_type  string, kernel type (distance | gauss | radial)
#   max_it       maximum iterations
#
RVM <- function(X, t, mode, kernel_type, max_it) {
  X <- as.matrix(X)
  t <- as.matrix(t)
  
  N <- nrow(X)
  d <- ncol(X)
  
  # Convert to kernelized basis
  PHI	 <-  Kernel(X, X, kernel_type)
  # add bias vector
  PHI	 <-  cbind(PHI, matrix(1L, nrow=N, ncol=1))
  
  # initial hyperparameter guess
  alpha <- (1/N)^2
  
  # obtain relevance vectors
  result <- RVM_train(PHI, t, mode, alpha, kernel, max_it)
  weights <- result$weights
  vectors <- result$vectors
  
  # strip bias
  bias <- 0
  indexBias	= N + 1
  vector_index = which(vectors == indexBias)
  
  if (length(vector_index) != 0) {
    bias <- weights[vector_index]
    vectors <- vectors[-vector_index]
    weights <- weights[-vector_index]
  }
  
  return(list("vectors" = vectors, "weights" = weights, "bias" = bias))
}

# Estimates the relevance vectors and their weights
#
# OUTPUTS
#   vectors      indices of the relevance vectors in the training set
#   weights      weights of the relevance vectors
#
# INPUTS
#   phi     data representation in kernel basis form
#   t       target values
#   mode    string, either (REGRESSION | CLASSIFICATION)
#   kernel  string, kernel type (distance | gauss | radial)
#   alpha   initial hyperparameters
#   max_it   maximum iterations
#
RVM_train <- function(PHI, t, mode, alpha, kernel, max_it) {
  
  N <- nrow(PHI)
  M <- ncol(PHI)
  
  w <- matrix(0L, nrow = M, ncol = 1)
  alpha <- matrix(alpha, nrow = M, ncol = 1)
  gamma <- matrix(1L, nrow = M, ncol = 1)
  PHIt <- t(PHI) %*% t
  
  # convergence parameter
  CONVERGENCE = 1E-3
  
  # maximum alpha size 
  # if alpha exceeds for a vector,
  # prune that vector
  PRUNING_THRESHOLD = 1e9
  
  LAST_IT <- FALSE
  REGRESSION <- strcmp(mode, "REGRESSION")
  
  epsilon <- std(t) * 10/100
  beta <- 1/epsilon^2
  
  for (i in 1:max_it) {
    
    # do pruning
    useful_indices <- which(alpha < PRUNING_THRESHOLD)
    useless_indices <- which(alpha >= PRUNING_THRESHOLD)
    useful <- alpha[useful_indices]
    M <- length(useful_indices)
    w[useless_indices] <- 0L
    PHI_used <- PHI[, useful_indices]
    
    if (REGRESSION){
      # determine gaussian likelihood
      Hessian <- (t(PHI_used) %*% PHI_used) * beta + diag(useful)
      U <- chol(Hessian)
      Ui <- solve(U)
      w[useful_indices] <- (Ui %*% (t(Ui) %*% PHIt[useful_indices]))*beta
      ED <- sum(t - (PHI_used %*% w[useful_indices])^2)
      dataLikelihood <- (N*log(beta) - beta*ED)/2
    } else {
      # for classification we use a bernoulli likelihood
      # we call the posterior mode finder
      post <- posterior(w[useful_indices], alpha[useful_indices], PHI_used, t, 75)
      w[useful_indices] <- post$weight
      Ui <- post$Ui
      dataLikelihood <- post$lMode
    }
    
    diagSig <- rowSums(Ui^2)
    gamma <- 1 - useful*diagSig
    
    if (!LAST_IT) {
      
      # MacKay updates of hyperparameters
      oldAlpha <- log(alpha[useful_indices])
      alpha[useful_indices] <- gamma / (w[useful_indices]^2)
      au <- alpha[useful_indices]
      
      # check for convergence of log hyperparameters
      maxDAlpha <- max(abs(oldAlpha[which(au != 0)]) - log(au[which(au != 0)]))
      
      if (maxDAlpha < CONVERGENCE) {
        print(paste("Reached convergence in ", i, "iterations"))
        LAST_IT <- TRUE
      }
      
      if (REGRESSION) {
        beta = (N - sum(gamma)/ED)
      }
      
    } else {
      print("Max iterations reached")
      break
    }
  }
  
  weights <- w[useful_indices]
  vectors <- useful_indices
  
  return(list("vectors" = vectors, "weights" = weights))
}

# Finds the mode of the posterior distribution
#
# OUTPUTS
#   w       weights for each vector, maximises the posterior
#   Ui      inverse of the cholesky factor of the hessian
#   lMode   log likelihood of the data at the mode
#
# INPUTS
#   phi     data representation in kernel basis form
#   t       target values
#   w       initial vector weights
#   alpha   initial hyperparameters
#   iters   maximum iterations
#
posterior <- function(w, alpha, phi, t, iters) {
  
  # convergence parameter
  GRAD_STOP <- 1e-6
  # minimum search step
  LAMBDA_MIN <- 2^(-8)
  
  N <- nrow(phi)
  d <- ncol(phi)
  M <- length(w)
  
  A <- diag(alpha)
  errs <- matrix(0L, nrow = iters, ncol = 1)
  PHIw <- phi %*% w
  y = sigmoid(as.matrix(PHIw))
  
  # initial log posterior value
  data_term <- -(sum(log(y[which(t == 1)])) + sum(log(1-y[which(t==0)])))/N
  regulariser <- (t(alpha) %*% (w^2))/(2*N)
  err_new <- data_term + regulariser
  
  for (i in 1:iters) {
    yvar <- y * (1 - y)
    PHIV <- phi * (yvar %*% matrix(1L, nrow = 1, ncol = d))
    e <- (t - y)
    
    # initial gradient vector and hessian
    g <- (t(phi) %*% e) - (alpha * w)
    Hessian <- (t(PHIV) %*% phi) + A
    
    if (i==1) {
      if (rcond(Hessian) < (2^-52)){
        stop("ill conditioned Hessian")
      }
    }
    
    errs[i] <- err_new
    
    # check for convergence
    if (i > 2 & norm(g)/M < GRAD_STOP) {
      errs <- errs[1:i]
      print(paste("posterior converged after ", i, " iterations"))
      break
    }
    
    # Newton step
    
    U <- chol(Hessian)
    delta_w <- mldivide(U, mldivide(t(U), g))
    lambda = 1
    
    while (lambda > LAMBDA_MIN) {
      w_new <- w + lambda*delta_w
      PHIw <- phi %*% w_new
      y <- sigmoid(PHIw)
      
      # compute new error
      if (any(y[which(t == 1)] == 0) | any(y[which(t == 0)] == 1)){
        err_new = Inf
      } else {
        data_term <- -(sum(log(y[which(t == 1)])) + sum(log(1 - y[which(t == 0)])))/N
        regulariser <- (t(alpha) %*% (w_new^2))/(2*N)
        err_new <- data_term + regulariser
      }
      
      # if the error has increased, reduce the step size
      if (err_new > errs[i]) {
        lambda = lambda/2
      } else {
      # error decreased, accept the step
        w = w_new
        lambda = 0
      }
    }
    
    # we have converged close to the optimum and cannot take any more steps
    if (lambda) {
      break
    }
  }
  
  # calculate final values and return
  Ui <- solve(U)
  lMode <- -N*data_term
  
  return(list("weight" = w, "Ui" = Ui, "lMode" = lMode))
}

