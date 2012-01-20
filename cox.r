# Chris DuBois
# January 10th, 2012

# Helper functions for simulating from and fitting a Cox model.  MLE, MAP, and MH supported.
# C++ backend: see cox.cpp.r and cox.cpp

# Explanation of parameters:
# y: vector of observed choices.  Element i is in {1,...,|omega.i|}.
# X: a list where each element is a |omega.i| x P matrix whose jth row contains the covariates for the jth option.
# beta: a vector of parameters of length P

# Simulate from the model.  
cox.sim <- function(X,beta) {
  sapply(1:length(X),function(i) {
    py <- exp(X[[i]] %*% beta)
    sample(1:length(py),1,prob=py)
  })
}

# Compute the likelihood
cox.llk <- function(value,grad=FALSE) {
  y <- value$y
  beta <- value$pv
  llk <- cox.llk.cpp(y,X,beta)
  if (grad) {
    g <- cox.grad.cpp(y,X,beta)
    attr(llk,"grad") <- list(pv=as.vector(g))
  }
  return(llk)
}

# Compute the log prior
cox.lprior <- function(value,grad=FALSE) {
  sigma <- value$sigma
  mu <- value$mu
  pv <- value$pv
  lprior <- sum(dnorm(pv,value$mu,value$sigma, log=TRUE))
  if (grad) {
    attr(lprior,"grad") <- list(pv = - (pv - value$mu)/value$sigma^2)
  }
  return(lprior)
}

# Compute the log posterior
cox.lpost <- function(value,grad=FALSE){
  llk <- cox.llk(value,grad)
  lprior <- cox.lprior(value,grad)
  lpost <- llk+lprior
  if (grad)
    attr(lpost,"grad") <- list(pv= attr(llk,"grad")$pv + attr(lprior,"grad")$pv)#,mu=0,sigma=0) 
  return(lpost)
}

# Compute the log posterior (but different interface)
cox.lpost2 <- function(y,X,beta,mu=0,sigma=5,grad=FALSE){
  cox.llk.cpp(y,X,beta) + sum(dnorm(beta,mu,sigma, log=TRUE))
}

# Compute the MLE by optimizing the likelihood
cox.mle <- function(y,X) {
  npar <- ncol(X[[1]])
  beta <- rep(0,npar)
  fn <- function(pv) {
    - cox.llk.cpp(y,X,pv) 
  }
  gr <- function(pv) {
    - cox.grad.cpp(y,X,pv) 
  }
  optim(beta,fn,gr)$par
}


# Find maximum a posteriori (MAP) estimates.
# hyper: list for Normal hyperprior parameters, mu and sigma.
cox.map <- function(y,X,init=NULL,hyper=list(mu=0,sigma=5)) {
  npar <- ncol(X[[1]])
  if (is.null(init)) beta <- rnorm(npar,0,1)#rep(0,npar)
  else beta <- init
  fn <- function(pv) {
    - cox.llk.cpp(y,X,pv) - sum(dnorm(pv,hyper$mu,hyper$sigma, log=TRUE))
  }
  gr <- function(pv) {
    - (cox.grad.cpp(y,X,pv) - (pv - hyper$mu)/hyper$sigma^2)
  }
  optim(beta,fn,gr)$par
}

# Run Metropolis-Hastings
# ndraw:   number of draws
# burnin:  number of draws to discard as burnin
# mcmc.sd: standard deviation to use in proposal distribution
cox.mh <- function (y,X,hyper=list(mu=0,sigma=5),ndraw=1000,burnin=500,mcmc.sd=.1,verbose=TRUE) {
  npar <- ncol(X[[1]])
  beta.draws <- matrix(NA,ndraw,npar)
  params <- hyper
  params$pv <- rnorm(npar)
  params$y <- y
  olp <- cox.lpost(params)
  for (i in 1:ndraw) {
   for (j in 1:npar) {
    cand <- params
    cand$pv[j]  <- cand$pv[j] + rnorm(1,0,mcmc.sd)
    #clp <- cox.lpost2(y,X,cand$pv)
    clp <- cox.lpost(cand)
    if (clp - olp > log(runif(1))) {
      params <- cand
      olp <- clp
    }
   }
   beta.draws[i,] <- params$pv
   if (verbose) cat("iter",i,"vals",round(params$pv[1:5],3),"\n")
  }
  return(beta.draws)
}
