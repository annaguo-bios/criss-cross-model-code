###############################################################################
#########################case 2: varying rho ##################################
###############################################################################

####################################################
#####Non-optimal GMM assuming known intercept====
####################################################
rm(list=ls())
# packages
library(stats4)
library(cubature)
library(R.utils)
library(reshape2)
library(ggplot2)
library(BB)  #BB: Solving and Optimizing Large-Scale Nonlinear Systems
library(gmm) #Computing Generalized Method of Moments and Generalized Empirical Likelihood with R
library(dplyr)
library(tidyr)

# define basic functions
expit <- function(x){
  return(1/(1+exp(-x)))
}


# number of simulation runs
t <- 100

# sample size
n=1000

# sample size
rho.vec <- c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)

# set up the parameters
mu.x <- 0.4
mu.y <- 2

sig.x <- 3
sig.y <- 1


# missing indicator

b0.rx <- -0.5
b1.rx <- 1

b0.ry <- 2
b1.ry <- -1
b2.ry <- 0.7


############define the estimating equation###############
f.est1 <- function(para,x) { #para: parameters to be estimated: alpha, beta; x: data
  
  # data
  X = x[,1]
  Y = x[,2]
  Rx = x[,3]
  Ry = x[,4]
  
  # parameters to be estimated
  alpha <- para[1]
  beta <- para[2]
  
  # fit logistic regression to get p(Ry=1|Rx=1, X)
  cond.ry = glm(Ry ~ Rx+X, data = x, family = "binomial")
  
  E.ry = rep(NA,n)
  E.ry[Rx==1] = predict(cond.ry,type="response")
  
  # E(x|y)
  
  E.x1y <- alpha+beta*Y
  
  # returned estimating function
  c1 = (Rx*Ry/E.ry)*(X-E.x1y)
  
  f = cbind(c1)
  
  # f=0 for Rx=0 or Ry=0 which is where the NA's located
  f[is.na(f)] <- 0
  return(f)
}


# empty data frame to hold the parameter estimate
columns = c("rho","beta")
estimated = data.frame(matrix(nrow = t*length(rho.vec), ncol = length(columns))) 

colnames(estimated) = columns

# empty list to hold the covariance estimate
covariance <- vector("list", t*length(rho.vec)) 


set.seed(7)
for (rho in rho.vec){
  for (i in 1:t){
    # generate data
    X = rnorm(n,mu.x,sig.x)
    Y = rnorm(n,mu.y+rho*sig.y/sig.x*(X-mu.x),sqrt((1-rho^2))*sig.y)
    Rx = rbinom(n,1,expit(b0.rx+b1.rx*Y))
    Ry = rbinom(n,1,expit(b0.ry+b1.ry*Rx+b2.ry*X))
    X[Rx==0] <- NA
    Y[Ry==0] <- NA
    
    dt <- data.frame(X=X, Y=Y, Rx=Rx, Ry=Ry)
    
    # estimation
    alpha=mu.x-rho*sig.x/sig.y*mu.y
    
    # be sure to set the upper and lower bound that covers the true value
    ans <- gmm(f.est1,x=dt,t0=c(alpha,1),eqConst=c(1),method="Brent",lower=-5, upper=5)
    
    index <- (which(rho.vec==rho)-1)*t+i
    
    estimated[index, ] <- c(rho, ans$coefficients)
    covariance[[index]] <- ans$vcov
  }
}

# save data
save(list = c("estimated","covariance"),file = "./Known intercept/local_rho_GMM_bivariate_dat.Rdata")

###########################################
#####Optimal GMM assuming known intercept====
###########################################
rm(list=ls())
# packages
library(stats4)
library(cubature)
library(R.utils)
library(reshape2)
library(ggplot2)
library(BB)  #BB: Solving and Optimizing Large-Scale Nonlinear Systems
library(gmm) #Computing Generalized Method of Moments and Generalized Empirical Likelihood with R
library(dplyr)
library(tidyr)

# define basic functions
expit <- function(x){
  return(1/(1+exp(-x)))
}


# number of simulation runs
t <- 100

# sample size
n=1000

# sample size
rho.vec <- c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)

# set up the parameters
mu.x <- 0.4
mu.y <- 2

sig.x <- 3
sig.y <- 1


# missing indicator

b0.rx <- -0.5
b1.rx <- 1

b0.ry <- 2
b1.ry <- -1
b2.ry <- 0.7


# use the estimated values obtained as medians over 100 simulation runs from the non-optimal GMM for constructing optimal GMM
load("./Known intercept/local_rho_GMM_bivariate_dat.Rdata")
beta.est=estimated %>% group_by(rho) %>% summarise(est=median(beta)) %>% pull(est)

############define the estimating equation###############
f.est1 <- function(para,x) { #para: parameters to be estimated: alpha, beta; x: data
  
  # data
  X = x[,1]
  Y = x[,2]
  Rx = x[,3]
  Ry = x[,4]
  
  # parameters to be estimated
  alpha <- para[1]
  beta <- para[2]
  
  # fit logistic regression to get p(Ry=1|Rx=1, X)
  cond.ry = glm(Ry ~ Rx+X, data = x, family = "binomial")
  
  E.ry = rep(NA,n)
  E.ry[Rx==1] = predict(cond.ry,type="response")
  
  # E(x|y)
  
  E.x1y <- alpha+beta*Y
  
  # f(Y)
  E <- alpha0+beta0*Y
  V <- (1-rho^2)*sig.x^2
  pi0 <- coef(cond.ry)[names(coef(cond.ry))=="(Intercept)"]
  pi1 <- coef(cond.ry)[names(coef(cond.ry))=="X"]
  fy2 <- Y/(V+exp(-0.5/V*(E^2+2*V*pi0-(E-V*pi1)^2))*(V+V^2*pi1^2))
  
  # returned estimating function
  c2 = (Rx*Ry/E.ry)*fy2*(X-E.x1y)
  
  f = cbind(c2)
  
  
  # f=0 for Rx=0 or Ry=0 which is where the NA's located
  f[is.na(f)] <- 0
  return(f)
}

# empty data frame to hold the parameter estimate
columns = c("rho","beta")
estimated = data.frame(matrix(nrow = t*length(rho.vec), ncol = length(columns))) 

colnames(estimated) = columns

# empty list to hold the covariance estimate
covariance <- vector("list", t*length(rho.vec)) 


set.seed(6)
for (rho in rho.vec){
  alpha0 = mu.x-rho*sig.x/sig.y*mu.y
  beta0 = beta.est[which(rho.vec==rho)]
  for (i in 1:t){
    # generate data
    X = rnorm(n,mu.x,sig.x)
    Y = rnorm(n,mu.y+rho*sig.y/sig.x*(X-mu.x),sqrt((1-rho^2))*sig.y)
    Rx = rbinom(n,1,expit(b0.rx+b1.rx*Y))
    Ry = rbinom(n,1,expit(b0.ry+b1.ry*Rx+b2.ry*X))
    X[Rx==0] <- NA
    Y[Ry==0] <- NA
    
    dt <- data.frame(X=X, Y=Y, Rx=Rx, Ry=Ry)
    
    # estimation
    # be sure to set the upper and lower bound that covers the true value
    ans <- gmm(f.est1,x=dt,t0=c(alpha0,1),eqConst=c(1),method="Brent",lower=-5, upper=5) 
    
    index <- (which(rho.vec==rho)-1)*t+i
    
    estimated[index, ] <- c(rho, ans$coefficients)
    covariance[[index]] <- ans$vcov
  }
}

# save data
save(list = c("estimated","covariance"),file = "./Known intercept/local_rho_GMM_opt_est_bivariate_dat.Rdata")

##############################
#####Pseudo-likelihood====
# better to run on cluster
##############################
rm(list=ls())
# packages
library(stats4)
library(cubature)
library(R.utils)
library(reshape2)
library(ggplot2)
library(BB)  #BB: Solving and Optimizing Large-Scale Nonlinear Systems
library(gmm) #Computing Generalized Method of Moments and Generalized Empirical Likelihood with R
library(dplyr)
library(tidyr)

# define basic functions
expit <- function(x){
  return(1/(1+exp(-x)))
}

# n of simulations
t <- 100

# sample size
n=1000

# sample size
rho.vec <- c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)

# set up the parameters
mu.x <- 0.4
mu.y <- 2

sig.x <- 3
sig.y <- 1


# missing indicator

b0.rx <- -0.5
b1.rx <- 1

b0.ry <- 2
b1.ry <- -1
b2.ry <- 0.7

# empty data frame to hold the parameter estimate
columns = c("rho","OR")
estimated = data.frame(matrix(nrow = t*length(rho.vec), ncol = length(columns))) 

colnames(estimated) = columns

# empty list to hold the covariance estimate
covariance <- vector("list", t*length(rho.vec)) 

set.seed(7)

for (rho in rho.vec){
  for (i in 1:t){
    # generate data
    X = rnorm(n,mu.x,sig.x)
    Y = rnorm(n,mu.y+rho*sig.y/sig.x*(X-mu.x),sqrt((1-rho^2))*sig.y)
    Rx = rbinom(n,1,expit(b0.rx+b1.rx*Y))
    Ry = rbinom(n,1,expit(b0.ry+b1.ry*Rx+b2.ry*X))
    X[Rx==0] <- NA
    Y[Ry==0] <- NA
    
    dt <- data.frame(X=X, Y=Y, Rx=Rx, Ry=Ry)
    dt.obs <- na.omit(dt)
    
    # get combination
    comb <- combn(nrow(dt.obs),2)
    
    # data  for logistics regression
    u <- vector(mode = "numeric", length = ncol(comb)) #outcome
    v <- vector(mode = "numeric", length = ncol(comb)) #covariate
    for (j in 1:ncol(comb)){
      # the smaller index
      s <- comb[1,j]
      
      # the larger index
      l <- comb[2,j]
      
      v[j] <- with(dt.obs, (X[s]-X[l])*abs(Y[s]-Y[l]))
      u[j] <- with(dt.obs, ifelse((Y[s]-Y[l])>0,1,0))
    }
    
    log.dt <- data.frame(u=u,v=v)
    
    # Logistic regression without intercept
    m1 <- glm(u ~ v - 1, family = binomial, data = log.dt)
    
    index <- (which(rho.vec==rho)-1)*t+i
    estimated[index, ] <- c(rho, summary(m1)$coefficients[1])
    covariance[[index]] <- summary(m1)$coefficients[2]
  }
}

# save data
save(list = c("estimated","covariance"),file = "local_rho_order_statistic_bivariate_dat.Rdata")