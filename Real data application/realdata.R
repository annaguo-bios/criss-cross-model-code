###############################################################################
##################### REAL DATA APPLICATION ##############################
###############################################################################

## MCRF====

##############################
#####Non-optimal GMM====
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
library(geepack)
library(dlookr) #EDA report
library(maxLik)
library(boot) # for bootstrap

# define basic functions
expit <- function(x){
  return(1/(1+exp(-x)))
}

dt <- muscatine %>% filter(occasion %in% c(1,3))

# convert long to wide
dt <- reshape(dt, idvar = c("id","gender","base_age"), timevar = "occasion", v.names=c("numobese","obese","age") , direction = "wide")
dt <- dt %>% mutate(m1=1-as.numeric(is.na(numobese.1)), m2=1-as.numeric(is.na(numobese.3))) %>% select(numobese.1,numobese.3,m1,m2,gender)

dt.girl <- dt %>% filter(gender=="F")
dt.boy <- dt %>% filter(gender=="M")


############define the estimating equation###############
f.est1 <- function(para,x) { #para: parameters to be estimated: theta01, theta10
  # data
  X = x[,1]
  Y = x[,2]
  Rx = x[,3]
  Ry = x[,4]
  
  # parameter that is known: theta00
  theta00
  
  # parameters to be estimated
  theta01 <- para[1]
  theta10 <- para[2]
  theta11 <- 1-theta00-theta10-theta01
  
  # p(x=1|y=1)
  p1.1 <- theta11/(theta01+theta11)
  
  # p(x=1|y=0)
  p1.0 <- theta10/(theta00+theta10)
  
  # fit logistic regression to get p(Ry=1|Rx=1, X)
  cond.ry = glm(Ry ~ Rx+X, data = x, family = "binomial")
  
  E.ry = rep(NA,nrow(x))
  E.ry[Rx==1] = predict(cond.ry,type="response")
  
  
  # E(x|y)
  E.x1y <- p1.0+(p1.1-p1.0)*Y
  
  # estimating function
  c1 = (Rx*Ry/E.ry)*(X-E.x1y)
  c2 = (Rx*Ry/E.ry)*Y*(X-E.x1y)
  
  f = cbind(c1,c2)
  
  # f=0 for Rx=0 or Ry=0 which is where the NA's located
  f[is.na(f)] <- 0
  return(f)
}

## bootstrap for the variance of log(OR)
set.seed(7)
f.boot.GMM <- function(d,i){
  d2 <- d[i,]
  ans <- gmm(f.est1,x=d2,t0=c(0.1,0.1))
  logOR <-log({theta00*(1-theta00-sum(ans$coefficients))}/(ans$coefficients[1]*ans$coefficients[2]))
}

# for boys
theta00 <- with(na.omit(dt.boy), mean(numobese.1==0 & numobese.3==0))
ans <- gmm(f.est1,x=dt.boy,t0=c(0.1,0.1))
output.boy <- data.frame(theta=c(theta00, ans$coefficients, 1-theta00-sum(ans$coefficients)), var=c(NA,diag(ans$vcov),sum(ans$vcov)))

# bootstrap
boot.GMM.boy <- boot(dt.boy, f.boot.GMM, R=1000)
logor.boy <- data.frame(logOR=with(output.boy,log(theta[1]*theta[4]/(theta[2]*theta[3]))), var=var(boot.GMM.boy$t))

# for girls
theta00 <-with(na.omit(dt.girl), mean(numobese.1==0 & numobese.3==0))
ans <- gmm(f.est1,x=dt.girl,t0=c(0.1,0.1))
output.girl <- data.frame(theta=c(theta00, ans$coefficients, 1-theta00-sum(ans$coefficients)), var=c(NA,diag(ans$vcov),sum(ans$vcov)))

# bootstrap
boot.GMM.girl <- boot(dt.girl, f.boot.GMM, R=1000)
logor.girl <- data.frame(logOR=with(output.girl,log(theta[1]*theta[4]/(theta[2]*theta[3]))), var=var(boot.GMM.girl$t))

save(list = c("output.boy","output.girl","logor.boy","logor.girl"),file = "realdata_GMM.Rdata")



##############################
#####Optimal GMM====
##############################

############define the estimating equation###############
f.est1 <- function(para,x) { #para: parameters to be estimated: alpha, beta; x: data
  
  # data
  X = x[,1]
  Y = x[,2]
  Rx = x[,3]
  Ry = x[,4]
  
  # parameter that is known: theta00
  theta00
  
  # parameters to be estimated
  theta01 <- para[1]
  theta10 <- para[2]
  theta11 <- 1-theta00-theta10-theta01
  
  # p(x=1|y=1)
  p1.1 <- theta11/(theta01+theta11)
  
  # p(x=1|y=0)
  p1.0 <- theta10/(theta00+theta10)
  
  # fit logistic regression to get p(Ry=1|Rx=1, X)
  cond.ry = glm(Ry ~ Rx+X, data = x, family = "binomial")
  
  E.ry = rep(NA,nrow(x))
  E.ry[Rx==1] = predict(cond.ry,type="response")
  
  # E(x|y)
  E.x1y <- p1.0+(p1.1-p1.0)*Y
  
  ## optimal f(Y) from Theorem 3
  # intercept and slope from regression Ry~Rx+X
  pi0 <- coef(cond.ry)[names(coef(cond.ry))=="(Intercept)"]
  pi1 <- coef(cond.ry)[names(coef(cond.ry))=="X"]
  
  # E{(x-h(y))^2/pi(x)|Y}
  E <- ((est.p1.0^2/pi0)*est.p0.0+{(1-est.p1.0)^2/(pi0+pi1)}*est.p1.0)*(1-Y)+
    ({(est.p1.1)^2/pi0}*est.p0.1+{(1-est.p1.1)^2/(pi0+pi1)}*est.p1.1)*Y
  
  # a(Y)
  aY <- data.frame(aY1=-theta[4]/(theta[2]+theta[4])^2*Y,
                   aY2=theta[1]/(theta[1]+theta[3])^2*(1-Y))
  
  fy1 <- E^(-1)*aY[,1]
  fy2 <- E^(-1)*aY[,2]
  
  # returned estimating function
  c1 = (Rx*Ry/E.ry)*fy1*(X-E.x1y)
  c2 = (Rx*Ry/E.ry)*fy2*(X-E.x1y)
  f = cbind(c1,c2)
  
  # f=0 for Rx=0 or Ry=0 which is where the NA's located
  f[is.na(f)] <- 0
  return(f)
}

## bootstrap for the variance of log(OR)
set.seed(7)
f.boot.GMM <- function(d,i){
  d2 <- d[i,]
  ans <- gmm(f.est1,x=d2,t0=c(0.1,0.1))
  logOR <-log({theta00*(1-theta00-sum(ans$coefficients))}/(ans$coefficients[1]*ans$coefficients[2]))
}

# use the estimated values obtained from the non-optimal GMM for constructing optimal GMM
load("realdata_GMM.Rdata")

## for boys
# estimates from non-optimal GMM
theta <- output.boy$theta # theta00, theta01, theta10, theta11
est.p0.0 <- theta[1]/(theta[1]+theta[3]) # p(x=0|y=0)
est.p0.1 <- theta[2]/(theta[2]+theta[4]) # p(x=0|y=1)
est.p1.0 <- theta[3]/(theta[1]+theta[3]) # p(x=1|y=0)
est.p1.1 <- theta[4]/(theta[2]+theta[4]) # p(x=1|y=1)

theta00 <- output.boy$theta[1]
ans <- gmm(f.est1,x=dt.boy,t0=c(0.1,0.1))
output.boy <- data.frame(theta=c(theta00, ans$coefficients, 1-theta00-sum(ans$coefficients)), var=c(NA,diag(ans$vcov),sum(ans$vcov)))

# bootstrap
boot.GMM.boy <- boot(dt.boy, f.boot.GMM, R=1000)
logor.boy <- data.frame(logOR=with(output.boy,log(theta[1]*theta[4]/(theta[2]*theta[3]))), var=var(boot.GMM.boy$t))

## for girls
# estimates from non-optimal GMM
theta <- output.girl$theta # theta00, theta01, theta10, theta11
est.p0.0 <- theta[1]/(theta[1]+theta[3]) # p(x=0|y=0)
est.p0.1 <- theta[2]/(theta[2]+theta[4]) # p(x=0|y=1)
est.p1.0 <- theta[3]/(theta[1]+theta[3]) # p(x=1|y=0)
est.p1.1 <- theta[4]/(theta[2]+theta[4]) # p(x=1|y=1)

theta00 <- output.girl$theta[1]
ans <- gmm(f.est1,x=dt.girl,t0=c(0.1,0.1))
output.girl <- data.frame(theta=c(theta00, ans$coefficients, 1-theta00-sum(ans$coefficients)), var=c(NA,diag(ans$vcov),sum(ans$vcov)))

# bootstrap
boot.GMM.girl <- boot(dt.girl, f.boot.GMM, R=1000)
logor.girl <- data.frame(logOR=with(output.girl,log(theta[1]*theta[4]/(theta[2]*theta[3]))), var=var(boot.GMM.girl$t))

save(list = c("output.boy","output.girl","logor.boy","logor.girl"),file = "realdata_optimal_GMM.Rdata")


##############################
#####Pseudo-likelihood====
# better to run on cluster
##############################

set.seed(7)
pselik <- function(x){
# boys
X = x[,1]
Y = x[,2]
Rx = x[,3]
Ry = x[,4]

X[Rx==0] <- NA
Y[Ry==0] <- NA

dt <- data.frame(X,Y,Rx,Ry)

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

estimated <- summary(m1)$coefficients[1]
covariance <- summary(m1)$coefficients[2]    

return(c(estimated, covariance))
}    

# boys
ans <- pselik(dt.boy)
logor.boy <- data.frame(logor=ans[1],sd=ans[2])

# girls
ans <- pselik(dt.girl)
logor.girl <- data.frame(logor=ans[1],sd=ans[2])

# save data
save(list = c("logor.boy","logor.girl"),file = "realdata_order_statistic.Rdata")


### make a table
rm(list=ls())
library(huxtable)
library(dplyr)
library(tidyr)
library(tidyverse)


# GMM
load("realdata_GMM.Rdata")
gmm.boy <- data.frame(est=c(output.boy$theta,logor.boy$logOR), sd=sqrt(c(output.boy$var, logor.boy$var)))%>% 
  mutate_all(funs(round(.,3))) %>%
  mutate(est.var=paste0(est," (",sd,")")) %>% select(est.var)

gmm.girl <- data.frame(est=c(output.girl$theta,logor.girl$logOR), sd=sqrt(c(output.girl$var, logor.girl$var)))%>% 
  mutate_all(funs(round(.,3))) %>%
  mutate(est.var=paste0(est," (",sd,")")) %>% select(est.var)

# optimal GMM
load("realdata_optimal_GMM.Rdata")
opt.gmm.boy <- data.frame(est=c(output.boy$theta,logor.boy$logOR), sd=sqrt(c(output.boy$var, logor.boy$var)))%>% 
  mutate_all(funs(round(.,3))) %>%
  mutate(est.var=paste0(est," (",sd,")")) %>% select(est.var)

opt.gmm.girl <- data.frame(est=c(output.girl$theta,logor.girl$logOR), sd=sqrt(c(output.girl$var, logor.girl$var)))%>% 
  mutate_all(funs(round(.,3))) %>%
  mutate(est.var=paste0(est," (",sd,")")) %>% select(est.var)

# pseudo-likelihood
load("realdata_order_statistic.Rdata")
pselik.boy <- data.frame(est=c(NA,NA,NA,NA,logor.boy$logor), sd=c(NA,NA,NA,NA,logor.boy$sd))%>% 
  mutate_all(funs(round(.,3))) %>% 
  mutate(est.var=paste0(est," (",sd,")")) %>% select(est.var)

pselik.girl <- data.frame(est=c(NA,NA,NA,NA,logor.girl$logor), sd=c(NA,NA,NA,NA,logor.girl$sd))%>% 
  mutate_all(funs(round(.,3))) %>%
  mutate(est.var=paste0(est," (",sd,")")) %>% select(est.var)

## data combined
# the results from optimal and non-optimal GEE are too similar, thus remove the optimal GEE results
est.dt <- rbind(cbind(gmm.girl,gmm.boy),cbind(pselik.girl,pselik.boy))
est.dt[est.dt==est.dt[6,1]] <- ""
est.dt <- matrix(unlist(lapply(est.dt, function(x) {gsub("NA", "-", x)})),nrow = 10)
est.dt <- est.dt[c(1:5,10),]

table1 <- as_hux(est.dt) %>% 
  insert_row("Non-optimal GEE","",after = 0) %>% 
  insert_row("Girls","Boys",after = 0) %>%
  insert_row("Pseudo-likelihood","",after = 7) %>% 
  insert_column("","","\\(\\theta_{11}\\)","\\(\\theta_{12}\\)","\\(\\theta_{21}\\)","\\(\\theta_{22}\\)", "log(OR)","","log(OR)",after = 0) %>%
  merge_cells(2, 2:3)%>%
  merge_cells(8, 2:3)%>%
  set_align(row=0:nrow(.), everywhere, "center") %>%
  set_align(col=0:ncol(.),everywhere,"center") %>%
  # set_align(col=2:ncol(.),everywhere,"center") %>%
  set_tb_padding(1, everywhere, 0) %>% 
  set_bold(1, everywhere)  %>%
  set_italic(row=8, everywhere)  %>%
  set_italic(row=2,everywhere) %>%
  # set_number_format(everywhere,col = 1,'%2.0f') %>%
  set_bottom_border(row = 1, col=everywhere ) %>% 
  # set_bottom_border(row = 2, col =everywhere) %>% 
  set_top_border(row=1,col=everywhere,brdr(1, "double")) %>% set_bottom_border(row = nrow(.),col = everywhere,brdr(1, "solid")) %>%
  set_font_size(9) %>% set_escape_contents(col=1, everywhere, FALSE) %>% set_caption("Estimates of girls’ and boys' obesity rates") %>%set_all_padding(1)

table1
quick_latex(table1)


## Income data====
library("xlsx")
dt <- read.xlsx("income.xlsx", 1) %>% 
  mutate(Y = log(ifelse(M4==1,income,NA)), X = log(aux), Rx=ifelse(is.na(aux),0,1),Ry=M4) %>% select(X,Y,Rx,Ry)

# parameter estimates from complete data
rho = with(na.omit(dt),cor(X,Y))
sig.x = sd(dt$X,na.rm = T)

##############################
#####Non-optimal GMM====
##############################
set.seed(7)
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
  
  E.ry = rep(NA,nrow(x))
  E.ry[Rx==1] = predict(cond.ry,type="response")
  
  # E(x|y)
  E.x1y <- alpha+beta*Y
  
  # estimating function
  c1 = (Rx*Ry/E.ry)*(X-E.x1y)
  c2 = (Rx*Ry/E.ry)*Y*(X-E.x1y)
  
  f = cbind(c1,c2)
  
  # f=0 for Rx=0 or Ry=0 which is where the NA's located
  f[is.na(f)] <- 0
  return(f)
}

initial.est <- lm(X~Y,data = dt)
ans <- gmm(f.est1,x=dt,t0=c(initial.est$coefficients[1],initial.est$coefficients[2]))

# bootstrap
f.boot <- function(d,i){
  d2 <- d[i,]
  boot.ans <- gmm(f.est1,x=d2,t0=c(initial.est$coefficients[1],initial.est$coefficients[2]))
  logOR <-boot.ans$coefficients[2]/((1-rho^2)*sig.x^2)
}

boot.GMM <- boot(dt, f.boot, R=1000)

gee <- data.frame(theta=c(ans$coefficients,ans$coefficients[2]/((1-rho^2)*sig.x^2)), sd=c(sqrt(diag(ans$vcov)),sd(boot.GMM$t)))


# save data
save(list = c("gee"),file = "income_GMM.Rdata")

##############################
#####Non-optimal GMM====
##############################

alpha0 = gee$theta[1]
beta0 = gee$theta[2]

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
  
  E.ry = rep(NA,nrow(x))
  E.ry[Rx==1] = predict(cond.ry,type="response")
  
  # E(x|y)
  E.x1y <- alpha+beta*Y
  
  # optimal f(Y) from Theorem 3
  E <- alpha0+beta0*Y
  V <- (1-rho^2)*sig.x^2
  pi0 <- coef(cond.ry)[names(coef(cond.ry))=="(Intercept)"]
  pi1 <- coef(cond.ry)[names(coef(cond.ry))=="X"]
  fy1 <- 1/(V+exp(-0.5/V*(E^2+2*V*pi0-(E-V*pi1)^2))*(V+V^2*pi1^2))
  fy2 <- Y/(V+exp(-0.5/V*(E^2+2*V*pi0-(E-V*pi1)^2))*(V+V^2*pi1^2))
  
  # returned estimating function
  c1 = (Rx*Ry/E.ry)*fy1*(X-E.x1y)
  c2 = (Rx*Ry/E.ry)*fy2*(X-E.x1y)
  f = cbind(c1,c2)
  
  # f=0 for Rx=0 or Ry=0 which is where the NA's located
  f[is.na(f)] <- 0
  return(f)
}

ans <- gmm(f.est1,x=dt,t0=c(initial.est$coefficients[1],initial.est$coefficients[2]))

# bootstrap
opt.boot.GMM <- boot(dt, f.boot, R=1000)

opt.gee <- data.frame(theta=c(ans$coefficients,ans$coefficients[2]/((1-rho^2)*sig.x^2)), sd=c(sqrt(diag(ans$vcov)),sd(opt.boot.GMM$t)))
# View(cbind(gee,opt.gee))

# save data
save(list = c("opt.gee"),file = "income_optimal_GMM.Rdata")

##############################
#####Pseudo-likelihood====
##############################

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

logor <- data.frame(logor=summary(m1)$coefficients[1],sd=summary(m1)$coefficients[2])

# save data
save(list=c("logor"),file="income_order_statistic.Rdata")

##############################
##### Make a table====
##############################

library(huxtable)
library(dplyr)
library(tidyr)
library(tidyverse)


# GMM
load("income_GMM.Rdata")
gmm <- gee %>% mutate_all(funs(round(.,3))) %>%
  mutate(est.var=paste0(theta," (",sd,")")) %>% select(est.var)

# optimal GMM
load("income_optimal_GMM.Rdata")
opt.gmm <- opt.gee %>% 
  mutate_all(funs(round(.,3))) %>%
  mutate(est.var=paste0(theta," (",sd,")")) %>% select(est.var)


# pseudo-likelihood
load("income_order_statistic.Rdata")
pselik <- c("","",paste0(round(logor$logor,3)," (",round(logor$sd,3),")"))



## data combined
est.dt <- t(cbind(gmm,opt.gmm,pselik))


table1 <- as_hux(est.dt) %>% 
  insert_column(c("Non-optimal GEE","Optimal GEE","Pseudo-likelihood")) %>%
  insert_row("","\\(\\alpha\\)","\\(\\beta\\)", "log(OR)",after = 0) %>% 
  set_align(row=1:nrow(.), everywhere, "center") %>%
  set_align(col=1,everywhere,"left") %>%
  set_align(col=2:ncol(.),everywhere,"center") %>%
  set_tb_padding(1, everywhere, 0) %>% 
  set_bold(1, everywhere)  %>%
  set_bold(col=1, everywhere)  %>%
  set_bottom_border(row = 1, col=everywhere ) %>% 
  set_top_border(row=1,col=everywhere,brdr(1, "double")) %>% set_bottom_border(row = nrow(.),col = everywhere,brdr(1, "solid")) %>%
  set_font_size(9) %>% set_escape_contents(1, everywhere, FALSE) %>% set_caption("Parameter estimates for KLIPS data") %>%set_all_padding(1)

table1
quick_latex(table1)
