########################################################################################################################
#########################TABLES&FIGURES#################################################################################
########################################################################################################################

##############################################################################
#####Figure 3: OR estimation with varying sample size====
##############################################################################
rm(list=ls())
library(ggplot2)
library(dplyr)
load("/Known intercept/local_GMM_bivariate_dat.Rdata")
estimated1 <- estimated

load("/Known intercept/local_GMM_opt_est_bivariate_dat.Rdata")
estimated2 <- estimated

load("cluster_order_statistic_bivariate_dat.Rdata")
estimated3 <- estimated

# set up the parameters
mu.x <- 0.4
mu.y <- 2

sig.x <- 3
sig.y <- 1

rho <- 0.3

# missing indicator

b0.rx <- -0.5
b1.rx <- 1

b0.ry <- 2
b1.ry <- -1
b2.ry <- 0.7

# sample size
n.vec <- c(500, 1000, 2000, 4000)

# function for trans
trasOR <- function(beta){
  rho <- beta/(sig.x/sig.y)
  OR <- exp(rho/(1-rho^2)/(sig.x*sig.y))
  return(OR)
}

compare <- data.frame(n=estimated1$n, OR1=unlist(lapply(estimated1$beta, trasOR)), OR2=unlist(lapply(estimated2$beta, trasOR)), OR3=exp(estimated3$OR))
l.compare <- reshape2::melt(compare, id.vars = c("n"))

### BOXPLOT====
ggplot(l.compare, aes(x=factor(n), y=value, fill=factor(variable))) + 
  geom_boxplot()+geom_hline(yintercept=OR0, linetype="dashed", color = "red",size=0.5) +
  xlab("Sample size") + ylab("Parameter estimate")+
  theme_bw()+
  scale_y_continuous(breaks = c(1,1.1,OR0, 1.2),
                     labels = c(1,1.1,"True OR", 1.2),limits =c(1.025,1.2))+
  scale_fill_discrete(name = "Estimation methods", labels = c("Non-optimal GEE", "Optimal GEE", "Pseudo-likelihood"))+
  theme(legend.title=element_text(size=16,face="bold"),legend.text=element_text(size=14), axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, hjust = 1), axis.text.y = element_text(size = 14, hjust = 1), 
        axis.title.x = element_text(size=16, face="bold"),axis.title.y = element_text(size=16, face="bold"),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank()) 

##############################################################################
#####Table 1: Parameter estimates with varying sample size====
##############################################################################
rm(list=ls())
library(huxtable)
library(dplyr)
library(tidyr)
library(tidyverse)

load("local_GMM_bivariate_dat.Rdata")
estimated1 <- estimated

load("local_GMM_opt_est_bivariate_dat.Rdata")
estimated2 <- estimated

# set up the parameters
mu.x <- 0.4
mu.y <- 2

sig.x <- 3
sig.y <- 1

rho <- 0.3

# true value of parameter of interest
alpha0 = mu.x-rho*sig.x/sig.y*mu.y
beta0 = rho*sig.x/sig.y
OR0 = exp(rho/(1-rho^2)/(sig.x*sig.y))

# add OR estimate for m1 and m2
# function for trans beta to OR
trasOR <- function(beta){
  rho <- beta/(sig.x/sig.y)
  OR <- exp(rho/(1-rho^2)/(sig.x*sig.y))
  return(OR)
}

# model outputs
m1 <- estimated1 %>% gather(var, est, alpha:beta, factor_key=TRUE) %>% group_by(n,var) %>% 
  summarise(bias=mean(est-eval(as.name(paste0(var,"0")))), 
            SD=sd(est), 
            MSE=mean(est-eval(as.name(paste0(var,"0"))))^2) %>% 
  gather(var2,est,bias:MSE) %>% spread(var,est)

m2 <- estimated2 %>% gather(var, est, alpha:beta, factor_key=TRUE) %>% group_by(n,var) %>% 
  summarise(bias=mean(est-eval(as.name(paste0(var,"0")))), 
            SD=sd(est), 
            MSE=mean(est-eval(as.name(paste0(var,"0"))))^2) %>% 
  gather(var2,est,bias:MSE) %>% spread(var,est)

table1 <- as_hux(cbind(m1,m2[,3:4])) %>% 
  insert_row("","","Non-optimal GEE","", "Optimal GEE","", after = 0) %>% 
  merge_cells(1, 3:4) %>% 
  merge_cells(1, 5:6) %>% insert_column(rep("",nrow(.)), after = "beta...4") %>%
  set_align(col=1, everywhere, "left") %>%
  set_align(col=2:ncol(.),everywhere,"center") %>%
  set_tb_padding(1, everywhere, 0) %>% 
  set_bold(1, everywhere) %>% set_number_format(everywhere,col=3:ncol(.),4) %>% 
  set_number_format(everywhere,col = 1,'%2.0f') %>%
  set_bottom_border(row = 1, col =c(3:4,6:7) ) %>% 
  set_bottom_border(row = 2, col =everywhere) %>% 
  set_top_border(row=1,col=everywhere,brdr(1, "double")) %>% set_bottom_border(row = nrow(.),col = everywhere,brdr(1, "solid")) %>%
  set_font_size(9) %>% set_escape_contents(2, c(3:4,6:7), FALSE) %>% set_caption("Estimation with varying sample size") %>%set_all_padding(1)

table1[2,] <- c("N","Statistics","\\(\\alpha\\)","\\(\\beta\\)","","\\(\\alpha\\)","\\(\\beta\\)")
table1[,1] <- c("","N","500","","","1000","","","2000","","","4000","","")
table1
quick_latex(table1)

#####################################################################################################
#####Appendix Table 1: Standard deviation of estimators with varying correlation between X and Y====
#####################################################################################################
load("/Known intercept/local_rho_GMM_bivariate_dat.Rdata")
m1 <- estimated %>% group_by(rho) %>% summarise(SD=sd(beta))

load("/Known intercept/local_rho_GMM_opt_est_bivariate_dat.Rdata")
m2 <- estimated %>% group_by(rho) %>% summarise(SD=sd(beta))

load("local_rho_order_statistic_bivariate_dat.Rdata")
m3 <- estimated %>% group_by(rho) %>% summarise(logOR=sd(OR))

# set up the parameters
# rho
rho.vec <- c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)

# set up the parameters
mu.x <- 0.4
mu.y <- 2

sig.x <- 3
sig.y <- 1

table1 <- as_hux(cbind(m1,m2[,2],m3[,2])) %>% 
  set_align(col=1, everywhere, "left") %>%
  set_align(col=2:ncol(.),everywhere,"center") %>%
  set_tb_padding(1, everywhere, 0) %>% 
  set_bold(1, everywhere) %>% set_number_format(everywhere,col=3:ncol(.),4) %>% 
  set_bottom_border(row = 1, col =everywhere ) %>%
  set_top_border(row=1,col=everywhere,brdr(1, "double")) %>% set_bottom_border(row = nrow(.),col = everywhere,brdr(1, "solid")) %>%
  set_font_size(9) %>% set_escape_contents(1, c(1,2,3), FALSE) %>% set_caption("Standard deviation of estimators with varying correlation between $X$ and $Y$") %>%set_all_padding(1)

table1[1,] <- c("\\(\\rho\\)","\\(\\beta\\)\n(non-optimal GEE)","\\(\\beta\\)\n(optimal GEE)","logOR\n(pseudo-likelihood)")
table1
quick_latex(table1)

#######################################################################
#####Appendix Figure 1: OR estimation under model misspecification====
#######################################################################
rm(list=ls())
library(ggplot2)
library(dplyr)
load("/Known intercept/local_mis_GMM_bivariate_dat.Rdata")
estimated1 <- estimated

load("/Known intercept/local_mis_GMM_opt_est_bivariate_dat.Rdata")
estimated2 <- estimated

load("cluster_mis_order_statistic_bivariate_dat.Rdata")
estimated3 <- estimated

# set up the parameters
mu.x <- 0.4
mu.y <- 2

sig.x <- 3
sig.y <- 1

rho <- 0.3

# missing indicator

b0.rx <- -0.5
b1.rx <- 1

b0.ry <- 2
b1.ry <- -1
b2.ry <- 0.7
b3.ry <- 0.2 # the coefficient for x^2

# sample size
n.vec <- c(500, 1000, 2000, 4000)

OR0 <- exp(rho/(1-rho^2)/(sig.x*sig.y))

# function for trans beta to OR
trasOR <- function(beta){
  rho <- beta/(sig.x/sig.y)
  OR <- exp(rho/(1-rho^2)/(sig.x*sig.y))
  return(OR)
}

compare <- data.frame(n=estimated1$n, OR1=unlist(lapply(estimated1$beta, trasOR)), OR2=unlist(lapply(estimated2$beta, trasOR)), OR3=exp(estimated3$OR))
l.compare <- reshape2::melt(compare, id.vars = c("n"))

### BOXPLOT====
ggplot(l.compare, aes(x=factor(n), y=value, fill=factor(variable))) + 
  geom_boxplot()+geom_hline(yintercept=OR0, linetype="dashed", color = "red") +
  xlab("Sample size") + ylab("Parameter estimate")+
  theme() + theme_bw()+
  scale_y_continuous(breaks = c(1.05,1.1,OR0,1.15, 1.2,1.25),
                     labels = c(1.05,1.1,"True OR",1.15, 1.2,1.25))+
  scale_fill_discrete(name = "Estimation methods", labels = c("Non-optimal GEE", "Optimal GEE", "Pseudo-likelihood"))+
  theme(legend.title=element_text(size=16,face="bold"),legend.text=element_text(size=14), axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, hjust = 1), axis.text.y = element_text(size = 14, hjust = 1), 
        axis.title.x = element_text(size=16, face="bold"),axis.title.y = element_text(size=16, face="bold"),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank())

#######################################################################
#####Appendix Table 2: Estimation under model misspecification====
#######################################################################
load("local_mis_GMM_bivariate_dat.Rdata")
estimated1 <- estimated

load("local_mis_GMM_opt_est_bivariate_dat.Rdata")
estimated2 <- estimated

load("/Users/apple/Library/CloudStorage/Dropbox/Missing data_Anna/Simulation for the draft/cluster_mis_order_statistic_bivariate_dat.Rdata")
estimated3 <- estimated

# set up the parameters
mu.x <- 0.4
mu.y <- 2

sig.x <- 3
sig.y <- 1

rho <- 0.3

# parameter of interest
alpha0 = mu.x-rho*sig.x/sig.y*mu.y
beta0 = rho*sig.x/sig.y
OR0 = exp(rho/(1-rho^2)/(sig.x*sig.y))

# function for trans beta to OR
trasOR <- function(beta){
  rho <- beta/(sig.x/sig.y)
  OR <- exp(rho/(1-rho^2)/(sig.x*sig.y))
  return(OR)
}

# model outputs
m1 <- estimated1 %>% gather(var, est, alpha:beta, factor_key=TRUE) %>% group_by(n,var) %>% 
  summarise(bias=mean(est-eval(as.name(paste0(var,"0")))), 
            SD=sd(est), 
            MSE=mean(est-eval(as.name(paste0(var,"0"))))^2) %>% 
  gather(var2,est,bias:MSE) %>% spread(var,est)

m2 <- estimated2 %>% gather(var, est, alpha:beta, factor_key=TRUE) %>% group_by(n,var) %>% 
  summarise(bias=mean(est-eval(as.name(paste0(var,"0")))), 
            SD=sd(est), 
            MSE=mean(est-eval(as.name(paste0(var,"0"))))^2) %>% 
  gather(var2,est,bias:MSE) %>% spread(var,est)

table1 <- as_hux(cbind(m1,m2[,3:4])) %>% 
  insert_row("","","Non-optimal GEE","", "Optimal GEE","", after = 0) %>% 
  merge_cells(1, 3:4) %>% 
  merge_cells(1, 5:6) %>% insert_column(rep("",nrow(.)), after = "beta...4") %>%
  set_align(col=1, everywhere, "left") %>%
  set_align(col=2:ncol(.),everywhere,"center") %>%
  set_tb_padding(1, everywhere, 0) %>% 
  set_bold(1, everywhere) %>% set_number_format(everywhere,col=3:ncol(.),4) %>% 
  set_number_format(everywhere,col = 1,'%2.0f') %>%
  set_bottom_border(row = 1, col =c(3:4,6:7) ) %>% 
  set_bottom_border(row = 2, col =everywhere) %>% 
  set_top_border(row=1,col=everywhere,brdr(1, "double")) %>% set_bottom_border(row = nrow(.),col = everywhere,brdr(1, "solid")) %>%
  set_font_size(9) %>% set_escape_contents(2, c(3:4,6:7), FALSE) %>% set_caption("Estimation under model misspecification") %>%set_all_padding(1)

table1[2,] <- c("N","Statistics","\\(\\alpha\\)","\\(\\beta\\)","","\\(\\alpha\\)","\\(\\beta\\)")
table1[,1] <- c("","N","500","","","1000","","","2000","","","4000","","")
table1
quick_latex(table1)
