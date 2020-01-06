rm(list=ls())
setwd("C:/Users/jrash/ubuntu_share/Google Drive/bioinformatics_phd_docs/written_prelim/reich/")

set.seed(1231)

expit<- function(x){1/(1+exp(-x))}

theta1 <- 0.1 + rexp(1,25)
theta2 <- theta1+rgamma(1,10,10)
theta3 <- rbeta(1,5,3)
theta4 <- rgamma(1,5,100)

a <- theta1 
b <- theta2
c <- -theta3/theta4
d <- 1/theta4

n <- 100
m.vec <- c(7, 10, 15)

sim.ls <- list()

for(j in 1:3) {
  m <- m.vec[j]
  x <- seq(0,1,length=m)
  Y <- matrix(0,n,m)
  
  for(i in 1:n){
    mu    <- a + b*expit(c + d*x)
    # constant error variance
    Y[i,] <- rnorm(m,mu,.05)
  }
  Y <- ifelse(Y<0,0,Y)
  
  real.p<- c(a, b, c, d, .05)
  names(real.p) <- c("a", "b", "c", "d", "res")
  
  sim.ls[[j]] <- list(x = x, y = Y, real.p = real.p, m = m, mu.true = mu)
  
}

rm(i,theta1,theta2,theta3,theta4,expit)
save.image("DRdata.RData")

# Larger data set to estimate standard errors

set.seed(1231)

expit<- function(x){1/(1+exp(-x))}

theta1 <- 0.1 + rexp(1,25)
theta2 <- theta1+rgamma(1,10,10)
theta3 <- rbeta(1,5,3)
theta4 <- rgamma(1,5,100)

a <- theta1 
b <- theta2
c <- -theta3/theta4
d <- 1/theta4

n <- 10000
m.vec <- c(7, 10, 15)

sim.ls <- list()
set.seed(1231)
for(j in 1:3) {
  m <- m.vec[j]
  x <- seq(0,1,length=m)
  Y <- matrix(0,n,m)
  
  for(i in 1:n){
    mu    <- a + b*expit(c + d*x)
    # constant error variance
    Y[i,] <- rnorm(m,mu,.05)
  }
  Y <- ifelse(Y<0,0,Y)
  
  real.p<- c(a, b, c, d, .05)
  names(real.p) <- c("a", "b", "c", "d", "res")
  
  sim.ls[[j]] <- list(x = x, y = Y, real.p = real.p, m = m, mu.true = mu)
  
}

rm(i,theta1,theta2,theta3,theta4,expit)
save.image("DRdata_big.RData")


