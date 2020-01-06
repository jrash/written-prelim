library("mclust")
library(MASS)
library(mclust)
library(pracma)

# K-means clustering function

kmpp <- function(X, k, iter.max = 100) { 
  n <- nrow(X) 
  C <- numeric(k) 
  C[1] <- sample(1:n, 1) 
  
  for (i in 2:k) { 
    dm <- distmat(X, X[C, ]) 
    pr <- apply(dm, 1, min); pr[C] <- 0 
    C[i] <- sample(1:n, 1, prob = pr) 
  } 
  # convergence criteria, repeat iterations until nothing changes?
  kmeans(X, X[C, ], algorithm = "Lloyd", iter.max = iter.max) 
} 

khm <- function(X, k, max.iter = 100) {
  
  ## ToDo: add pre-conditions
  
  centers <- as.matrix(X[sample(1:nrow(X), k), ])
  labels <- rep(0, nrow(X))
  
  K <- nrow(centers)
  N <- nrow(X)
  d <- q <- matrix(0, K, N)
  
  for (run in 1:max.iter) {
    for (i in 1:K) {
      for (j in 1:N) {
        d[i,j] <- sqrt(sum((centers[i,] - X[j,])^2))
      }
    }
    labels <- apply(d, 2, which.min) ## get the cluster labels
    d <- d + 1e-6                    ## avoid numerical issues
    
    ## calculate weights
    for (i in 1:K) {
      for (j in 1:N) {
        q[i,j] <- (d[i,j]^3)*((sum(1/d[,j]^2))^2)
      }
    }
    
    ## calculate weighted mean
    for (i in 1:K) {      
      centers[i,] <- colSums(q[i,labels == i]*X[labels == i,]) / sum(q[i,labels == i])
    }    
    
  }
  
  return(list(cluster=labels, centers=centers))  
}

euclid <- function(points1, points2) {
  distanceMatrix <- matrix(NA, nrow=dim(points1)[1], ncol=dim(points2)[1])
  for(i in 1:nrow(points2)) {
    distanceMatrix[,i] <- sqrt(rowSums(t(t(points1)-points2[i,])^2))
  }
  distanceMatrix
}

K_means <- function(x, centers, distFun, nItter) {
  clusterHistory <- vector(nItter, mode="list")
  centerHistory <- vector(nItter, mode="list")
  
  for(i in 1:nItter) {
    distsToCenters <- distFun(x, centers)
    clusters <- apply(distsToCenters, 1, which.min)
    centers <- apply(x, 2, tapply, clusters, mean)
    # Saving history
    clusterHistory[[i]] <- clusters
    centerHistory[[i]] <- centers
  }
  
  list(clusters=clusterHistory, centers=centerHistory)
}

kmpp_slow <- function(X, k, iter.max = 100) { 
  n <- nrow(X) 
  C <- numeric(k) 
  C[1] <- sample(1:n, 1) 
  
  for (i in 2:k) { 
    dm <- distmat(X, X[C, ]) 
    pr <- apply(dm, 1, min); pr[C] <- 0 
    C[i] <- sample(1:n, 1, prob = pr) 
  } 
  # convergence criteria, repeat iterations until nothing changes?
  K_means(X, X[C, ], euclid, nItter = iter.max) 
} 


