library(mclust)
library(MASS)
library(mclust)
library(pracma)
library(clues)
library(NMI)
library(mcclust)
library(fpc)
library(MeanShift)
library(knitr)
library(formattable)
library(magrittr)
library(kableExtra)
library(mixtools)  #for ellipse
library(LPCM)

asw_cluster <- function(dat, met = "kmp", imax = 100, ks, cov_model = "EII") {
  
  if (met == "kmpp") {
    ASW <- sapply(ks, FUN=function(k) {
      cluster.stats(d = dist(dat), clustering = kmpp(dat, k, iter.max = imax)$cluster)$avg.silwidth
    })
    nClus <- ks[which.max(ASW)]
    kclus <- kmpp(dat, nClus)$cluster
  } else if (met == "khm") {
    ASW <- sapply(ks, FUN=function(k) {
        tryCatch({
          cluster.stats(d = dist(dat), clustering = khm(dat, k, max.iter = imax)$cluster)$avg.silwidth
        }, error =  function(e) {
          NA
        })
    })
    nClus <- ks[which.max(ASW)]
    kclus <- khm(dat, nClus)$cluster
  } else if (met == "gmm-em") {
    randPairs <- randomPairs(dat)
    ASW <- sapply(ks, FUN=function(k) {
      mixclust <- Mclust(dat, G = k, modelNames = cov_model,
                         initialization = list(hcPairs = randPairs))
      if(length(unique(mixclust$classification)) == 1) {
        0
      } else {
        cluster.stats(d = dist(dat), clustering = mixclust$classification)$avg.silwidth
      }
    })
    nClus <- ks[which.max(ASW)]
    mixclust <- Mclust(dat, G = nClus, modelNames = cov_model,
                       initialization = list(hcPairs = randPairs))
    kclus <- mixclust$classification
  } else if (met == "ms") {
    # options(mc.cores=2)
    kclus <- ms(dat, plotms=0, scaled = T)$cluster.label
  }
  
  # plot of ASW
  if (met != "ms") {
    plot(y=ASW, x=ks, type="l", ylab="Average Silhouette Width", xlab="Number of Clusters",
         xaxt = "n", main = paste0("Best ", met, " clustering: ASW index"), cex.main = .75)
    axis(side = 1, at = ks)
  }
  
  if (met == "ms") {
    print("Selection of the bandwidth parameter")
  }
  
  # plot of clustering
  plot(dat, col = kclus + 1, pch = 19, main = paste0("Best ", met, " clustering: ASW index"), cex.main = .75)
  
  list(clusters = kclus)
}

ch_cluster <- function(dat, met = "kmpp", imax = 100, ks, cov_model = "EII") {
  
  # CH for all clusterings and plot
  if (met == "kmpp") {
    CH <- sapply(ks, FUN=function(k) {
      cluster.stats(d = dist(dat), clustering = kmpp(dat, k, iter.max = imax)$cluster)$ch
    })
    nClus <- ks[which.max(CH)]
    kclus <- kmpp(dat, nClus)$cluster
  } else if (met == "khm") {
    CH <- sapply(ks, FUN=function(k) {
      tryCatch({
        cluster.stats(d = dist(dat), clustering = khm(dat, k, max.iter = imax)$cluster)$ch
      }, error =  function(e) {
        NA
      })
    })
    nClus <- ks[which.max(CH)]
    kclus <- khm(dat, nClus)$cluster
  } else if (met == "gmm-em") {
    randPairs <- randomPairs(dat)
    CH <- sapply(ks, FUN=function(k) {
      mixclust <- Mclust(dat, G = k, modelNames = cov_model,
                         initialization = list(hcPairs = randPairs))
      if(length(unique(mixclust$classification)) == 1) {
        0
      } else {
        cluster.stats(d = dist(dat), clustering = mixclust$classification)$ch
      }
    })
    nClus <- ks[which.max(CH)]
    mixclust <- Mclust(dat, G = nClus, modelNames = cov_model,
                       initialization = list(hcPairs = randPairs))
    kclus <- mixclust$classification
  } else if (met == "ms") {
    # kclus <- ms(dat, plotms=0, scaled = T)$cluster.label
    result <- ms.self.coverage(dat)
    kclus <- ms(dat, h = result$select$select[2], plotms=0, scaled = T)$cluster.label
  }
  
  # plot of CH
  if (met != "ms") {
    plot(y=CH, x=ks, type="l", ylab="Calinski and Harabasz index", xlab="Number of Clusters",
         xaxt = "n", main = paste0("Best ", met, " clustering: CH index"), cex.main = .75)
    axis(side = 1, at = ks)
  }
  
  if (met == "ms") {
    print("Selection of the bandwidth parameter")
  }
  
  # plot of clustering
  plot(dat, col = kclus + 1, pch = 19, main = paste0("Best ", met, " clustering: CH index"), cex.main = .75)
  
  list(clusters = kclus)
  
}

simulate <- function(dat, truth, run.seed, imax = 100, ks, cov_model = "EII") {
  
  ext.val.mat <- matrix(ncol = 2*(5 + 2), nrow = 4)
  colnames(ext.val.mat) <- rep(c("RI", "HA", "MA", "FM",
                                 "JI", "VI", "NMI"), 2)
  
  asw.clusters <- matrix(ncol = 4, nrow = nrow(dat))
  colnames(asw.clusters) <- c("kmpp", "khm", "gmm-em", "ms")
  
  par(mfrow = c(1, 2))
  set.seed(run.seed)
  asw.clusters[, 1] <- asw_cluster(dat, met = "kmpp", imax, ks, cov_model)$clusters
  asw.clusters[, 2] <- asw_cluster(dat, met = "khm", imax, ks, cov_model)$clusters
  asw.clusters[, 3] <- asw_cluster(dat, met = "gmm-em", imax, ks, cov_model)$clusters
  par(mfrow = c(1, 1))
  asw.clusters[, 4] <- asw_cluster(dat, met = "ms", imax, ks)$clusters
  
  for(i in 1:4) {
    ext.val.mat[i, 1:5] <- adjustedRand(asw.clusters[, i],
                                        truth,
                                        randMethod = c("Rand", "HA",
                                                       "MA", "FM", "Jaccard"))
    ext.val.mat[i, 6] <- vi.dist(asw.clusters[, i], truth)
    ext.val.mat[i, 7] <- NMI(cbind(1:nrow(dat), asw.clusters[, i]),
                             cbind(1:nrow(dat), truth))$value
  }
  
  ch.clusters <- matrix(ncol = 4, nrow = nrow(dat))
  colnames(ch.clusters) <- c("kmpp", "khm", "gmm-em", "ms")

  par(mfrow = c(1, 2))
  set.seed(run.seed)
  ch.clusters[, 1] <- ch_cluster(dat, met = "kmpp", imax, ks, cov_model)$clusters
  ch.clusters[, 2] <- ch_cluster(dat, met = "khm", imax, ks, cov_model)$clusters
  ch.clusters[, 3] <- ch_cluster(dat, met = "gmm-em", imax, ks, cov_model)$clusters
  par(mfrow = c(1, 1))
  ch.clusters[, 4] <- ch_cluster(dat, met = "ms", imax, ks, cov_model)$clusters

  adjustedRand(ch.clusters[, 1], truth, randMethod = c("Rand", "HA", "MA", "FM", "Jaccard"))

  for(i in 1:4) {
    ext.val.mat[i, 7+1:5] <- adjustedRand(ch.clusters[, i],
                                          truth,
                                          randMethod = c("Rand", "HA",
                                                         "MA", "FM", "Jaccard"))
    ext.val.mat[i, 7+6] <- vi.dist(ch.clusters[, i], truth)
    ext.val.mat[i, 7+7] <- NMI(cbind(1:nrow(dat), ch.clusters[, i]),
                               cbind(1:nrow(dat), truth))$value
  }
  ext.val.mat
  
}

simulate_short <- function(dat, truth, run.seed, imax = 100, ks, cov_model = "EII") {
  
  ext.val.mat <- matrix(ncol = 2*(5 + 2), nrow = 4)
  colnames(ext.val.mat) <- rep(c("RI", "HA", "MA", "FM",
                                 "JI", "VI", "NMI"), 2)
  
  asw.clusters <- matrix(ncol = 2, nrow = nrow(dat))
  colnames(asw.clusters) <- c("kmpp", "khm")
  
  par(mfrow = c(1, 2))
  set.seed(run.seed)
  asw.clusters[, 1] <- asw_cluster(dat, met = "kmpp", imax, ks, cov_model)$clusters
  asw.clusters[, 2] <- asw_cluster(dat, met = "khm", imax, ks, cov_model)$clusters
  
  for(i in 1:2) {
    ext.val.mat[i, 1:5] <- adjustedRand(asw.clusters[, i],
                                        truth,
                                        randMethod = c("Rand", "HA",
                                                       "MA", "FM", "Jaccard"))
    ext.val.mat[i, 6] <- vi.dist(asw.clusters[, i], truth)
    ext.val.mat[i, 7] <- NMI(cbind(1:nrow(dat), asw.clusters[, i]),
                             cbind(1:nrow(dat), truth))$value
  }
  
  ch.clusters <- matrix(ncol = 2, nrow = nrow(dat))
  colnames(ch.clusters) <- c("kmpp", "khm")
  
  par(mfrow = c(1, 2))
  set.seed(run.seed)
  ch.clusters[, 1] <- ch_cluster(dat, met = "kmpp", imax, ks, cov_model)$clusters
  ch.clusters[, 2] <- ch_cluster(dat, met = "khm", imax, ks, cov_model)$clusters

  
  adjustedRand(ch.clusters[, 1], truth, randMethod = c("Rand", "HA", "MA", "FM", "Jaccard"))
  
  for(i in 1:2) {
    ext.val.mat[i, 7+1:5] <- adjustedRand(ch.clusters[, i],
                                          truth,
                                          randMethod = c("Rand", "HA",
                                                         "MA", "FM", "Jaccard"))
    ext.val.mat[i, 7+6] <- vi.dist(ch.clusters[, i], truth)
    ext.val.mat[i, 7+7] <- NMI(cbind(1:nrow(dat), ch.clusters[, i]),
                               cbind(1:nrow(dat), truth))$value
  }
  ext.val.mat
  
}


run_time_simulate <- function(dat) {
  
  if (met == "kmpp") {
    ASW <- sapply(ks, FUN=function(k) {
      cluster.stats(d = dist(dat), clustering = kmpp(dat, k)$cluster)$avg.silwidth
    })
    nClus <- ks[which.max(ASW)]
    kclus <- kmpp(dat, nClus)$cluster
  } else if (met == "khm") {
    ASW <- sapply(ks, FUN=function(k) {
      cluster.stats(d = dist(dat), clustering = khm(dat, k)$cluster)$avg.silwidth
    })
    nClus <- ks[which.max(ASW)]
    kclus <- khm(dat, nClus)$cluster
  } else if (met == "gmm-em") {
    randPairs <- randomPairs(dat)
    ASW <- sapply(ks, FUN=function(k) {
      mixclust <- Mclust(dat, G = k, modelNames = "EEV",
                         initialization = list(hcPairs = randPairs))
      cluster.stats(d = dist(dat), clustering = mixclust$classification)$avg.silwidth
    })
    nClus <- ks[which.max(ASW)]
    mixclust <- Mclust(dat, G = nClus, modelNames = "EEV",
                       initialization = list(hcPairs = randPairs))
    kclus <- mixclust$classification
  } else if (met == "ms") {
    # options(mc.cores=2)
    kclus <- msClustering(t(dat), h = NULL, kernel = "gaussianKernel")$labels
  }
}




