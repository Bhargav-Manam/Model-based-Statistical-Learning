# Tutorial EM for GMMs
#
# Your report must be a pdf, a jupyter notebook, or a R markdown.
#
# 1) Implementing the EM
#
# Implement (from scratch) the EM for a GMM on the variables 2 and 4 of the wine data set. Cluster the data and compare your results with k-means.
# An R file called "useful_functions.R" can be useful for EM. Apart from that, try not to use packages to implement EM.
# To assess the quality of the clustering, you may use the function classError and/or adjustedRandIndex from the Mclust package.
#
# 2) Model selection
#
# Try to find a relevant number of clusters using at least three methods among the ones seen in class: AIC, BIC, ICL, and (cross-)validated likelihood.
#
# 3) Towards higher dimensional spaces
#
# Try to model more than just two variables of the same data set. Do you find the same clusters, the same number of clusters.
#


#### First we load the data and look at it

library(pgmm)
library("mvtnorm")
library("mclust")
source("useful_functions.R")

data(wine)


X <- wine[, c(2, 4)]

y <- wine[, 1]
plot(X, col = y)
mu <- rmvnorm(3, colMeans(X), cov(X))
mu

EM_MV_GMM <- function(X, K, maxit = 50, ...) {
  n_variables <- ncol(X)
  Comp_Lik <- rep(0, maxit)
  n_obs <- nrow(X)
  T <- matrix(NA, nrow = n_obs, ncol = K)

  # Initialization of \theta
  prop <- rep(1 / K, K)
  mu <- rmvnorm(K, colMeans(X), cov(X))
  sigma <- array(cov(X), dim = c(n_variables, n_variables, K))
  # The main loop of the EM algo

  for (it in 1:maxit) {
    # The E step
    for (k in 1:K) {
      lgd <- log_gaussian_density(X, mean = mu[k,], sigma = sigma[, , k])
      T[, k] <- log(prop[k]) + lgd
    }
    T <- T / rowSums(T) %*% matrix(1, nrow = 1, ncol = K) # for normalizing

    # The M step
    for (k in 1:K) {
      nk <- sum(T[, k])
      prop[k] <- nk / n_obs
      mu[k,] <- colSums(T[, k] * X) / nk
      head(T[,k])
      head(X)
      head(T[,k] * X)

      mu_to_substract <- mu[k,] %x% rep(1, n_obs)
      diff <- as.matrix(X - mu_to_substract)
      cov <- matrix(rep(0, n_variables * n_variables), nrow = n_variables, ncol = n_variables)
      for (obs in 1:n_obs) {
        cov <- cov + T[obs, k] * (diff[obs,] %*% t(diff[obs,]))
      }
      sigma[, , k] <- cov / nk

    }

    grp <- max.col(T) # get class

    # Likelihood evaluation

    for (k in 1:K) {
      lgd <- log_gaussian_density(X, mean = mu[k,], sigma = sigma[, , k])
      log_propk_norm <- log(prop[k]) + lgd
      mult <- T[, k] * log_propk_norm
      sums <- sum(mult)
      Comp_Lik[it] <- Comp_Lik[it] + sums
    }

  }
  # Returning the results
  list(prop = prop, mu = mu, sigma = sigma, T = T, grp = grp, Comp_Lik = Comp_Lik)
}

# EM clustering plot
# em_result <- EM_MV_GMM(X, 3)
em_result <- Mclust(X, 3)
# plot(X[,1],X[,2],col=em_result$grp)
plot(X[, 1], X[, 2], col = em_result$classification)
# points(em_result$mu[, 1], em_result$mu[, 2], col='black', pch=19)
points(em_result$parameters$mean[1,], em_result$parameters$mean[2,], col = 'black', pch = 19)

# K-means clustering plot
km <- kmeans(X, 3)
plot(X[, 1], X[, 2], col = km$cluster)
points(km$centers[, 1], km$centers[, 2], col = 'black', pch = 19)

# EM clustering log-likelihood plot
#plot(em_result$Comp_Lik)

# Quality of the clustering
#classError(em_result$grp,wine[,1])
classError(em_result$classification, wine[, 1])

# 2) Model selection
#
# Try to find a relevant number of clusters using at least three methods among the ones seen in class: AIC, BIC, ICL, and (cross-)validated likelihood.
#


get_free_params <- function(X, k) {
  n_variables <- ncol(X)
  free_params <- (k - 1) +
    k * n_variables +
    (k * n_variables * (n_variables + 1) / 2)
  free_params
}

aic <- function(X, k) {
  # em_result <- EM_MV_GMM(X, k)
  em_result <- Mclust(X, k)
  #max_likelihood <- max(em_result$Comp_Lik)
  max_likelihood <- em_result$loglik

  aic <- max_likelihood - get_free_params(X, k)
  aic
}

bic <- function(X, k) {
  # em_result <- EM_MV_GMM(X, k)
  em_result <- Mclust(X, k)
  #max_likelihood <- max(em_result$Comp_Lik)
  max_likelihood <- em_result$loglik
  n_obs <- nrow(X)

  bic_addend <- -(get_free_params(X, k) / 2) * log(n_obs)
  bic <- max_likelihood + bic_addend
  bic
}

icl <- function(X, k) {
  # em_result <- EM_MV_GMM(X, k)
  em_result <- Mclust(X, k)
  # T <- em_result$T
  T <- em_result$z


  icl_penalty <- sum(T * log(T)) / 2
  icl <- bic(X, k) - icl_penalty
  icl
}

aic_values <- c()
bic_values <- c()
icl_values <- c()
clusters <- c(1:10)
for (k in clusters) {
  aic_values <- append(aic_values, aic(X, k))
  bic_values <- append(bic_values, bic(X, k))
  icl_values <- append(icl_values, icl(X, k))
}

max_value <- max(aic_values, bic_values, icl_values)
min_value <- min(aic_values, bic_values, icl_values)

plot(clusters, aic_values, type = "l", pch=19, ylim = c(min_value,max_value), col="coral", ylab = "method value")
lines(clusters, bic_values, type = "l", pch=19, col="cornflowerblue")
lines(clusters, icl_values, type = "l", pch=19, col="darkolivegreen3")
legend(x=clusters[1], y=min_value+40, legend = c('aic','bic','icl'), col=c("coral", "cornflowerblue", "darkolivegreen3"),  lty = 1, lwd = 1,)

## AIC at k=3 seems to be the best model


# 3) Towards higher dimensional spaces
#
# Try to model more than just two variables of the same data set. Do you find the same clusters, the same number of clusters.
#

X3 <- wine[,c('pH', 'Flavanoids', 'Hue')]
X4 <- wine[,c('pH', 'Ash', 'Flavanoids', 'Hue')]
X5 <- wine[,c('pH', 'Ash', 'Flavanoids', 'Hue', 'Glycerol')]

clusters <- c(2:5)
for (k in clusters) {
  em_result <- EM_MV_GMM(X3, k)
  points(em_result$parameters$mean[1,], em_result$parameters$mean[2,], col = 'black', pch = 19)
  plot(X[, 1], X[, 2], col = em_result$classification)

  # em_result <- EM_MV_GMM(X4,k)
  # points(em_result$parameters$mean[1,], em_result$parameters$mean[2,], col = 'black', pch = 19)
  # plot(X[, 1], X[, 2], col = em_result$classification)
  #
  # em_result <- EM_MV_GMM(X5,k)
  # points(em_result$parameters$mean[1,], em_result$parameters$mean[2,], col = 'black', pch = 19)
  # plot(X[, 1], X[, 2], col = em_result$classification)


}
em_result <- Mclust(X,3)

# plot(X[,1],X[,2],col=em_result$grp)
plot(X[,1],X[,2],col=em_result$classification)
# points(em_result$mu[, 1], em_result$mu[, 2], col='black', pch=19)
points(em_result$parameters$mean[1,], em_result$parameters$mean[2,], col='black', pch=19)