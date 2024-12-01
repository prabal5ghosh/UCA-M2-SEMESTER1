# Co-clustering Algorithm

# x: data
# z: row partition
# w: column partition

# Co-clustering based on continuous data

# Structure of block clustering model parameters
prop.r <- c(0.5, 0.5) # proportions for rows
prop.c <- c(0.5, 0.5) # proportions for columns
a <- vector("list", 2)
a[[1]][[1]] <- list(mu = 0, var = 1)
a[[1]][[2]] <- list(mu = 0.5, var = 1)
a[[2]][[1]] <- list(mu = -0.5, var = 1)
a[[2]][[2]] <- list(mu = -0.25, var = 1)
# Specific parameter for co-clustering (mu, sigma)
param <- list(prop.r = prop.r, prop.c = prop.c, a = a)

param$a

# Calculation of row group probabilities given column groups (default), 
# or vice versa
posterior <- function(x, w, param, row = TRUE) {
  n <- nrow(x)
  d <- ncol(x)
  K <- length(param$prop.r)
  G <- length(param$prop.c)
  # Rows' groups with fixed column groups
  if (row) {
    lnf <- matrix(0, n, K)
    for (k in 1:K) {
      for (g in 1:G) {
        data <- x[, w == g, drop = FALSE]
        lnf[, k] <- lnf[, k] + rowSums(dnorm(data, param$a[[k]][[g]]$mu,
                                             sqrt(param$a[[k]][[g]]$var), log = TRUE))
      }
      lnf[, k] <- lnf[, k] + log(param$prop.r[k])
    }
  }
  # Columns' groups with fixed row groups
  else {
    z <- w # for more coherence
    lnf <- matrix(0, d, G)
    for (g in 1:G) {
      for (k in 1:K) {
        data <- x[z == k, , drop = FALSE]
        lnf[, g] <- lnf[, g] + colSums(dnorm(data, param$a[[k]][[g]]$mu,
                                             sqrt(param$a[[k]][[g]]$var), log = TRUE))
      }
      lnf[, g] <- lnf[, g] + log(param$prop.c[g])
    }
  }
  # Calculation of posterior probabilities
  t <- prop.table(exp(sweep(lnf, 1, apply(lnf, 1, max), "-")), 1)
  return(t)
}

# SEM algorithm in block clustering:
mstep <- function(x, z, w) {
  K <- max(z)
  G <- max(w)
  param <- NULL
  param$prop.r <- prop.table(table(z))
  param$prop.c <- prop.table(table(w))
  param$a <- vector("list", K)
  for (k in 1:K) {
    param$a[[k]] <- vector("list", G)
    for (g in 1:G) {
      data <- as.vector(x[z == k, w == g])
      param$a[[k]][[g]]$mu <- mean(data)
      param$a[[k]][[g]]$var <- var(data) * (length(data) - 1) / length(data)
    }
  }
  return(param)
}

loglikelihood <- function(x, z, w, param) {
  K <- length(param$prop.r)
  G <- length(param$prop.c)
  ll <- 0
  for (k in 1:K) {
    for (g in 1:G) {
      data <- as.vector(x[z == k, w == g])
      ll <- ll + sum(dnorm(data, param$a[[k]][[g]]$var,
                           sqrt(param$a[[k]][[g]]$var), log = TRUE))
    }
  }
  return(ll)
}

# Beware of the disappearance of modality
# Version S1 - M1 - S2 - M2 (another possible version S1 - S2 - M)
# Possibility of semi-supervised learning
sem <- function(x, z, w, param, niter = 100, missz = TRUE, missw = TRUE) {
  ll <- NULL
  for (i in 1:niter) {
    K <- length(param$prop.r)
    G <- length(param$prop.c)
    if (missz) {
      tik <- posterior(x, w, param)
      if (any(is.na(tik))) {
        break
      }
      z <- apply(tik, 1, function(y) sample(K, 1, prob = y))
      z <- as.numeric(as.factor(z)) # Handling modalities disappearance
    }
    param <- mstep(x, z, w)
    if (missw) {
      sjg <- posterior(x, z, param, row = FALSE)
      if (any(is.na(sjg))) {
        break
      }
      w <- apply(sjg, 1, function(y) sample(G, 1, prob = y))
      w <- as.numeric(as.factor(w)) # Handling modalities disappearance
    }
    param <- mstep(x, z, w)
    # Completion of log-likelihood tracking:
    ll <- c(ll, loglikelihood(x, z, w, param))
  }
  return(list(z = z, w = w, param = param, ll = ll))
}

# Initialisation by random partition
init <- function(x, K, G) {
  n <- nrow(x)
  d <- ncol(x)
  z <- sample(K, n, replace = TRUE)
  w <- sample(G, d, replace = TRUE)
  param <- mstep(x, z, w)
  return(list(z = z, w = w, param = param))
}

# Example usage:
x <- matrix(rnorm(625), 25, 25)
heatmap(x)
init.val <- init(x, 3, 2)
res <- sem(x, init.val$z, init.val$w, init.val$param, 1000)


# ToDo 
# 1. Try to simulate a data matrix from the block clustering model with parameter
# param 
# 2. Try to recover the parameters with the function sem
# 3. Try to recover the parameter with the package blockmodels
install.packages("blockmodels")

