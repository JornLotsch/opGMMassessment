# A parallel version of cluster::clusGap, modified from the R package 'cluster'
# Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2021).
# cluster: Cluster Analysis Basics and Extensions. R package version 2.1.2.
#' @importFrom stats dist runif var
#' @importFrom parallel detectCores mclapply
#' @importFrom cluster maxSE
clusGapP <- function(x, FUNcluster, K.max, B = 100, d.power = 1, spaceH0 = c("scaledPCA", "original"), nProc, ...) {
  # Check input parameters
  stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >= 2, (n <- nrow(x)) >= 1, ncol(x) >= 1)
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) {
    stop("'B' has to be a positive integer")
  }

  # Get the call
  cl. <- match.call()

  # Convert x to matrix if it's a data frame
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  # Create a sequence of row indices
  ii <- seq_len(n)

  # Define the W.k function
  W.k <- function(X, kk) {
    clus <- if (kk > 1) {
      FUNcluster(X, kk)$cluster
    } else {
      rep.int(1L, nrow(X))
    }
    0.5 * sum(vapply(split(ii, clus), function(I) {
      xs <- X[I, , drop = FALSE]
      sum(stats::dist(xs)^d.power / nrow(xs))
    }, 0))
  }

  # Initialize the vectors
  logW <- E.logW <- SE.sim <- numeric(K.max)

  # Calculate logW
  for (k in 1:K.max) {
    logW[k] <- log(W.k(x, k))
  }

  # Prepare the data for bootstrapping
  spaceH0 <- match.arg(spaceH0)
  xs <- scale(x, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  switch(spaceH0,
         scaledPCA = {
           V.sx <- svd(xs, nu = 0)$v
           xs <- xs %*% V.sx
         },
         original = {
         },
         stop("invalid 'spaceH0':", spaceH0))
  rng.x1 <- apply(xs, 2L, range)
  logWks <- matrix(0, B, K.max)
  logWks1 <- matrix(0, 1, K.max)

  # Perform the bootstrapping in parallel
  list.of.bootstraps <- 1:B
  if (nProc > 1 && !Sys.info()[["sysname"]] == "Windows") {
    BootGap <- parallel::mclapply(list.of.bootstraps, function(b) {
      z1 <- apply(rng.x1, 2, function(M, nn) stats::runif(nn, min = M[1], max = M[2]), nn = n)
      z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx), original = z1) + m.x
      for (k in 1:K.max) {
        logWks1[, k] <- log(W.k(z, k))
      }
      return(logWks1)
    }, mc.cores = nProc)
  } else {
    BootGap <- lapply(list.of.bootstraps, function(b) {
      z1 <- apply(rng.x1, 2, function(M, nn) stats::runif(nn, min = M[1], max = M[2]), nn = n)
      z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx), original = z1) + m.x
      for (k in 1:K.max) {
        logWks1[, k] <- log(W.k(z, k))
      }
      return(logWks1)
    })
  }

  # Combine the bootstrapping results
  logWks <- matrix(unlist(BootGap), ncol = K.max, nrow = B, byrow = TRUE)

  # Calculate the final results
  E.logW <- colMeans(logWks)
  SE.sim <- sqrt((1 + 1 / B) * apply(logWks, 2, stats::var))

  # Return the result
  structure(class = "clusGap", list(
    Tab = cbind(logW, E.logW, gap = E.logW - logW, SE.sim),
    call = cl., spaceH0 = spaceH0, n = n, B = B, FUNcluster = FUNcluster
  ))
}
