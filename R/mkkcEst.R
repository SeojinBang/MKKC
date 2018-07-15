# Part of MKKC package
# (c) 2018 by Seo-Jin Bang, Wei Wu, and Carnegie Mellon University
# See LICENSE for licensing.

#' An internal function called \code{mkkcEst}.
#'
#' It performs multiple kernel K-means clustering on a multi-view data.
#'
#' @param K \eqn{N x N x P} array containing \eqn{P} kernel matrices with size \eqn{N x N}.
#' @param centers The number of clusters, say \eqn{k}.
#' @param iter.max The maximum number of iterations allowed. The default is 10.
#' @param A Linear constraint matrix.
#' @param bc Lower and upper constraint bounds.
#' @param epsilon Convergence threshold. The default is \eqn{10^{-4}}.
#' @return \code{mkkcEst} returns the following components:
#' \describe{
#'   \item{cluster}{A vector of integers (from \code{1:k}) indicating the cluster to which each point is allocated.}
#'   \item{totss}{The total sum of squares.}
#'   \item{withinss}{Matrix of within-cluster sum of squares by cluster, one row per view.}
#'   \item{withinsscluster}{Vector of within-cluster sum of squares, one component per cluster.}
#'   \item{withinssview}{Vector of within-cluster sum of squares, one component per view.}
#'   \item{tot.withinss}{Total within-cluster sum of squares, i.e. \code{sum(withinsscluster)}.}
#'   \item{betweenssview}{Vector of  between-cluster sum of squares, one component per view.}
#'   \item{tot.betweenss}{The between-cluster sum of squares, i.e. \code{totss-tot.withinss}.}
#'   \item{clustercount}{The number of clusters, say \code{k}.}
#'   \item{coefficients}{The kernel coefficients}
#'   \item{size}{The number of points, one component per cluster.}
#'   \item{iter}{The number of iterations.}
#' }
#' @keywords internal
#' @import assertthat Rmosek
#' @importFrom Matrix Matrix
#' @importFrom stats kmeans
mkkcEst = function(K, centers, iter.max = 10, A = NULL, bc = NULL, epsilon = 1e-04) {

  assert_that(centers > 0, msg = "the number of cluster should be an integer larger than 0.")
  assert_that(round(centers) == centers, msg = "the number of cluster should be an integer larger than 1.")
  assert_that(is.array(K), msg = "K should be a N x N x P array.")
  assert_that(are_equal(length(dim(K)), 3), msg = "K should be a N x N x P array.")

  Km <- K
  state <- list()
  P <- dim(Km)[3]

  assert_that(is.numeric(P), msg = "K should have at least one view.")

  ## initialize theta
  theta0 <- theta <- rep(1 / P, P)

  ## initialize combined kernel matrix
  Ktheta <- matrix(0, nrow(Km), ncol(Km))
  for (m in 1:P) {
    assert_that(noNA(Km[,,m]), msg = "NA/NaN in argument.")
    Ktheta <- Ktheta + theta[m] * Km[,,m]
  }

  ## initialize H
  H.svd <- eigen(Ktheta, symmetric = TRUE)
  H <- H.svd$vectors[, 1:centers] # normalized to unit length already

  ## iteration start
  for (iter in 1:iter.max) {

    # problem setting
    problem <- list()
    problem$sense <- "max"
    problem$c <- c(sapply(1:P, function(m) WithinClusterSS(K = Km[,,m], H = H)), 0)
    problem$A <- Matrix(0, ncol = P+1, byrow = T, sparse = TRUE)
    problem$bc <- rbind(blc = c(0), buc =  c(0))
    if (is.null(A)==FALSE & is.null(bc) == FALSE) {
      problem$A <- A
      problem$bc <- bc
    }
    assert_that(is.null(A) == is.null(bc), msg = "Error: both A and bc should be assigned.")
    problem$bx <- rbind(blx = c(rep(0, P), 1),
                        bux = c(rep(Inf, P), 1))
    problem$cones <- cbind(list("QUAD", c(P+1, 1:P)))
    rownames(problem$cones) <- c("type", "sub")

    opts <- list()
    opts$verbose <- 0

    # optimization
    result <- mosek(problem, opts)

    # update theta
    if (iter == 1) {
      theta0 <- theta
    } else {
      theta0 <- theta
      theta <- result$sol$itr$xx[1:P]
    }

    # update combined kernel matrix
    Ktheta <- matrix(0, nrow(Km), ncol(Km))
    for (m in 1:P) {
      Ktheta <- Ktheta + theta[m] * Km[,,m]
    }

    # Update H
    H.svd <- eigen(Ktheta, symmetric = TRUE)
    H <- H.svd$vectors[, 1:centers] # normalized to unit length already

    # Break if converges
    if (norm(theta0-theta, "2") < epsilon & iter > 1) {
      break()
    }

    if (iter == iter.max) {
      message(paste0("did not converge in ", iter, " iterations"))
    }
  }

  ## Recover cluster index
  Hnorm <- H / matrix(sqrt(rowSums(H^2, 2)), nrow(H), centers, byrow = FALSE)
  res <- kmeans(Hnorm, centers = centers, iter.max = 500, nstart = 1000)
  cluster <- res$cluster
  Z <- LabelToBinaryMat(cluster) # binary cluster indicator matrix
  state$cluster <- cluster

  ## Within Cluster Sum of Square
  withinSStotal = WithinClusterSS(K = Ktheta, H = Z)
  withinSSviews = sapply(1:P, function(m) WithinClusterSS(K = Km[,,m], H = Z))
  names(withinSSviews) = paste0("view", 1:P)
  withinSScluster = sapply(1:centers, function(cl) WithinClusterSS(K = Ktheta[which(Z[,cl]>0), which(Z[,cl]>0)], H = Z[which(Z[,cl]>0),cl]))
  names(withinSScluster) = paste0("cluster", 1:centers)
  withinSS = sapply(1:P, function(m) sapply(1:centers, function(cl) WithinClusterSS(K = Km[which(Z[,cl]>0),which(Z[,cl]>0),m], H = Z[which(Z[,cl]>0),cl])))
  withinSS = t(withinSS)
  colnames(withinSS) = paste0("cluster", 1:centers)
  rownames(withinSS) = paste0("view", 1:P)

  ## Between Cluster Sum of Squar
  btwSStotal = BtwClusterSS(K = Ktheta, H = Z)
  btwSSviews = sapply(1:P, function(m) BtwClusterSS(K = Km[,,m], H = Z))
  names(btwSSviews) = paste0("view", 1:P)

  ## Summary
  state$totss <- withinSStotal + btwSStotal
  state$withinss <- withinSS
  state$withinsscluster <- withinSScluster
  state$withinssview <- withinSSviews
  state$tot.withinss <- withinSStotal
  state$betweenssview <- btwSSviews
  state$tot.betweenss <- btwSStotal
  state$clustercount <- centers
  state$coefficients <- theta
  state$size <- table(state$cluster)
  state$iter <- iter

  return(state)
}

#' @keywords internal
WithinClusterSS <- function(K, H) {
  sum(diag(K)) - sum(diag(t(H) %*% K %*% H))
}

#' @keywords internal
BtwClusterSS <- function(K, H) {
  sum(diag(t(H) %*% K %*% H))
}

#' @keywords internal
LabelToBinaryMat <- function(x){

  x = as.matrix(as.character(x))
  n = dim(x)[1]
  tbl = table(x)
  mylabel.list = names(tbl)
  K = length(tbl)

  res = matrix(0, ncol = K, nrow = n)
  for (lb in 1:K) {
    mylabel = mylabel.list[lb]
    mylabel.count = tbl[lb]
    res[which(x == mylabel),lb] = 1/sqrt(mylabel.count)
  }

  return(res)
}
