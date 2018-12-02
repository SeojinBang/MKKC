# Part of MKKC package
# (c) 2018 by Seo-Jin Bang, Wei Wu, and Carnegie Mellon University
# See LICENSE for licensing.

#' @title Multiple Kernel K-means Clustering
#'
#' @name mkkc
#' @aliases MultipleKernelKmeans.default
#' @description Performs multiple kernel K-means clustering on a multi-view data.
#' @param K \eqn{N x N x P} array containing \eqn{P} kernel matrices with size \eqn{N x N}.
#' @param centers The number of clusters, say \eqn{k}.
#' @param iter.max The maximum number of iterations allowed. The default is 10.
#' @param A \eqn{m x P} linear constraint matrix where \eqn{P} is the number of views and \eqn{m} is the number of constrints.
#' @param bc \eqn{2 x m} numeric matrix with the two rows representing the lower and upper constraint bounds.
#' @param epsilon Convergence threshold. The default is \eqn{10^{-4}}.
#' @param theta intial values for kernel coefficients. The default is 1/P for all views.
#' @param x Object of class inheriting from \code{mkkc}.
#' @param ... Additional arguments passed to \code{print}.
#' @export
#' @import assertthat Rmosek
#' @importFrom Matrix Matrix
#' @importFrom stats kmeans
#' @details The optimization problem is described with an array of the multiple
#' kernels (\code{K}), the number of clusters (\code{centers}), and the linear
#' constraint on the kernel coefficient (\eqn{\theta}). The constraint
#' is defined as \eqn{blc \le A \theta \le buc} where the lower bound (\eqn{blc})
#' is the first row of \code{bc} and the upper bound (\eqn{buc}) is the second row of
#' \code{bc}.
#'
#' Our method put more weight on views having weak signal for cluster information so
#' that we can utilize important, complementary information collected from all the
#' views. That is, the larger unexplained variance in a view is, the larger kernel
#' coefficient \eqn{\theta} on the view will be assigned. After combining multiple views,
#' it minimizes un-explained variance in the combined view by optimizing continuous
#' cluster assignments. Discrete cluster assignments are recovered by performing
#' K-means clustering on the (normalized) continuous cluster assignments. The K-means
#' algorithm is performed with 1000 random starts and the best result minimizing the
#' objective function is reported.
#'
#' We recommend to standardize all the original features to have zero-mean and unit-
#' variance before processing with kernel functions. The kernel matrices should be
#' constructed ahead of the algorithm. We recommend to normalized the kernels using
#' \code{\link{StandardizeKernel}} which makes the multiple views are comparable to each
#' other.
#' @references \insertRef{bang2018mkkc}{MKKC}
#' @keywords \code{\link[stats]{kmeans}}
#' @return \code{mkkc} returns an object of class "\code{MultipleKernelKmeans}" which has a \code{print} and a \code{coef} method. It is a list with at least the following components:
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
#' @seealso \code{\link[kernlab]{kernelMatrix}}, \code{\link{StandardizeKernel}}
#' @examples
#' require(kernlab)
#'
#' # define kernel
#' rbf <- rbfdot(sigma = 0.5)
#'
#' # construct kernel matrices
#' n.noise <- 3
#' dat1 <- kernelMatrix(rbf, simCnoise$view1[,1:(2 + n.noise)])
#' dat2 <- kernelMatrix(rbf, simCnoise$view2)
#' dat3 <- kernelMatrix(rbf, simCnoise$view3)
#'
#' # construct multiview data
#' K = array(NA, dim = c(nrow(dat1), ncol(dat1), 3))
#' K[,,1] = StandardizeKernel(dat1, center = TRUE, scale = TRUE)
#' K[,,2] = StandardizeKernel(dat2, center = TRUE, scale = TRUE)
#' K[,,3] = StandardizeKernel(dat3, center = TRUE, scale = TRUE)
#'
#' # perform multiple kernel k-means
#' res <- mkkc(K = K, centers = 3)
#'
#' coef(res) # kernel coefficients of the three views
#' res$cluster
#'
#' # perfom multiple kernel k-means with constraint (2 * theta3 <= theta1 and theta3 <= theta 2)
#' require(Matrix)
#' myA <- Matrix(c(1, 0, -2,
#'                 0, 1, -1), ncol = 3, byrow = TRUE, sparse = TRUE)
#' mybc <- rbind(blc = c(0, 0), buc =  c(Inf, Inf))
#' res <- mkkc(K = K, centers = 3, A = myA, bc = mybc)
#'
#' coef(res) # kernel coefficients of the three views
#' res$cluster
mkkc <- function(K, centers, iter.max = 10, A = NULL, bc = NULL, epsilon = 1e-04, theta = rep(1/dim(Km)[3], dim(Km)[3])) UseMethod("MultipleKernelKmeans")

#' @export
MultipleKernelKmeans.default <- function(K, centers, iter.max = 10, A = NULL, bc = NULL, epsilon = 1e-04, theta = rep(1/dim(Km)[3], dim(Km)[3])) {

  assert_that(is.null(A) == is.null(bc), msg = "both A and bc should be assigned.")

  if (is.null(A) & is.null(bc)) {
    myA <- A
  } else {
    P <- dim(K)[3]
    assert_that(P == ncol(A), msg = "the number of columns of A should be the same as the number of views.")
    assert_that(nrow(A) == ncol(bc), msg = "the number of rows of A should be the same as the number of columns of bc.")
    assert_that(is.matrix(bc) & is.numeric(bc) & nrow(bc) == 2, msg = "bc shuold be a numeric matrix with 2 rows")
    myA <- Matrix(0, ncol = P+1, nrow = nrow(A), byrow = T, sparse = TRUE)
    myA[,1:P] <- A
  }

  state <- mkkcEst(K = K, centers = centers, iter.max = iter.max, A = myA, bc = bc, epsilon = epsilon, theta = theta)
  state$call <- match.call()

  class(state) <- "MultipleKernelKmeans"
  return(state)

}

#' @rdname mkkc
#' @export
print.MultipleKernelKmeans <- function(x, ...) {

  cat("\nMultiple kernel K-means clustering with", x$clustercount, "clusters of sizes ", paste(x$size, collapse = ", "), "\n")
  cat("\nKernel coefficients of views:\n")
  print(x$coefficients)
  cat("\nClustering vector:\n")
  print(x$cluster)
  cat("\nWithin cluster sum of squares by cluster:\n")
  print(x$withinsscluster)
  cat("(between_SS / total_SS =  ", round(x$tot.betweenss / x$totss, 3) * 100, " %)\n")
  cat("\nWithin cluster sum of squares by cluster for each view:\n")
  print(x$withinss)
  cat("\nAvailable components:\n")
  print(names(x))

}
