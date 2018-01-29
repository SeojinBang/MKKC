# Part of MKKC package
# (c) 2018 by Seo-Jin Bang, Wei Wu, and Carnegie Mellon University
# See LICENSE for licensing.

#' @title Calculate Adjusted Rand Index
#'
#' @description Calculate Adjusted Rand Index measuring the similarity between two clustering labels.
#' @param x  A vector of class labels
#' @param y  A vector of class labels to be compared with \code{x}. The length of \code{y} should be the same as \code{x}.
#' @return This function returns Adjusted Rand Index between two label vectors \code{x} and \code{y}.
#' @references \insertRef{hubert1985comparing}{MKKC}
#' @seealso \code{\link{NormMI}}, \code{\link{Purity}}
#' @export
#' @import assertthat
#' @examples
#' # true label
#' x <- rep(1:3, each = 10)
#'
#' # predicted label
#' y <- sample(x)
#'
#' # compare x and y
#' AdjRandIndex(x, y)
#' NormMI(x, y)
#' Purity(x, y)

AdjRandIndex <- function(x, y) {

  assert_that(noNA(x))
  assert_that(noNA(y))
  x = as.vector(x)
  y = as.vector(y)
  assert_that(length(x) == length(y), msg = "x and y must have the same length!\n")

  tbl = table(x, y)
  tbl.choose = choose(tbl, 2)
  tbl.rsum.choose = choose(apply(tbl, 1, sum), 2)
  tbl.csum.choose = choose(apply(tbl, 2, sum), 2)
  n.choose = choose(sum(tbl), 2)

  adjRI = (sum(tbl.choose) - (sum(tbl.rsum.choose)*sum(tbl.csum.choose))/n.choose) / (0.5*(sum(tbl.rsum.choose)+sum(tbl.csum.choose))-sum(tbl.rsum.choose)*sum(tbl.csum.choose)/n.choose)

  return(adjRI)

}
