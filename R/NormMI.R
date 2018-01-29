# Part of MKKC package
# (c) 2018 by Seo-Jin Bang, Wei Wu, and Carnegie Mellon University
# See LICENSE for licensing.

#' @title Calculate Normalized Mutual Information
#'
#' @description Calculate Normalized Mutual Information measuring the mutual dependence between two clustering labels.
#' @param x  A vector of class labels.
#' @param y  A vector of class labels to be compared with \code{x}. The length of \code{y} should be the same as \code{x}.
#' @return \code{NormMI} returns Normalized Mutual Information between two label vectors \code{x} and \code{y}.
#'
#' \code{MI} returns Mutual Information between two label vectors.
#'
#' \code{Entropy} returns Entropy of the label vector \code{x}.
#' @seealso \code{\link{AdjRandIndex}}, \code{\link{Purity}}
#' @export
#' @import assertthat
#' @references \insertRef{strehl2002cluster}{MKKC}
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
NormMI <- function(x, y) {

  assert_that(noNA(x))
  assert_that(noNA(y))
  x = as.vector(x)
  y = as.vector(y)
  assert_that(length(x) == length(y), msg = "x and y must have the same length!\n")

  entx = Entropy(x)
  enty = Entropy(y)
  if (entx + enty == 0) {
    normMI = 0
  } else {
    normMI = 2 * MI(x, y) / (entx + enty)
  }
  return(normMI)

}

#' @rdname NormMI
#' @export
MI <- function(x, y) {

  assert_that(all(!is.na(x)))
  assert_that(all(!is.na(y)))
  x = as.vector(x)
  y = as.vector(y)
  assert_that(length(x) == length(y), msg = "x and y must have the same length!\n")

  tbl = table(x, y)
  tbl.rsum = apply(tbl, 1, sum)
  tbl.csum = apply(tbl, 2, sum)
  n = sum(tbl)

  prob.tbl = tbl/sum(tbl)
  mi = sum(ifelse(tbl>0, prob.tbl * log(n * tbl / (tbl.rsum %o% tbl.csum), 2), 0))
  return(mi)

}

#' @rdname NormMI
#' @export
Entropy <- function(x) {

  assert_that(all(!is.na(x)))
  x = as.vector(x)
  tbl = table(x)
  prob.tbl = tbl/sum(tbl)

  etrp = - sum(ifelse(prob.tbl > 0, prob.tbl * log(prob.tbl), 0))
  return(etrp)

}
