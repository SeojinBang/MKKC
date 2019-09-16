# Part of MKKC package
# (c) 2018 by Seojin Bang
# See LICENSE for licensing.

#' @title Calculate Purity
#'
#' @description Calculate Purity measuring the mutual dependence between two clustering labels.
#' @param x  A vector of the predicted class labels
#' @param y  A vector of the true class labels. The length of \code{y} should be the same as \code{x}.
#' @return This function returns Purity between two label vectors \code{x} and \code{y}.
#' @references \insertRef{Manning2008information}{MKKC}
#' @seealso \code{\link{AdjRandIndex}}, \code{\link{NormMI}}
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
Purity = function(x, y) {

  assert_that(noNA(x))
  assert_that(noNA(y))
  x = as.vector(x)
  y = as.vector(y)
  assert_that(length(x) == length(y), msg = "x and y must have the same length!\n")

  tbl = table(x, y)
  tbl.major = apply(tbl, 1, max)
  n = sum(tbl)

  purity = sum(tbl.major)/ n
  return(purity)

}
