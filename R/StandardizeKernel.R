# Part of MKKC package
# (c) 2018 by Seojin Bang
# See LICENSE for licensing.

#' @title Standardize a Kernel Matrix
#'
#' @description \code{StandardizeKernel} is used to standardize kernel matrices obtained from multiple views so that the views are comparable to each other. A kernel matrix \code{x} is centered and then scaled by dividing by its trace.
#' @param x A kernel matrix to be standardized.
#' @param center Logical. If \code{TRUE}, the kernel matrix \code{x} will be centered. The default is \code{TRUE}
#' @param scale Logical. If \code{TRUE}, the kernel matrix \code{x} will be scaled. The default is \code{TRUE}.
#' @export
#' @import assertthat
#' @references \insertRef{scholkopf1998nonlinear}{MKKC}
#'
#' \insertRef{bang2018mkkc}{MKKC}
#' @return This function returns the standardized kernel matrix. If both arguments \code{center} and \code{scale} are \code{TRUE}, the kernel matrix will be centered first and then scaled.
#' @examples
#' x <- diag(3); x
#' StandardizeKernel(x, center = TRUE, scale = FALSE)   # centered
#' StandardizeKernel(x, center = FALSE, scale = TRUE)   # scaled
#' StandardizeKernel(x, center = TRUE, scale = TRUE)    # centered and scaled
StandardizeKernel <- function(x, center = TRUE, scale = TRUE) {

  assert_that(is.logical(center), msg = "center is not logical.")
  assert_that(is.logical(scale), msg = "scale is not logical.")
  assert_that(noNA(x), msg = "NA/NaN in argument.")
  assert_that(all(is.numeric(x)), msg = "argument is not numeric or logical.")
  assert_that(all(is.finite(x)), msg = "Inf/-Inf in argument.")
  assert_that(are_equal(t(x), x), msg = "argument is not symmetric.")
  assert_that(all(eigen(x, symmetric = TRUE)$values > -1e-9), msg = "argument is not positive semi-definite.")

  res = as.matrix(x)
  res.size = ncol(x)
  assert_that(res.size > 1, msg = "argument should have more than one variable (column).")

  if (center == TRUE) {
    unit1 = matrix(1/res.size, ncol = res.size, nrow = res.size)
    res = res - unit1 %*% res - res %*% unit1 + unit1 %*% res %*% unit1
  }
  if (scale == TRUE) {
    res = res/sum(diag(res))
  }

  return(res)

}
