#' Combine two p-values via weighted Cauchy
#'
#' @param p1 Numeric, first p-value (SAIGE)
#' @param p2 Numeric, second p-value (GAUSS)
#' @param w1 Weight for p1 (default 0.9)
#' @param w2 Weight for p2 (default 0.1)
#' @return Numeric, combined p-value
#' @export
cauchy_combine <- function(p1, p2, w1 = 0.9, w2 = 0.1) {
  adj <- function(x) pmin(pmax(x, 1e-16), 1 - 1e-16)
  t1  <- w1 * tan((0.5 - adj(p1)) * pi)
  t2  <- w2 * tan((0.5 - adj(p2)) * pi)
  T   <- (t1 + t2) / (w1 + w2)
  return(0.5 - atan(T) / pi)
}
