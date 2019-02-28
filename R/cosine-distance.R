
cosine.similarity.from.matrix <- function(v, m) {
  m <- as.matrix(m[, names(v), drop = FALSE])
  ret <- apply(m, 1, function(x, v) { return(crossprod(x, v)/sqrt(crossprod(x) * crossprod(v))) }, v)
  return(ret)
}

cosine.similarity.matrix <- function(m) {
  ret <- t(apply(m, 1, function(x, m) {cosine.similarity.from.matrix(x, m)}, m = m))
  return(ret)
}
