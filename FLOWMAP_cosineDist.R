
# cosineDist <- function(x){
#   as.dist(1 - x %*% t(x) / (sqrt(rowSums(x ^ 2) %*% t(rowSums(x ^ 2))))) 
# }

# do clusterFCS on each time series from each individual treat
# make list of all FLOWMAPclusters objects


cosine_similarity_from_matrix <- function(v, m)
{
  m <- as.matrix(m[, names(v), drop = F])
  ret <- apply(m, 1, function(x, v) {return(crossprod(x, v)/sqrt(crossprod(x) * crossprod(v)))}, v)
  return(ret)
}

cosine_similarity_matrix <- function(m)
{
  ret <- t(apply(m, 1, function(x, m) {cosine_similarity_from_matrix(x, m)}, m = m))
  return(ret)
}

# The first calculates the cosine distance of a vector from all
# the row of a matrix, while the second one, given a matrix,
# calculates a matrix of similarities between all the rows

