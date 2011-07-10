checkIfInDistSet <- function(vec1, distSet) {
  d <- nrow(distSet)
  Amat <- distSet[,1:(d-1)]
  dVec <- distSet[,d]

  sum(round(as.matrix(Amat) %*% vec1[-d], 9) <= round(dVec, 9)) == d
  #sum(round(Amat %*% vec1[-d], 9) <= round(dVec, 9)) == d
  # sum(Amat %*% vec1[-d] <= dVec) == d
}

verticesFromDistSet <- function(distSet) {
  d <- nrow(distSet)
  Amat <- distSet[,1:(d-1)]
  dVec <- distSet[,d]
  
  vertices <- matrix(0, nrow=d, ncol=d)

  vertices[1,1:(d-1)] <- -dVec[1:d-1]
  vertices[1,d] <- 1 - sum(vertices[1,1:(d-1)])

  for(i in 2:d) {
    currId <- 1:(d-1)
    currId <- currId[-(i-1)]
    vertices[i,currId] <- -dVec[currId]
    vertices[i, i-1] <- dVec[d] - sum(vertices[i, 1:(d-1)])
    vertices[i, d] <- 1 - sum(vertices[i, 1:(d-1)])
  }

  vertices
}

smallFnComponent <- function(coeff, vertices) {
  min(coeff %*% t(log(vertices)))
}
