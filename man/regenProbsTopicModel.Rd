\name{regenProbsTopicModel}
\alias{regenProbsTopicModel}
\title{
Compute Regeneration Probability Given Current and Future States in Topic Model
}
\description{
For the LDA model (Topic Model) minorized with a distinguished point 
and set as described in the reference below, this function computes the value 
of the small measure, which is needed when computing regeneration
probabilities.
}
\usage{
regenProbsTopicModel(x, y, xTil, distSetThetaAll, distSetBetaAll, K, V, wid, did, totalWords, numDocs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of length totalWords + K*numDocs + K*V, which is divided 
  as follows. The first totalWords elements are categorical, and they 
  represent the latent assignment of words to topics.}
  \item{y}{The next state of the chain.}
  \item{xTil}{ The distinguished point in the minorization.}
  \item{distSetThetaAll}{The distinguished set for theta. This should be of 
  the form Ax <= d, with the input argument containing the columns of A and d. 
  This input argument should have K columns and rows.}
  \item{distSetBetaAll}{The distinguished set for theta. This should be of 
  the form Ax <= d, with the input argument containing the columns of A and d. 
  This input argument should have V columns and rows.}
  \item{K}{The number of topics.}
  \item{V}{The number of words in the vocabulary.}
  \item{wid}{A vector representing the words, as indexed by 1:V - the words in 
  the vocabulary.}
  \item{did}{A vector demarcating the documents that the words belong to.}
  \item{totalWords}{The total number of words in the corpus.}
  \item{numDocs}{The number of documents in the corpus.}
}
\details{
See the demo file TopicModelsDemo.R for more information on constructing the 
distinguished sets and points.
}
\value{
A real value, denoting the regeneration probability for this minorization when 
moving from state x to y.  }
\references{
  Blei, D.M. and Ng, A.Y. and Jordan, M.I. (2003) Latent dirichlet allocation. 
  _Journal of Machine Learning Research_, *3*, 993-1022.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function(x, y, xTil, distSetThetaAll, distSetBetaAll, K, V, 
#         wid, did, totalWords, numDocs)
#{
## check if points are in distinguished set. Otherwise return 0
#  for(ii in 1:numDocs) {
#    theta1 <- y[(totalWords + 1 + (ii-1)*K):(totalWords + ii*K)]
#    #tmpVert <- verticesFromDistSet(distSetThetaAll[(1+(ii-1)*K):(ii*K),])
#    tmp <- checkIfInDistSet(theta1, distSetThetaAll[(1+(ii-1)*K):(ii*K),])
#    if(!tmp) return(0)
#  }
#  for(ii in 1:K) {
#    beta1 <- y[(totalWords + 1 + numDocs*K + (ii-1)*V):
#               (totalWords + numDocs*K + ii*V)]
#    tmp <- checkIfInDistSet(beta1, distSetBetaAll[(1+(ii-1)*V):(ii*V),])
#    if(!tmp) return(0)
#  }
#
#  smallFn <- 0
#  zidX <- x[1:totalWords]
#  zidTilde <- xTil[1:totalWords]
#  njZidX <- as.vector(unlist(tapply(zidX, did, tabulate, nbins=K)))
#  njZidTilde <- as.vector(unlist(tapply(zidTilde, did, tabulate, nbins=K)))
#  njDiff <- ifelse(njZidX - njZidTilde > 0, njZidX - njZidTilde, 0)
#  for(j in 1:numDocs) {
#    coeff1 <- njDiff[((j-1)*K + 1):(j*K)]
#    smallFn <- smallFn + smallFnComponent(coeff1, 
#                verticesFromDistSet(distSetThetaAll[(1+(j-1)*K):(j*K),]))
#  }
# 
#  MjiZidX <- matrix(unlist(tapply(zidX, wid, tabulate, nbins=K)), nrow=K)
#  MjiZidTilde <- matrix(unlist(tapply(zidTilde, wid, tabulate, nbins=K)), 
#                        nrow=K)
#  MjiDiff <- ifelse(MjiZidX - MjiZidTilde > 0, MjiZidX - MjiZidTilde, 0)
#  for(j in 1:K) {
#    coeff1 <- MjiDiff[j,]
#    smallFn <- smallFn + smallFnComponent(coeff1, 
#                verticesFromDistSet(distSetBetaAll[(1+(j-1)*V):(j*V),]))
#  }
#  betaY <- matrix(y[(1+totalWords+K*numDocs):length(y)], 
#                  nrow=K, ncol=V, byrow=TRUE)
#  A1 <- sum((MjiZidTilde  - MjiZidX) * log(betaY)) # log of the p2a ratio
#
#  thetaY <- matrix(y[(totalWords+1):(totalWords+K*numDocs)], nrow=K)
#  A2 <- sum((njZidTilde  - njZidX) * log(thetaY)) # log of the p2a ratio
#
#  exp(smallFn + A1 + A2)
#  }
}
\keyword{ htest }
