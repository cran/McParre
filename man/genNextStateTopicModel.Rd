\name{genNextStateTopicModel}
\alias{genNextStateTopicModel}
\title{
One-step generation function for Topic Models Gibbs Sampler
}
\description{
For the Latent Dirichlet Allocation model defined below in the reference,
this function generates one iteration of the Gibbs sampler to sample 
from the posterior distribution.
}
\usage{
genNextStateTopicModel(x, wid, did, K, numDocs, V, totalWords, alpha, betaPrior)
}
\arguments{
  \item{x}{A vector of length totalWords + K*numDocs + K*V, which is divided 
  as follows. The first totalWords elements are categorical, and they 
  represent the latent assignment of words to topics.}
  \item{wid}{A vector representing the words, as indexed by 1:V - the words in 
  the vocabulary.}
  \item{did}{A vector demarcating the documents that the words belong to.}
  \item{K}{The number of topics.}
  \item{numDocs}{The number of documents in the corpus.}
  \item{V}{The number of words in the vocabulary.}
  \item{totalWords}{The total number of words in the corpus.}
  \item{alpha}{The prior parameter on the mixture across topics. A vector of 
  length K.}
  \item{betaPrior}{The prior parameter for the Dirichlet distribution on the 
  vocabulary. A vector of length V.}
}
\details{
Topic models are used to describe the structure of a corpus, by attempting to 
classify each of the documents as a mixture of topics. The full data 
generating hierarchy can be found in the paper below.
}
\value{
A vector that represents the next state of the Gibbs sampler, 
when started at state x.
}
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
#function(x, wid, did, K, numDocs, V, totalWords,
#                                   alpha, betaPrior) {
# y <- rep(0, length(x))
#
# .C("genNextStateTopicModelInC", as.double(x), as.double(y), as.integer(wid),
#    as.integer(did), as.integer(K), as.integer(numDocs), as.integer(V), 
#    as.integer(totalWords), as.double(alpha), as.double(betaPrior),
#    PACKAGE="McParre")[[2]]
#  }
}
\keyword{ htest }
