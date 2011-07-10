genNextStateTopicModel <- function(x, wid, did, K, numDocs, V, totalWords,
                                   alpha, betaPrior) {
 y <- rep(0, length(x))

 .C("genNextStateTopicModelInC", as.double(x), as.double(y), as.integer(wid),
    as.integer(did), as.integer(K), as.integer(numDocs), as.integer(V), 
    as.integer(totalWords), as.double(alpha), as.double(betaPrior),
    PACKAGE="McParre")[[2]]
}
