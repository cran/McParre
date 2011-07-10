# function to plot scaled regeneration quantile plot
srq <- function(outputObj, subset, ...) {

  if(missing(subset))
    tmp <- outputObj[,"N"] else 
  {
    if(max(subset) > nrow(outputObj))
      stop("Subset index out of range.\n")
    tmp <- outputObj[subset, "N"]
  }

  R <- length(tmp)
  plot(1:R/R, cumsum(tmp)/sum(tmp), xlab=expression(i/n), 
       ylab=expression(paste(T[i]/T[n])), 
       main="Scaled Regeneration Quantile Plot", ...)

  abline(a=0, b=1, lty=2)

}
