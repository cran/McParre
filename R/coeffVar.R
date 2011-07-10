coeffVar <-
function(outObj) {
  tmp2 <- outObj[,"N"]

  sum((tmp2 - mean(tmp2))^2)/sum(tmp2)^2

}

