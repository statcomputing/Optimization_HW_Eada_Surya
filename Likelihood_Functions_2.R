# argument x refers to theta. It is named so for ease of using the curve function.
sample2 <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
             2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)

log_likelihood2 <- function(x){
  log_likelihood2 <- NA
  if (-pi <= x && x <= pi) {
    log_likelihood2 <- -length(sample2)*log(2*pi)
    for (y in sample2) {
      log_likelihood2 <- log_likelihood2 + log(1 - cos(y - x))
    }
  }
  return(log_likelihood2)
}

score2  <- function(x){
  score2 <- 0
  for (y in sample2) {
    score2 <- score2 + sin(x - y)/(1 - cos(y - x)) 
  }
  return(score2)
}

observed_information2 <- function(x) {
  observed_information2 <- 0
  for (y in sample2) {
    observed_information2 <- observed_information2 - 1/(1 - cos(y - x))
  }
  return(observed_information2)
}


