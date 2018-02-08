##log-likelihood function
log_likelihood <- function(x) {
  sample <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
              3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)  
  log_likelihood <- -length(sample)*log(pi)
  for (y in sample) { 
    log_likelihood <- log_likelihood - log(1 + (y - x)^2)
  }
  return(log_likelihood)
}

## score refers to the first derivative of log-likelihood function.
score <- function(x) {
  sample <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
              3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
  score <- 0
  for (y in sample) {
    score <- score - 2*(y - x)/(1 + (y - x)^2)
  }
  return(score)  
}

## observed_information is the second derivative of the log-likelihood function.
observed_information <- function(x) {
  sample <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
              3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
  observed_information <- 0
  for (y in sample) {
    observed_information <- observed_information - 2*(1 - (y - x)^2)/(1 + (y - x)^2)
  }
  return(observed_information)
}

#start refers to the initial points
start <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
table_of_newton_roots <- data.frame()
for (z in start) {
  newton <- newton_raphson(z, score, observed_information, maxiter = 1000)
  table_of_newton_roots <- rbind(table_of_newton_roots, c(z, newton$root, log_likelihood(newton$root), newton$iter))
}
colnames(table_of_newton_roots) <- c("Initial Point", "MLE Estimate", "Log-likelihood_of_estimate", "Iterations to Converge")


theta_fixed <- function(p0, alpha = 1 ,tol = 1E-6, max.iter = 100){
  sample <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
              3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
  theta_old <- p0
  theta_new <- alpha*-2*sum((theta_old - sample)/(1 + (theta_old - sample)^2))  + theta_old
  iter <- 1
  while ((abs(theta_new - theta_old) > tol) && (iter < max.iter)){
    theta_old <- theta_new
    theta_new <- -2*alpha*sum((theta_old - sample)/(1 + (theta_old - sample)^2)) + theta_old
    iter <- iter + 1
  }
  if (abs(theta_new - theta_old) > tol) {
    return(NULL)
  }
  else {
    return(c(theta_new, iter))
  }
}


start <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
table_of_fixed_points <- data.frame()
for (z in start) {
  fixed_point_1 <- theta_fixed(z, alpha = 1, tol = 1E-6, max.iter = 1000)
  fixed_point_2 <- theta_fixed(z, alpha = 0.64, tol = 1E-6, max.iter = 1000)
  fixed_point_3 <- theta_fixed(z, alpha = 0.25, tol = 1E-6, max.iter = 1000)
  if (is.null(fixed_point_1)) {fixed_point_1[1:2] = "Failed"} 
  else {fixed_point_1[1:2] = signif(fixed_point_1[1:2],7)}
  if (is.null(fixed_point_2)) {fixed_point_2[1:2] = "Failed"}
  else {fixed_point_2[1:2] = signif(fixed_point_2[1:2],7)}
  if (is.null(fixed_point_3)) {fixed_point_3[1:2] = "Failed to Converge"}
  else {fixed_point_3[1:2] = signif(fixed_point_3[1:2],7)}
  table_of_fixed_points <- rbind(table_of_fixed_points, c(z, fixed_point_1[1], fixed_point_1[2], 
                                                         fixed_point_2[1], fixed_point_2[2],
                                                         fixed_point_3[1], fixed_point_3[2]))
}
colnames(table_of_fixed_points) <- c("Starting Point", "Estimate_alpha_1", "Iterations_alpha_1", "Estimate_alpha_0.64", "Iterations_alpha_0.64",
                                     "Estimate_alpha_0.25", "Iterations_alpha_0.25")


## Constructing negative log_likelihood, and negative gradient as it should be input to the nlminb function which minimizes
start <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
sample <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
            3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
objective <- function(x){ length(sample)*log(pi) + sum(log(1 + (x - sample)^2)) }
gradient <- function(x) {2* sum((x - sample)/(1 + (x - sample)^2))}
fisher_information <- function(x){return(matrix(length(sample)/2,nrow = 1, ncol = 1))}
table_fisher <- data.frame()
for (z in start) {
fisher_conv <- nlminb(z, objective, gradient, fisher_information)
newton_fisher <- newton_raphson(fisher_conv$par, score, observed_information, maxiter = 1000)
table_fisher <- rbind(table_fisher,c(z, fisher_conv$par, fisher_conv$iterations, newton_fisher$root, newton_fisher$iter))
}
colnames(table_fisher) <- c("Initial Point", "Fisher Estimate", "Fisher_Iterations", "Newton_Fisher Estimate", "Newton_Fisher Iterations")


