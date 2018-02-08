
## Find starting points of K and r (K is assumed to be 1500 initially, as it 
## has to be more than the maximum population at any time), and r is estimated 
## using the logistic function of population/K as a linear model of days.

beetles <- data.frame(
  days = c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154),
  beetles = c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024))

# to obtain a starting estimate of r, assuming K = 1500
K_0 = 1500
r_0_vec <- log(beetles$beetles*(K_0 -2)/(K_0 - beetles$beetles)*2)/beetles$days
r_0 = mean(r_0_vec[2:10])

pop.model <- nls(beetles$beetles ~ K*2/(2 + (K - 2)*exp(-r*beetles$days)),
            
            start = list(K = K_0, r = r_0),
            data = beetles,
            trace = TRUE )


## SSE Function

SSE <- function(K, r) {
  estimated_beetles <-  K*2/(2 + (K - 2)*exp(-r*beetles$days)) 
  beetles_diff <- estimated_beetles - beetles$beetles
  SSE <- t(beetles_diff)%*%beetles_diff 
  return(SSE)
}

K_vec = seq(500, 1500, length.out = 200)
r_vec = seq(0.1, 0.2, length.out = 200)
SSE_matrix <- matrix(nrow = length(K_vec), ncol = length(r_vec))
for (i in 1:length(K_vec)) {
  for (j in 1:length(r_vec)) {
    SSE_matrix[i,j] = SSE(K_vec[i],r_vec[j])
  }
}

contour(x = K_vec, y = r_vec, z = SSE_matrix, xlab = "Population Capacity (K)", ylab = "growth parameter (r)", 
        plot.title = title("Contour Plot of Sum of Squared Errors"))



