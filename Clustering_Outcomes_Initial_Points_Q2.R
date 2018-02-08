newton_raphson2 <- function (initial, f, fdash, maxiter, give = TRUE, tol = .Machine$double.eps) 
{
  old.guess <- initial
  for (i in seq_len(maxiter)) {
    new.guess <- old.guess - f(old.guess)/fdash(old.guess)
    jj <- f(new.guess)
    if (is.na(jj) | is.infinite(jj)) {
      break
    }
    if (near.match(new.guess, old.guess) | abs(jj) < tol) {
      if (give) {
        return(list(root = new.guess, f.root = jj, iter = i))
      }
      else {
        return(new.guess)
      }
    }
    old.guess <- new.guess
  }
  return("did not converge")
}
    
clustering <- data.frame(initial_points = seq(-pi, pi, length.out = 200))
for (i in 1:nrow(clustering)) {
  temp = newton_raphson2(clustering[i,1], score2, observed_information2, maxiter = 1000)
  if (temp == "did not converge") {
    clustering[i, 2] = "No Convergence"
  }
  else {
    clustering[i, 2] = temp$root
  }
}

unique_outputs <- unique(strtrim(clustering[,2],9))
final_cluster <- data.frame(unique_outputs)

for ( i in 1:length(unique_outputs)) {
    indices <- which(strtrim(clustering[,2], 9) == unique_outputs[i])
    inputs_in_group <- clustering[indices,1]
    final_cluster[i,2] = min(inputs_in_group)
    final_cluster[i,3] = max(inputs_in_group)
  eval(parse(text = paste0("group",i," <- cbind(inputs_in_group, rep(unique_outputs[",i,"], length(inputs_in_group)))")))
}

colnames(final_cluster) <- c("Outcome", "Min_Initial", "Max_Initial")


