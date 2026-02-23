
#Figure 3(b)

library(expm)   # For matrix exponentiation
library(pracma) # For numerical integration
library(matrixcalc)

a <- function(m){matrix(c(1, rep(0, m - 1)), 1, m)}
e <- function(m){rep(1, m)} # Vector of ones
I <- function(m){diag(m)} # Identity matrix of size m

A <- function(m){
  out <- -m * diag(m)
  out[cbind(1:(m - 1), 2:m)] <- m
  out}

alpha_times_matinv <- function(m, t, l){
  matrix((m / (m + l * t))^(0:(m - 1)) / (1 + m / (l * t)), 1, m)}

matinv <- function(m, t, l){
  out <- (1 + m / (l * t))^(-1) * (m / (m + l * t))^abs(outer(1:m, 1:m, "-"))
  out[lower.tri(out)] <- 0
  out}

rel <- function(t, m, H = 0.8, lambda.x = 3, l = 2, 
                size1 = 1, g1 = matrix(1, 1, 1), G1 = matrix(-3, 1, 1),
                size2 = 1, g2 = matrix(1, 1, 1), G2 = matrix(-2, 1, 1)){
  
  r1 <- kronecker(g1, alpha_times_matinv(m, t, l))
  R1 <- kronecker(G1, I(m)) +
    kronecker(-1 * rowSums(as.matrix(G1)) %*% g1, matinv(m, t, l))
  
  r2 <- function(x){kronecker(g2, alpha_times_matinv(m, t - x, l))}
  
  R2 <- function(x){kronecker(G2, I(m)) +
      kronecker(-1 * rowSums(as.matrix(G2)) %*% g2, matinv(m, t - x, l))}
  
  r3 <- function(x){kronecker(g1, alpha_times_matinv(m, x, l))}
  
  R3 <- function(x){kronecker(G1, I(m)) +
      kronecker(-1 * rowSums(as.matrix(G1)) %*% g1, matinv(m, x, l))}
  
  P_X_greater_than_t <- function(t){exp(-lambda.x * t)}
  
  f_X <- function(x){lambda.x * exp(-lambda.x * x)}
  
  q2 <- function(x){1 - sum(r3(x))}
  
  s <- function(x){cbind(r3(x), q2(x) * r2(x))}
  
  q3 <- function(x){-1 * rowSums(R3(x)) %*% r2(x)}
  
  zermat <- function(x){array(0, dim = dim(R3(x)))}
  
  S <- function(x){rbind(cbind(R3(x), q3(x)), cbind(zermat(x), R2(x)))}
  
  term_1 <- c(r1 %*% rowSums(expm::expm(R1 * H))) * P_X_greater_than_t(t)
  
  integrand <- function(x){c(s(x) %*% rowSums(expm::expm(S(x)*H))) * f_X(x)}
  
  integrand_vec <- Vectorize(integrand)
  
  result <- integral(integrand_vec, 0, t)
  
  updated_result <- term_1 + result
  
  sytem_rel <- 1 - updated_result
  
  sytem_rel}


approxrel <- function(t,M){
  rel(t, M)}



M_values <- 1:50

library(parallel)
library(doParallel)

ncores <- 4

cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterExport(cl, varlist = c("rel", "matinv", "alpha_times_matinv", "expm",
                              "integral","A", "a", "e", "I",
                              "M_values", "approxrel"))

rel_values_1 <- parSapply(cl, M_values, approxrel, t = 1)

stopCluster(cl)





rel <- function(t, H, lambda = 2, n = 1e5, x1 = 1/3, x2 = 1/2){
  
  set.seed(1)
  
  poisson_values <- rpois(n, lambda * t)
  change_value <- rexp(n, 1/x1)
  
  samps <- sapply(1:n, function(i){
    ifelse(t <= change_value[i],
           sum(rexp(poisson_values[i], 1 / x1)),
           sum(rexp(rpois(1, lambda * change_value[i]), 1 / x1)) +
             sum(rexp(rpois(1, lambda * (t - change_value[i])), 1 / x2)))})
  
  mean(samps < H)}

truerel_1 <- rel(1, H = 0.8)


err_rel_1 <- rel_values_1 - truerel_1


library(expm)   # For matrix exponentiation
library(pracma) # For numerical integration
library(matrixcalc)

a <- function(m){matrix(c(1, rep(0, m - 1)), 1, m)}
e <- function(m){rep(1, m)} # Vector of ones
I <- function(m){diag(m)} # Identity matrix of size m

A <- function(m){
  out <- -m * diag(m)
  out[cbind(1:(m - 1), 2:m)] <- m
  out}

alpha_times_matinv <- function(m, t, l){
  matrix((m / (m + l * t))^(0:(m - 1)) / (1 + m / (l * t)), 1, m)}

matinv <- function(m, t, l){
  out <- (1 + m / (l * t))^(-1) * (m / (m + l * t))^abs(outer(1:m, 1:m, "-"))
  out[lower.tri(out)] <- 0
  out}

rel <- function(t, m, H = 0.8, lambda.x = 3, l = 2,
                size1 = 1, g1 = matrix(1, 1, 1), G1 = matrix(-3, 1, 1),
                size2 = 1, g2 = matrix(1, 1, 1), G2 = matrix(-2, 1, 1)){
  
  r1 <- kronecker(g1, alpha_times_matinv(m, t, l))
  R1 <- kronecker(G1, I(m)) +
    kronecker(-1 * rowSums(as.matrix(G1)) %*% g1, matinv(m, t, l))
  
  r2 <- function(x){kronecker(g2, alpha_times_matinv(m, t - x, l))}
  
  R2 <- function(x){kronecker(G2, I(m)) +
      kronecker(-1 * rowSums(as.matrix(G2)) %*% g2, matinv(m, t - x, l))}
  
  r3 <- function(x){kronecker(g1, alpha_times_matinv(m, x, l))}
  
  R3 <- function(x){kronecker(G1, I(m)) +
      kronecker(-1 * rowSums(as.matrix(G1)) %*% g1, matinv(m, x, l))}
  
  P_X_greater_than_t <- function(t){exp(-lambda.x * t)}
  
  f_X <- function(x){lambda.x * exp(-lambda.x * x)}
  
  q2 <- function(x){1 - sum(r3(x))}
  
  s <- function(x){cbind(r3(x), q2(x) * r2(x))}
  
  q3 <- function(x){-1 * rowSums(R3(x)) %*% r2(x)}
  
  zermat <- function(x){array(0, dim = dim(R3(x)))}
  
  S <- function(x){rbind(cbind(R3(x), q3(x)), cbind(zermat(x), R2(x)))}
  
  term_1 <- c(r1 %*% rowSums(expm::expm(R1 * H))) * P_X_greater_than_t(t)
  
  integrand <- function(x){c(s(x) %*% rowSums(expm::expm(S(x)*H))) * f_X(x)}
  
  integrand_vec <- Vectorize(integrand)
  
  result <- integral(integrand_vec, 0, t)
  
  updated_result <- term_1 + result
  
  sytem_rel <- 1 - updated_result
  
  sytem_rel}


approxrel <- function(t,M){
  rel(t, M)}


M_values <- 1:50

library(parallel)
library(doParallel)

ncores <- 4

cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterExport(cl, varlist = c("rel", "matinv", "alpha_times_matinv", "expm",
                              "integral","A", "a", "e", "I",
                              "M_values", "approxrel"))

rel_values_2 <- parSapply(cl, M_values, approxrel, t = 2)

stopCluster(cl)


rel <- function(t, H, lambda = 2, n = 1e5, x1 = 1/3, x2 = 1/2){
  
  set.seed(1)
  
  poisson_values <- rpois(n, lambda * t)
  change_value <- rexp(n, 1/x1)
  
  samps <- sapply(1:n, function(i){
    ifelse(t <= change_value[i],
           sum(rexp(poisson_values[i], 1 / x1)),
           sum(rexp(rpois(1, lambda * change_value[i]), 1 / x1)) +
             sum(rexp(rpois(1, lambda * (t - change_value[i])), 1 / x2)))})
  
  mean(samps < H)}

truerel_2 <- rel(2, H = 0.8)


err_rel_2 <- rel_values_2 - truerel_2

err_rel_2


plot(M_values, err_rel_1, col = "blue", lty = 1, lwd = 2.5,
     xlab = "Number of iterations", ylab = "Absolute error", main = "",
     ylim = c(0, 0.20), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, bty = "n",
     col.lab = "black", col.axis = "black", font.main = 2, pch = 1)


points(M_values, err_rel_2, col = "red", pch = 8)


grid(col = "lightgray", lty = "dotted")
legend("topright", legend = c(expression(t == 1), expression(t == 2)),
       col = c("blue", "red"), pch = c(1, 8),
       lty = c(NA, NA),
       lwd = 2.5, bty = "n", cex = 1.5)


