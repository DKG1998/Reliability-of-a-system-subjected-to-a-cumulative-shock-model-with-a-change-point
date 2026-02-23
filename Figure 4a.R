
#Figure 4(a): With respect to parameter \lambda


#1

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


rel <- function(t, m, H = 0.8, lambda.x = 3, l = 1, #lambda.y = 3, lambda.z = 2,
                # size1 = 1, #g1 = matrix(c(1,0), 1, 2), G1 =  matrix(c(-3,3,0,-3), 2, 2, byrow = TRUE),
                # size2 = 1, #g2 = matrix(c(1,0), 1, 2), G2 =  matrix(c(-2,2,0,-2), 2, 2, byrow = TRUE)
                #){
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


approxrel <- function(t){
  z <- 1
  m <- 2
  
  x.old <- rel(t, m)
  x.new <- rel(t, m + 1)
  z <- abs(x.new - x.old)
  
  while(z > 1e-3){
    x.old <- x.new
    m <- m + 1
    x.new <- rel(t, m + 1)
    z <- abs(x.new - x.old)
  }
  x.old}

t_values <- seq(0, 10, by = 0.2)

library(parallel)
library(doParallel)

ncores <- 4

cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterExport(cl, varlist = c("rel", "matinv", "alpha_times_matinv", "expm",
                              "integral","A", "a", "e", "I",
                              "t_values", "approxrel"))

rel_values1 <- parSapply(cl, t_values, approxrel)

stopCluster(cl)


#2

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

rel <- function(t, m, H = 0.8, lambda.x = 3, l = 2, #lambda.y = 3, lambda.z = 2,
                # size1 = 1, #g1 = matrix(c(1,0), 1, 2), G1 =  matrix(c(-3,3,0,-3), 2, 2, byrow = TRUE),
                # size2 = 1, #g2 = matrix(c(1,0), 1, 2), G2 =  matrix(c(-2,2,0,-2), 2, 2, byrow = TRUE)
                #){
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


approxrel <- function(t){
  z <- 1
  m <- 2
  
  x.old <- rel(t, m)
  x.new <- rel(t, m + 1)
  z <- abs(x.new - x.old)
  
  while(z > 1e-3){
    x.old <- x.new
    m <- m + 1
    x.new <- rel(t, m + 1)
    z <- abs(x.new - x.old)
  }
  x.old}

t_values <- seq(0, 10, by = 0.2)

library(parallel)
library(doParallel)

ncores <- 4

cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterExport(cl, varlist = c("rel", "matinv", "alpha_times_matinv", "expm",
                              "integral","A", "a", "e", "I",
                              "t_values", "approxrel"))

rel_values2 <- parSapply(cl, t_values, approxrel)

stopCluster(cl)


#3

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

rel <- function(t, m, H = 0.8, lambda.x = 3, l = 3, #lambda.y = 3, lambda.z = 2,
                # size1 = 1, #g1 = matrix(c(1,0), 1, 2), G1 =  matrix(c(-3,3,0,-3), 2, 2, byrow = TRUE),
                # size2 = 1, #g2 = matrix(c(1,0), 1, 2), G2 =  matrix(c(-2,2,0,-2), 2, 2, byrow = TRUE)
                #){
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


approxrel <- function(t){
  z <- 1
  m <- 2
  
  x.old <- rel(t, m)
  x.new <- rel(t, m + 1)
  z <- abs(x.new - x.old)
  
  while(z > 1e-3){
    x.old <- x.new
    m <- m + 1
    x.new <- rel(t, m + 1)
    z <- abs(x.new - x.old)
  }
  x.old}

t_values <- seq(0, 10, by = 0.2)

library(parallel)
library(doParallel)

ncores <- 4

cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterExport(cl, varlist = c("rel", "matinv", "alpha_times_matinv", "expm",
                              "integral","A", "a", "e", "I",
                              "t_values", "approxrel"))

rel_values3 <- parSapply(cl, t_values, approxrel)

stopCluster(cl)


plot(t_values, rel_values1, type = "l", col = "blue", lty = 1, lwd = 2.2,
     xlab = "Time (t)", ylab = "P(L > t)", main = "",
     cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.4, bty = "n",
     col.lab = "black", col.axis = "black", font.main = 2)
lines(seq(0, 10, 0.2), rel_values2, col = "red", lwd = 2.2, lty = 2)
lines(seq(0, 10, 0.2), rel_values3, col = "black", lwd = 2.2, lty = 3)
#grid(col = "lightgray", lty = "dotted")
legend("topright", legend = c(expression(lambda == 1),
                              expression(lambda == 2),
                              expression(lambda == 3)),
       col = c("blue", "red", "black"), lty = c(1, 2, 3), lwd = 2.2, bty = "n", cex = 1.5)