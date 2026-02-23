


# Figure 2(a): 


library(mcompanion)
library(matrixcalc)
library(expm)
library(pracma)
library(doParallel)



a <- function(m)
{
  
  v = c()
  
  v[1] = 1
  
  for (i in 2:m) {
    
    v[i] = 0
    
  }
  
  return(t(v))
  
}

e <- function(m)
{
  
  v = c()
  
  v[1] = 1
  
  for (i in 2:m) {
    
    v[i] = 1
    
  }
  
  return(v)
  
}

G <- function(m)
{
  library(mcompanion)
  A  = m*(Jordan_matrix(m,m) - m*diag(m))-m*diag(m)
  return(A)
}


hpppmf <- function(t,l,n)
{
  
  x = exp((-l)*t) * ((l*t)^(n)/factorial(n))
  
  return(x)
}


pphppmf <- function(m,t,n,l)
{
  
  library(matrixcalc)
  
  x = l*t*diag(m) - G(m)
  
  xin = matrix.inverse(x)
  
  y = matrix.power(xin, n+1)
  
  z = (-1)*((l*t)^(n))*a(m)%*%y%*%G(m)%*%e(m)
  
  w = z[1]
  
  return(w)
  
}

approxhpp <- function(t,l,n)
{
  z = 1
  m = 2
  
  while(z > 0.0000008)
  {
    
    x = pphppmf(m,t,n,l)
    
    z = abs(pphppmf(m+1,t,n,l) - pphppmf(m,t,n,l))
    
    m = m+1
  }
  
  return(x)
  
}


library(ggplot2)

# Define n, t, and l (lambda)
t <- 1
l <- 1

# Define the range for n
n_values <- seq(1, 10, by = 1)

# Calculate hpppmf and approxhpp for each n
data1 <- data.frame(
  n = n_values,
  hpppmf = sapply(n_values, function(n) hpppmf(t, l, n)),
  approxhpp = sapply(n_values, function(n) approxhpp(t, l, n))
)




p = ggplot(data1, aes(x = approxhpp, y = hpppmf)) +
  geom_line(linetype = "dotted", color = "blue") +  
  geom_point(shape = 3, size = 3, color = "red") +  
  labs(
    title = "",
    x = "Approximated probability",
    y = "Exact probability"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 18),  
    axis.title.y = element_text(size = 18), 
    axis.text.x = element_text(size = 18),   
    axis.text.y = element_text(size = 18)    
  )

print(p)


ggsave(p, filename = "fig31.eps", height = 6, width = 6)