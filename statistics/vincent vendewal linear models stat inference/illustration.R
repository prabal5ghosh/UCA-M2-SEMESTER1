# Logistic regression
# Newton Raphson algorithm + different test statistics 

# Compute pi(x,b)
# x : the design matrix includes a column of 1 
pi <- function(x,b){
  p = exp(x %*% b)/(1 + exp(x %*% b))
  p
}

# Basic test on a data point
x = matrix(c(1,2,2),1,3)
b = matrix(c(0.5,2,2),3,1)
pi(x,b)

# Simple test on simulated data
x = matrix(rnorm(100),100,1)
x = cbind(1,x)
x = as.matrix(x)
b = matrix(c(0.5,2),2,1)
plot(x[,2],pi(x,b), xlab = "x") # Fonctionne en multi-individus
y = rbinom(100,1,pi(x,b))
points(x[,2],y,col = "red", pch=2)
legend("bottomright",legend = c("pi(x;b)","y"), col = 1:2, pch = c(1,2))

# Number of columns : 
d = ncol(x)
# Initialisation of b with 0 :
b = matrix(0,d,1)
# L vector of log-likelihood :
L = rep(0,10)

# Newton Raphson algorithm
for (i in 1:10){ # Nb iteration fixed by advance (to improve by setting a given tolerance)
  # Computation of the gradient
  Gradient = t(x) %*% matrix(y - pi(x,b),ncol = 1)
  # Computation of Vtilde : matrix with pi(x)*(1-pi(x)) on the diagonal
  Vtilde = diag(as.vector(pi(x,b)*(1-pi(x,b))))
  # Computation of the Hessian matrix
  Hessian = -t(x) %*% Vtilde %*% x
  # Calcul computation of the variance covairance matrix
  Variance = solve(-Hessian)
  # Updating b
  b = b - solve(Hessian) %*% Gradient
  # Log-likelihood
  L[i] = sum(y*log(pi(x,b)) + (1-y)*log(1-pi(x,b)))
}
L
diff(L)
plot(L,main = "Evolution of the log-likelihood",type = "l")


b # Estimated value
Variance # Estimated variance
sqrt(diag(Variance)) # Standard value of the coefficient of each variable
Wald = b^2/diag(Variance) # Wald statitistic for each variable
Wald
# Or equivalently b/(standard deviation)
b/sqrt(diag(Variance))
# Computation of p-values for each coefficient with a chi-square distribution
pchisq(Wald,1,lower.tail = FALSE)
# Computation of AIC criterion
AIC = -2*L[10] + 2*d
AIC

# Comparison of the result with the dedicated R function
don = cbind.data.frame(y=y,x=x[,-1])
res.glm = glm(y ~ x, family = "binomial", data = don)
res.glm
temp = summary(res.glm)
temp
temp$cov.unscaled

# Conclusion: It's ok! 

### Computation of likelihood ratio test and score test

# Newton-Raphson algorithm
Newton <- function(x,y){
  d = ncol(x)
  b = matrix(0,d,1)
  L = rep(0,10) # Value of 10 arbitrarly fixed 
  for (i in 1:10){
    # Gradient
    Gradient = t(x) %*% matrix(y - pi(x,b),ncol = 1)
    # Vtilde
    Vtilde = diag(as.vector(pi(x,b)*(1-pi(x,b))))
    # Hessian matrix
    Hessian = -t(x) %*% Vtilde %*% x
    # Covariance matrix
    Variance = solve(- Hessian)
    # Calcul de l'intéré suivant
    b = b - solve(Hessian) %*% Gradient
    L[i] = sum(y*log(pi(x,b)) + (1-y)*log(1-pi(x,b)))
  }
  return(list(L = L, b = b, Variance = Variance))
}

#### Computation of LRT + score statistic: 

## LRT : Maximum likelihood under H0
# We remove the related columns in the design matrix x
mod.H0 = Newton(as.matrix(x[,-2]),y)
LH0 = mod.H0$L[10] # log-likelihood under H0
LRT = -2*(LH0 - L[10]) #
pchisq(LRT,1,lower.tail = FALSE) # Reject of null hypothesis for alpha = 0.05

## Score test
# Parameters bH0
bH0 = matrix(c(mod.H0$b,0),ncol = 1)
bH0

# Gradient in bH0
Gradient = t(x) %*% matrix(y - pi(x,bH0),ncol = 1)
Gradient

# Variance in bH0
Vtilde = diag(as.vector(pi(x,bH0)*(1-pi(x,bH0))))
Hessian = -t(x) %*% Vtilde %*% x
Variance = solve(- Hessian)
Variance
# Score statistic
score = t(Gradient) %*% Variance %*% Gradient
# p-value
pchisq(score,1,lower.tail = FALSE)

### Confidence intervals on the parameter 


alpha = 0.05
modele = Newton(x,y)  
u = qnorm(1 - alpha/2) 
b0 = modele$b[1] 
b1 = modele$b[2] 
s0 = sqrt(modele$Variance[1,1])
s1 = sqrt(modele$Variance[2,2])
ICb0 = c(b0 - u*s0, b0 + u*s0) # IC on b0
ICb1 = c(b1 - u*s1, b1 + u*s1) # IC on b1

## IC sur the odds-ratio
# Let focus on the odds-ration of  x = 3 vs x = 2 : OR(3/2)
Odds3 = pi(matrix(c(1,3),1,2),modele$b)/(1-pi(matrix(c(1,3),1,2),modele$b))
# or
exp(matrix(c(1,3),1,2) %*% modele$b)

Odds2 = pi(matrix(c(1,2),1,2),modele$b)/(1-pi(matrix(c(1,2),1,2),modele$b))
# or
exp(matrix(c(1,2),1,2) %*% modele$b)

# OR(3/2)
OR3_2 = Odds3/Odds2
OR3_2
# Ce qu'on aurait pu calculer directement par exp(b1*(60 - 50)) 
exp(modele$b[2]*(3 - 2))

# Then we deduce the confidence interval for the odds-ration
exp(ICb1*(3-2))










