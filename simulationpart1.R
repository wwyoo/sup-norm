#Author: William Weimin Yoo
###############################################################################
#load libraries
library(splines)
library(MASS)
par(mai = c(0.5, 0.5, 0.1, 0.1))

##################################################################################
#true function
f0 = function(x){
  sum(sqrt(2) * I ^ (- 3 / 2) * sin(I) * cos((I - 0.5) * pi * x))
}

#----------------------------------------------------------------------------------------
#function to calculate derivatives of b-spline and form its design matrix, taken from bs()
bsprime <- function(x, derivs, df = NULL, knots = NULL, degree = 3, intercept = TRUE, Boundary.knots = range(x))
{
    nx <- names(x)
    x <- as.vector(x)
    Boundary.knots <- sort(Boundary.knots)
    ord <- 1L + (degree <- as.integer(degree))
    if(!is.null(df) && is.null(knots)){
      nIknots <- df - ord + (1L - intercept)
      knots <- if(nIknots > 0L) {
           knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
             2L)[-c(1L, nIknots + 2L)]
           stats::quantile(x, knots)
      }
    }
    Aknots <- sort(c(rep(Boundary.knots, ord), knots))

    basis <- splineDesign(Aknots, x, ord, derivs)
    n.col <- ncol(basis)
    dimnames(basis) <- list(nx, 1L:n.col)
    basis
}

#------------------------------------------------------------------------------------------
####################################################################################
#Bayes approach
####################################################################################
#calculate sigma2hat for empirical Bayes
sigma2emp <- function(y, cmiddle, B, eta){
  ydiff <- y - B %*% eta
  sigma2 = (crossprod(ydiff) - crossprod(forwardsolve(cmiddle, crossprod(B, ydiff)))) / n
  return(sigma2)
}

#returns posterior mean and variance
pfmeanf <- function(y, cmiddle, B, b, Omegainv, eta){
    ans = forwardsolve(cmiddle, crossprod(B, y) + Omegainv %*% eta)
    pmean = crossprod(t(b), backsolve(t(cmiddle), ans))
    return(pmean)
}

pfvarf <- function(cmiddle, b){ 
    pSigma = crossprod(forwardsolve(cmiddle, t(b)))
    return(pSigma)
}

#---------------------------------------------------------------------------------
###################################################################################
#frequentist approach
###################################################################################
sigma2hat <- function(y, B, cBB){
  sigma2 <- (crossprod(y) - crossprod(forwardsolve(cBB, crossprod(B, y)))) / (n - J)
  return(sigma2)
}

fhatmeanf <- function(y, B, cBB, b){
  ans <- forwardsolve(cBB, crossprod(B, y))
  fmean <- crossprod(t(b), backsolve(t(cBB), ans)) 
  return(fmean)
}

fhatvarf <- function(cBB, b){ 
  variance <- crossprod(forwardsolve(cBB, t(b)))
  return(variance)
}

#--------------------------------------------------------------------------------------
#determine the best number of splines
optimal <- function(y, x, a, Omegainv, eta){
  mse <- matrix(0, n, 2)
  for(i in 1:n){
    testy <- y[i]
    testx <- x[i]
    trainy <- y[-i]
    trainx <- x[-i]

    B <- bs(trainx, knots = a, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
    BB = crossprod(B)
    V = BB + Omegainv
    cM <- t(chol(V))
    cBB <- t(chol(BB))
    bnewx <- predict(B, testx)

    pfmean = pfmeanf(trainy, cmiddle = cM, B = B, Omegainv = Omegainv, eta = eta, b = bnewx)
    fhat = fhatmeanf(trainy, B = B, cBB = cBB, b = bnewx)

    mse[i, 1] <- (pfmean - testy) ^ 2  #1st col Bayes

    mse[i, 2] <- (fhat - testy) ^ 2  #2nd col Frequentist
  }

return(colMeans(mse))
}

#-------------------------------------------------------------------------------------------
#calculate Monte Carlo l2 norm
l2norm <- function(x){
  sqrt(mean(abs(x) ^ 2))
}

#----------------------------------------------------------------------------------------------
#functions to compute quantiles for confidence bands
gammaf <- function(x){
  bprimex <- bsprime(x, derivs = rep(1, length(x)), knots = a, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
  bx <- bs(x, knots = a, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
  A = crossprod(bprimex, bx) - crossprod(bx, bprimex)
  com = A %*% solve(BB, t(bx))
  num = crossprod(com, solve(BB, com))
  denom = crossprod(t(bx), solve(BB, t(bx)))
  gamma = (num ^ (1 / 2)) / (denom ^ (3 / 2))
  return(gamma)
}

cinfty <- function(x){
  gamma * exp(-x ^ 2 / 2)/(2 * pi) + 1 - pnorm(x) - w / 2
}

######################################################################################
#data generation
######################################################################################
set.seed(1000)
n = 2000  #number of data points
nsam = 1000 #number of samples from posterior
nrep = 1000 #number of data replicates for coverage calculations 
x <- seq(0, 1, length = n)
I = 1:9999

y0 = sapply(x, f0)
sigma0 <- sqrt(0.1)  #true sigma
y = y0 + sigma0 * rnorm(n)  #generate some data

q = 4  #cubic splines
N <- 16 #optimal in terms of mse
J <- N + q

tau2 = 1  #prior variance
eta <- rep(0, J)
Omegainv <- (1 / tau2) * diag(rep(1, J))

a = seq(from = 0, to = 1, length.out = N + 2)[-c(1, N + 2)]
B <- bs(x, knots = a, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
BB <- crossprod(B)
V = BB + Omegainv
cmiddle <- t(chol(V))
cBB <- t(chol(BB))

################################################################################
#credible band
################################################################################
w = 0.05

nnew <- 1000
xnew <- seq(from = 0, to = 1, length.out = nnew)
bnewx <- predict(B, xnew)

pfvar <- bnewx %*% solve(V, t(bnewx))

pfhat <- mvrnorm(n = 1000, mu = rep(0, nnew), Sigma = pfvar)
prmax = apply(abs(pfhat), 1, max)
prmaxquan = quantile(prmax, 1 - w) 

pfmean = pfmeanf(y, cmiddle = cmiddle, B = B, Omegainv = Omegainv, eta = eta, b = bnewx)
psigma2 <- as.numeric(sigma2emp(y, cmiddle = cmiddle, B = B, eta = eta))

radiusmax <- sqrt(psigma2) * prmaxquan / 2

#3rd row of Table 1
round(radiusmax, 3)

pfbandup <- pfmean + radiusmax
pfbandlower <- pfmean - radiusmax

f0newx <- sapply(xnew, f0)

#First row of Table 1
mean(pfbandlower <= f0newx & f0newx <= pfbandup)  #band actual coverage

#Figure 2 upper row in paper
plot(xnew, f0newx, type = "l", lty = 1, ylim = c(-1, 3), lwd = 5)
#lines(x, y, type = "p")
lines(xnew, pfmean, type = "l", lty = 3, lwd = 5)
lines(xnew, pfbandup, type = "l", lty = 5, lwd = 5)
lines(xnew, pfbandlower, type = "l", lty = 5, lwd = 5)

################################################################################
#confidence band
################################################################################
gamma <- mean(sapply(xnew, gammaf))  #Monte Carlo integration

calpha <- uniroot(cinfty, c(0, 10))[[1]]

fhatmean <- fhatmeanf(y, B = B, cBB = cBB, b = bnewx)
fhatvar <- fhatvarf(cBB = cBB, b = bnewx)
sigma2 <- sigma2hat(y, B = B, cBB = cBB)

fhatradius <- calpha * sqrt(sigma2 * diag(fhatvar))

fhatbandup <- fhatmean + fhatradius
fhatbandlower <- fhatmean - fhatradius

#4th and 5th row of Table 1
round(max(fhatradius), 3)
round(mean(fhatradius), 3)

#2nd row of Table 1
mean(fhatbandlower <= f0newx & f0newx <= fhatbandup)

#Figure 2 bottom row in paper 
plot(xnew, f0newx, type = "l", lty = 1, ylim = c(-1, 3), lwd = 5)
#lines(x, y, type = "p")
lines(xnew, fhatmean, type = "l", lty = 3, lwd = 5)
lines(xnew, fhatbandup, type = "l", lty = 5, lwd = 5)
lines(xnew, fhatbandlower, type = "l", lty = 5, lwd = 5)

##################################################################################
#pointwise credible and confidence intervals
##################################################################################
set.seed(1000)
wn <- 5 / n
z = qnorm(wn / 2, lower.tail = FALSE)
xnew <- seq(from = 0, to = 1, length = 100)
bnew <- predict(B, xnew)

#Bayes
pfmean <- pfmeanf(y = y, cmiddle = cmiddle, B = B, Omegainv = Omegainv, eta = eta, b = bnew)
psigma2 <- as.numeric(sigma2emp(y = y, cmiddle = cmiddle, B = B, eta = eta))
pfvar <- pfvarf(cmiddle = cmiddle, b = bnew)
pfsd <- z * sqrt(psigma2 * diag(pfvar))

coverageB <- rep(0, length(xnew))
coverageF <- rep(0, length(xnew))

for(i in 1:length(xnew)){
  bnew <- predict(B, xnew[i])
  yrep = y0 + sigma0 * matrix(rnorm(n * nrep), n, nrep)  #data replicates

  #Bayes
  pfmeanrep <- apply(yrep, 2, pfmeanf, cmiddle = cmiddle, B = B, Omegainv = Omegainv, eta = eta, b = bnew)
  psigma2rep <- apply(yrep, 2, sigma2emp, cmiddle = cmiddle, B = B, eta = eta)
  pfvarrep <- pfvarf(cmiddle = cmiddle, b = bnew)
  pfsd <- z * sqrt(psigma2rep * pfvarrep)

  pu = pfmeanrep + pfsd
  pl = pfmeanrep - pfsd

  #Frequentist
  fhatmeanrep <- apply(yrep, 2, fhatmeanf, B = B, cBB = cBB, b = bnew)
  sigma2rep <- apply(yrep , 2, sigma2hat, B = B, cBB = cBB)
  fhatvarrep <- fhatvarf(cBB = cBB, b = bnew)
  fsd <- z * sqrt(sigma2rep * fhatvarrep)

  fhatupper = fhatmeanrep + fsd
  fhatlower = fhatmeanrep - fsd
  
  #actual pointwise coverage 
  true <- f0(xnew[i])
  coverageB[i] <- mean(pl <= true & true <= pu) 
  coverageF[i] <- mean(fhatlower <= true & true <= fhatupper)
  print(i)
}

#Figure 1 in paper
plot(xnew, coverageB, type = "l", lwd = 2)
plot(xnew, coverageF, type = "l", lwd = 2)

####################################################################
#determine optimal number of B-splines 
####################################################################
Nmax = 30
mseBF = matrix(0, Nmax, 2)
tau2 <- 1
tau2inv <- 1 / tau2

for(i in 1:Nmax){
  set.seed(1000)
  N <- i  #optimal in terms of mse
  J <- N + q 

  eta <- rep(0, J)
  Omegainv <- tau2inv * diag(rep(1, J))

  a = seq(from = 0, to = 1, length.out = N + 2)[-c(1, N + 2)]
  mseBF[i, ] <- optimal(y = y, x = x, a = a, Omegainv = Omegainv, eta = eta)
  print(i)
}
