#------------------------------------------------------------------------------#
#                                                                              #
#                 Laplace-P-splines (LPS) for fast estimation                  #
#                       of the reproduction number                             #
#                Copyright,Oswaldo Gressani. All rights reserved.              #
#                                                                              #
#------------------------------------------------------------------------------#

library("blapsr")
library("EpiEstim")
library("numDeriv")
library("incidence")
library("COVID19")

K <- 10                  # Number of (cubic) B-splines in basis
dimlat <- K              # Dimension of latent vector
penorder <- 3            # Order of the penalty
alpha <- 5/100           # Level for credible intervals
toc_total <- 0           # Total time required for estimation
delay_struct <- "twoday" # Assumed delay structure

# Composition matrix for one day delay
C_oneday <- function(deltas){
  dMTu <- deltas[1]; dTuW <- deltas[2]; dWTh <- deltas[3]; dThF <- deltas[4]
  dFSa <-deltas[5]; dSaSu <- deltas[6]; dSuM <- deltas[7]
  Cmat <- matrix(0, nrow = 7, ncol = 7)
  Cmat[1,1] <- 1-dMTu; Cmat[1,2] <- dMTu
  Cmat[2,2] <- 1-dTuW; Cmat[2,3] <- dTuW
  Cmat[3,3] <- 1-dWTh; Cmat[3,4] <- dWTh
  Cmat[4,4] <- 1-dThF; Cmat[4,5] <- dThF
  Cmat[5,5] <- 1-dFSa; Cmat[5,6] <- dFSa
  Cmat[6,6] <- 1-dSaSu; Cmat[6,7] <- dSaSu
  Cmat[7,1] <- dSuM; Cmat[7,7] <- 1-dSuM
  return(Cmat)
}

# Composition matrix for a two day delay
C_twoday <- function(deltas){
  dMTu <- deltas[1]; dMW <- deltas[2];     # Monday 
  dTuW <- deltas[3]; dTuTh <- deltas[4]    # Tuesday
  dWTh <- deltas[5]; dWF <- deltas[6]      # Wednesday
  dThF <- deltas[7]; dThS <- deltas[8]     # Thursday
  dFSa <- deltas[9]; dFSu <- deltas[10]    # Friday
  dSaSu <- deltas[11]; dSaM <- deltas[12]  # Saturday
  dSuM <- deltas[13]; dSuT <- deltas[14]   # Sunday
  C <- matrix(0, nrow = 7, ncol = 7)
    C[1,1] <- 1-dMTu-dMW; C[1,2] <- dMTu; C[1,3] <- dMW
      C[2,2] <- 1-dTuW-dTuTh; C[2,3] <- dTuW; C[2,4] <- dTuTh
        C[3,3] <- 1-dWTh-dWF; C[3,4] <- dWTh; C[3,5] <- dWF
          C[4,4] <- 1-dThF-dThS; C[4,5] <- dThF; C[4,6] <- dThS
            C[5,5] <- 1-dFSa-dFSu; C[5,6] <- dFSa; C[5,7] <- dFSu
              C[6,1] <- dSaM; C[6,6] <- 1-dSaSu-dSaM; C[6,7] <- dSaSu
                C[7,1] <- dSuM; C[7,2] <- dSuT; C[7,7] <- 1-dSuM-dSuT
  return(C)
}

# Composition matrix for a weekend delay
C_weekend <- function(deltas){
  dMTu <- deltas[1]; dTuW <- deltas[2]; dWTh <- deltas[3]
  dThF <- deltas[4]; dFM <- deltas[5]; dSaM <- deltas[6]; dSuM <- deltas[7]
  C <- matrix(0, nrow = 7, ncol = 7)
    C[1,1] <- 1-dMTu; C[1,2] <- dMTu
      C[2,2] <- 1-dTuW; C[2,3] <- dTuW
        C[3,3] <- 1-dWTh; C[3,4] <- dWTh
          C[4,4] <- 1-dThF; C[4,5] <- dThF
            C[5,1] <- dFM; C[5,5] <- 1- dFM
              C[6,1] <- 1-dSaM; C[6,2] <- dSaM
                C[7,1] <- 1-dSuM; C[7,2] <- dSuM
 return(C)
}

#------------------------- Load  data -----------------------------------
#--- Last time the data below was obtained: 2022-09-02 @08h30
# dataAUS <- covid19(country="Australia", level=1,
#                    start="2020-05-15",
#                    end = "2020-09-15")
# 
# save(dataAUS,file = "dataAUSTRALIA.RData")
load("dataAUSTRALIA.RData")


Ot <- dataAUS$hosp

# Serial interval over 5 days
p <- discr_si(k = seq(1, 5), mu = 4, sigma = 2)
p <- p/sum(p)

# For the above serial interval calibration, see for instance:
# Ganyani et al. (2020). Estimating the generation interval for coronavirus
# disease (COVID-19) based on symptom onset data, March 2020. Eurosurveillance,
# 25(17), 2000257.

kmax <- length(p)                         
Ttime <- length(Ot)
tl <- 0
tr <- Ttime
if(delay_struct=="oneday" || delay_struct=="weekend"){
dim_eta <- 9
} else {dim_eta <- 16}

#-------------------------- LPS for fast estimation of Rt ----------------------

tic <- proc.time()                               # Start tic-tac of the clock!
B <- cubicbs(seq(Ttime), tl, tr, K = K)$Bmatrix  # B-spline basis 

#--- Penalty matrix (P) and precision matrix (Q) of spline parameters

D <- diag(K)
for(j in 1:penorder){D <- diff(D)}
P <- t(D) %*% D
P <- P + diag(1e-05, K)                        # Penalty matrix
Q <- function(lambda) lambda * P               # Precision matrix
a_lambda <- 1e-05; b_lambda <- 1e-05           # Hyperparameters

pmat <- matrix(rep(p, Ttime), ncol = kmax, byrow = TRUE) # Generation interval

# Function to moove the position of vector entries backwards
moove <- function(x, nmooves){
  # x is a vector
  # nmooves is an integer indicating the number of mooves
  if(nmooves == 0){
    xval <- x
  }else{
    xval <- c(tail(x,nmooves),x[1:(7-nmooves)])
  }
  return(xval)
}

#------------------------ Approximation of latent Ms 

Mtestim <- round((1/0.75)*Ot)

# plot(as.incidence(Ot))

#----- Estimation of Rt and hyperparameter vector (lambda, rho, delta)

sR <- function(theta, eta){
  
  Mmat <- matrix(0, nrow = Ttime, ncol = length(p))
  Mmat[1, 1] <- Mtestim[1]
  for(t in 2:Ttime){
    Mmat[t,] <- c(Mtestim[t], Mmat[(t - 1), 1], Mmat[(t - 1), 2],
                  Mmat[(t - 1), 3], Mmat[(t - 1), 4])
  }
  Mmat2 <- rbind(rep(1, kmax), Mmat[1:(Ttime - 1), ])
  Mpvec <- rowSums(Mmat2 * pmat)
  
  st_theta <- c()
  Rtvec <- as.numeric(exp(B %*% theta))
  mutvec <- Rtvec * Mpvec # Vector of means mut
  Tmod7 <- floor(Ttime / 7)
  
  # Define composition matrix
  if (delay_struct == "oneday") {
    Cm <- C_oneday(eta[3:dim_eta])
  } else if (delay_struct == "twoday") {
    Cm <- C_twoday(eta[3:dim_eta])
  } else if (delay_struct == "weekend") {
    Cm <- C_weekend(eta[3:dim_eta])
  }
  
  st_theta[1:7] <-  diag(apply(Cm * mutvec[1:7], 2, cumsum))
  
  if(Tmod7 == 1) {
    # if 8 <= Ttime <= 13 we have days 8,9,...,13
    Cm_augmented <- as.matrix(Cm[, 1:(Ttime %% 7)])
    mutmat <- matrix(0, nrow = 7, ncol = (Ttime %% 7))
    for (t in 1:ncol(mutmat)) {
   mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
    }
  } else if (Ttime %% 7 == 0) {
    # days 14,21,28,...
    Cm_augmented <- matrix(rep(Cm, Tmod7 - 1), ncol = (Tmod7 - 1) * 7)
    mutmat <- matrix(0, nrow = 7, ncol = (Tmod7 - 1) * 7)
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
    }
  } else if ((Tmod7 != 1) && (Ttime %% 7 != 0)) {
    # days > 14 and for which Ttime %% 7 != 0 
    Cm_augmented <- cbind(matrix(rep(Cm, Tmod7 - 1), ncol = (Tmod7 - 1) * 7),
                          Cm[, 1:(Ttime %% 7)])
    mutmat <- matrix(0, nrow = 7, ncol = ((Tmod7 - 1) * 7 + (Ttime %% 7)))
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
    }
  }
  
  st_theta[8:Ttime] <- colSums(Cm_augmented * mutmat)
  st_theta[which(st_theta <= 0)] <- 1e-8 # To avoid to input 0 in log()
  val <- eta[2] * st_theta
  return(val)
}

dsR<- function(k, theta, eta){
  
  Mmat <- matrix(0, nrow = Ttime, ncol = length(p))
  Mmat[1, 1] <- Mtestim[1]
  for(t in 2:Ttime){
    Mmat[t,] <- c(Mtestim[t], Mmat[(t - 1), 1], Mmat[(t - 1), 2],
                  Mmat[(t - 1), 3], Mmat[(t - 1), 4])
  }
  Mmat2 <- rbind(rep(1, kmax), Mmat[1:(Ttime - 1), ])
  Mpvec <- rowSums(Mmat2 * pmat)
  
  dstk <- c()
  Rtvec <- as.numeric(exp(B %*% theta))
  Tmod7 <- floor(Ttime / 7)
  mutvec <- Rtvec * Mpvec # Vector of means mut
  
  # Define composition matrix
  if (delay_struct == "oneday") {
    Cm <- C_oneday(eta[3:dim_eta])
  } else if (delay_struct == "twoday") {
    Cm <- C_twoday(eta[3:dim_eta])
  } else if (delay_struct == "weekend") {
    Cm <- C_weekend(eta[3:dim_eta])
  }
  
  dstk[1:7] <- diag(apply(Cm * mutvec[1:7] * B[, k][1:7], 2, cumsum))
  
  if(Tmod7 == 1) {
    # if 8 <= Ttime <= 13 we have days 8,9,...,13
    Cm_augmented <- as.matrix(Cm[, 1:(Ttime %% 7)])
    mutmat <- matrix(0, nrow = 7, ncol = (Ttime %% 7))
    Bkmat <- mutmat
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
      Bkmat[, t] <- moove(B[(2 + t - 1):(2 + t - 1 + 6), k], nmooves = t %% 7)
    }
  } else if (Ttime %% 7 == 0) {
    # days 14,21,28,...
    Cm_augmented <- matrix(rep(Cm, Tmod7 - 1), ncol = (Tmod7 - 1) * 7)
    mutmat <- matrix(0, nrow = 7, ncol = (Tmod7 - 1) * 7)
    Bkmat <- mutmat
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
      Bkmat[, t] <- moove(B[(2 + t - 1):(2 + t - 1 + 6), k], nmooves = t %% 7)
    }
  } else if ((Tmod7 != 1) && (Ttime %% 7 != 0)) {
    # days > 14 and for which Ttime %% 7 != 0 
    Cm_augmented <- cbind(matrix(rep(Cm, Tmod7 - 1), ncol = (Tmod7 - 1) * 7),
                          Cm[, 1:(Ttime %% 7)])
    mutmat <- matrix(0, nrow = 7, ncol = ((Tmod7 - 1) * 7 + (Ttime %% 7)))
    Bkmat <- mutmat
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
      Bkmat[, t] <- moove(B[(2 + t - 1):(2 + t - 1 + 6), k], nmooves = t %% 7)
    }
  }
  
  dstk[8:Ttime] <- colSums(Cm_augmented * mutmat * Bkmat)
  deriv_res <- eta[2] * dstk
  
  return(deriv_res)
  
}

d2sR <- function(k, l, theta, eta){
  
  Mmat <- matrix(0, nrow = Ttime, ncol = length(p))
  Mmat[1,1] <- Mtestim[1]
  for(t in 2:Ttime){
    Mmat[t,] <- c(Mtestim[t], Mmat[(t - 1), 1], Mmat[(t - 1), 2],
                  Mmat[(t - 1), 3], Mmat[(t - 1), 4])
  }
  Mmat2 <- rbind(rep(1, kmax), Mmat[1:(Ttime - 1), ])
  Mpvec <- rowSums(Mmat2 * pmat)
  
  d2stkl <- c()
  Rtvec <- as.numeric(exp(B %*% theta))
  Tmod7 <- floor(Ttime / 7)
  mutvec <- Rtvec * Mpvec # Vector of means mut
  
  # Define composition matrix
  if (delay_struct == "oneday") {
    Cm <- C_oneday(eta[3:dim_eta])
  } else if (delay_struct == "twoday") {
    Cm <- C_twoday(eta[3:dim_eta])
  } else if (delay_struct == "weekend") {
    Cm <- C_weekend(eta[3:dim_eta])
  }
  
  d2stkl[1:7] <- diag(apply(Cm * mutvec[1:7] * B[,k][1:7] * B[,l][1:7], 
                            2, cumsum))
  
  if(Tmod7 == 1) {
    # if 8 <= Ttime <= 13 we have days 8,9,...,13
    Cm_augmented <- as.matrix(Cm[, 1:(Ttime %% 7)])
    mutmat <- matrix(0, nrow = 7, ncol = (Ttime %% 7))
    Bkmat <- mutmat
    Blmat <- mutmat
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
      Bkmat[, t] <- moove(B[(2 + t - 1):(2 + t - 1 + 6), k], nmooves = t %% 7)
      Blmat[, t] <- moove(B[(2 + t - 1):(2 + t - 1 + 6), l], nmooves = t %% 7)
    }
  } else if (Ttime %% 7 == 0) {
    # days 14,21,28,...
    Cm_augmented <- matrix(rep(Cm, Tmod7 - 1), ncol = (Tmod7 - 1) * 7)
    mutmat <- matrix(0, nrow = 7, ncol = (Tmod7 - 1) * 7)
    Bkmat <- mutmat
    Blmat <- mutmat
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
      Bkmat[, t] <- moove(B[(2 + t - 1):(2 + t - 1 + 6), k], nmooves = t %% 7)
      Blmat[, t] <- moove(B[(2 + t - 1):(2 + t - 1 + 6), l], nmooves = t %% 7)
    }
  } else if ((Tmod7 != 1) && (Ttime %% 7 != 0)) {
    # days > 14 and for which Ttime %% 7 != 0 
    Cm_augmented <- cbind(matrix(rep(Cm, Tmod7 - 1), ncol = (Tmod7 - 1) * 7),
                          Cm[, 1:(Ttime %% 7)])
    mutmat <- matrix(0, nrow = 7, ncol = ((Tmod7 - 1) * 7 + (Ttime %% 7)))
    Bkmat <- mutmat
    Blmat <- mutmat
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
      Bkmat[, t] <- moove(B[(2 + t - 1):(2 + t - 1 + 6), k], nmooves = t %% 7)
      Blmat[, t] <- moove(B[(2 + t - 1):(2 + t - 1 + 6), l], nmooves = t %% 7)
    }
  }
  
  d2stkl[8:Ttime] <- colSums(Cm_augmented * mutmat * Bkmat *  Blmat)
  deriv_res <- eta[2] * d2stkl
  
  return(deriv_res)

}

#------ Log-likelihood function 

loglikR <- function(theta, eta){
  value <- sum(Ot * log(sR(theta, eta)) - sR(theta, eta))
  return(value)
}

dloglikR <- function(theta, eta){
  dsklist <- lapply(seq(K), dsR, theta=theta, eta=eta) 
  res <- sapply(Map("*",dsklist,list(((Ot/sR(theta,eta))-1))),"sum")
  return(res)
}

d2loglikR <- function(theta, eta){
  
  Hessloglik <- matrix(0, nrow = K, ncol = K)
  st <- sR(theta,eta)
  
  for(k in 1:K){
    dsk <- dsR(k,theta, eta)
    for(l in 1:K){
      dsl <- dsR(l,theta, eta)
      d2tkl <- d2sR(k,l,theta, eta)
      Hessloglik[k,l] <- sum((Ot * (d2tkl*st - dsk*dsl)*(st^(-2)))-d2tkl)
    }
  }
  return(Hessloglik)
}

# ------ gradient and hessians

#Check gradient and Hessian for R
# theta <- rnorm(K)
# eta <- runif(dim_eta)
# 
# round(grad(loglikR, theta, eta=eta),3)
# round(dloglikR(theta,eta),3)
# round(grad(loglikR, theta, eta=eta)-dloglikR(theta,eta))
# 
# round(hessian(loglikR, theta, eta=eta),3)
# round(d2loglikR(theta,eta),3)
# round(hessian(loglikR, theta, eta=eta)-d2loglikR(theta,eta),3)


#----- Laplace approximations
Laplace_approxR <- function(theta0,eta){
  
  dist <- 3
  epsilon <- 0.001
  Qlamb <- Q(eta[1])
  iter <- 0
  
  while(dist > epsilon){
    gradient <- dloglikR(theta0,eta=eta)
    Hess <-  d2loglikR(theta0,eta=eta)
    Sigma_new <- solve(Qlamb - Hess)
    theta_new <- as.numeric(Sigma_new %*% (gradient-Hess%*%theta0))
    lvec <- (theta_new-theta0)
    dist <- sqrt(sum((theta_new-theta0)^2))
    theta0 <- theta_new
    iter <- iter + 1
    #print(dist)
  }
  
  return(list(thetastar=theta0,
              Sigmastar=Sigma_new))
 
}

#----------------------- Hyperparameter posterior ------------------------------

# Function that back transforms the hyperparameters in original scale
invtransform <- function(x){
  backtransform <- c(exp(x[1]),exp(-exp(x[2:length(x)])))
  return(backtransform)
}

# From here on, eta is transformed to etaR where all parameters live in R

sR2 <- function(etaR, theta){
  
  Mmat <- matrix(0, nrow = Ttime, ncol = length(p))
  Mmat[1,1] <- Mtestim[1]
  for(t in 2:Ttime){
    Mmat[t,] <- c(Mtestim[t], Mmat[(t - 1), 1], Mmat[(t - 1), 2],
                  Mmat[(t - 1), 3], Mmat[(t - 1), 4])
  }
  Mmat2 <- rbind(rep(1, kmax), Mmat[1:(Ttime - 1), ])
  Mpvec <- rowSums(Mmat2 * pmat)
  
  st_theta <- c()
  Rtvec <- as.numeric(exp(B %*% theta))
  mutvec <- Rtvec * Mpvec # Vector of means mut
  Tmod7 <- floor(Ttime / 7)
  
  # Define composition matrix
  if (delay_struct == "oneday") {
    Cm <- C_oneday(invtransform(etaR)[3:dim_eta])
  } else if (delay_struct == "twoday") {
    Cm <- C_twoday(invtransform(etaR)[3:dim_eta])
  } else if (delay_struct == "weekend") {
    Cm <- C_weekend(invtransform(etaR)[3:dim_eta])
  }
  
  st_theta[1:7] <-  diag(apply(Cm * mutvec[1:7], 2, cumsum))
  
  if(Tmod7 == 1) {
    # if 8 <= Ttime <= 13 we have days 8,9,...,13
    Cm_augmented <- as.matrix(Cm[, 1:(Ttime %% 7)])
    mutmat <- matrix(0, nrow = 7, ncol = (Ttime %% 7))
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
    }
  } else if (Ttime %% 7 == 0) {
    # days 14,21,28,...
    Cm_augmented <- matrix(rep(Cm, Tmod7 - 1), ncol = (Tmod7 - 1) * 7)
    mutmat <- matrix(0, nrow = 7, ncol = (Tmod7 - 1) * 7)
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
    }
  } else if ((Tmod7 != 1) && (Ttime %% 7 != 0)) {
    # days > 14 and for which Ttime %% 7 != 0 
    Cm_augmented <- cbind(matrix(rep(Cm, Tmod7 - 1), ncol = (Tmod7 - 1) * 7),
                          Cm[, 1:(Ttime %% 7)])
    mutmat <- matrix(0, nrow = 7, ncol = ((Tmod7 - 1) * 7 + (Ttime %% 7)))
    for (t in 1:ncol(mutmat)) {
    mutmat[, t] <- moove(mutvec[(2 + t - 1):(2 + t - 1 + 6)], nmooves = t %% 7)
    }
  }
  
  st_theta[8:Ttime] <- colSums(Cm_augmented * mutmat)
  st_theta[which(st_theta <= 0)] <- 1e-8 # To avoid to input 0 in log()
  val <- exp(-exp(etaR[2])) * st_theta
  return(val)
}

logpetaR <- function(etaR,theta,Sigma){
  st <- sR2(etaR,theta)
  val <- sum(Ot * log(st) - st) -
    as.numeric(0.5 * t(theta)%*%(exp(etaR[1]) * P) %*% theta) + 
    0.5*log(det(Sigma))+
    (0.5 * K + a_lambda) * etaR[1] - b_lambda * exp(etaR[1]) +
    sum(etaR[2:dim_eta]-exp(etaR[2:dim_eta]))
  return(val)
} 


#-- (New) Gradient of logpetaR

is_even <- function(x){
  if(x%%2 == 0){
    even_number <- TRUE
  } else(even_number <- FALSE)
  return(even_number)
} # function to check if a number is even

dlogpetaR <- function(etaR, theta){

  Mmat <- matrix(0, nrow = Ttime, ncol = length(p))
  Mmat[1,1] <- Mtestim[1]
  for(t in 2:Ttime){
    Mmat[t,] <- c(Mtestim[t], Mmat[(t - 1), 1], Mmat[(t - 1), 2],
                  Mmat[(t - 1), 3], Mmat[(t - 1), 4])
  }
  Mmat2 <- rbind(rep(1, kmax), Mmat[1:(Ttime - 1), ])
  Mpvec <- rowSums(Mmat2 * pmat)
  
  Rtvec <- as.numeric(exp(B %*% theta))
  mutvec <- Rtvec * Mpvec # Vector of means mut
  st <- sR2(etaR, theta)

  # Partial derivative with respect to v
  deriv_v <- as.numeric((-0.5 * exp(etaR[1]) * (t(theta) %*% P %*% theta)) +
                          (0.5 * K + a_lambda) - b_lambda * exp(etaR[1]))

  # Partial derivative with respect to rho
  deriv_rho <- exp(etaR[2]) * sum(st - Ot) + (1 - exp(etaR[2]))

  # Partial derivative with respect to delta

  dOtst <- (-1) * diff(Ot / st)

  if(delay_struct=="oneday"){
    deriv_delta <- function(t_delta){

      daysvec <- seq(Ttime)

      if(t_delta == 1){
    idx <- daysvec[daysvec %% 7 == (t_delta + 5) | daysvec %% 7 == (t_delta-1)]
        eta_idx <- t_delta + 7
      }else if(t_delta == 0){
    idx <- daysvec[daysvec %% 7 == t_delta | daysvec %% 7 == (t_delta+1)][-1]
        eta_idx <- t_delta + 9
      }else{
        idx <- daysvec[daysvec %% 7 == (t_delta - 1) | daysvec %% 7 == t_delta]
        eta_idx <- t_delta + 1
      }
      nidx <- length(idx)
      idx_first <-  idx[seq(1, nidx, by = 2)]
      nidx_first <- length(idx_first)

      if(is_even(nidx)){
        ddelta <- exp(-exp(etaR[2])) * exp(etaR[eta_idx]- exp(etaR[eta_idx])) *
          sum(dOtst[idx_first] * mutvec[idx_first]) + (1-exp(etaR[eta_idx]))
      }else{
        ddelta <- exp(-exp(etaR[2])) * exp(etaR[eta_idx]-exp(etaR[eta_idx])) *
          (sum(dOtst[idx_first[-nidx_first]] * mutvec[idx_first[-nidx_first]]) +
    (Ot[tail(idx_first,1)]/st[tail(idx_first,1)]-1)*mutvec[tail(idx_first,1)])+
          (1-exp(etaR[eta_idx]))
      }

      return(ddelta)
    }
  }else if(delay_struct=="twoday"){
    dOtst2 <- (-1) * diff(Ot / st,lag = 2)

    deriv_delta <- function(t_delta){

      daysvec <- seq(Ttime)

      if(t_delta == 1){
    idx <- daysvec[daysvec %% 7 == (t_delta + 5) | daysvec %% 7 == (t_delta-1)]
    idx2 <- daysvec[daysvec %% 7 == (t_delta +5) | daysvec %% 7 == t_delta][-1]
    eta_idx <- 13; eta_idx2 <- 14
      }else if(t_delta == 0){
      idx <- daysvec[daysvec %% 7 == t_delta | daysvec %% 7 == (t_delta+1)][-1]
      idx2 <- daysvec[daysvec %% 7 == t_delta | daysvec %% 7 == t_delta+2][-1]
      eta_idx <- 15; eta_idx2 <- 16
      }else {
      idx <- daysvec[daysvec %% 7 == (t_delta - 1) | daysvec %% 7 == t_delta]
      idx2 <- daysvec[daysvec %% 7 == (t_delta - 1) | daysvec %% 7 == t_delta+1]
        if(t_delta == 2){
          eta_idx <- 3 ; eta_idx2 <- 4
        }else if (t_delta == 3){
          eta_idx <- 5; eta_idx2 <- 6
        }else if (t_delta == 4){
          eta_idx <- 7; eta_idx2 <- 8
        }else if (t_delta == 5){
          eta_idx <- 9; eta_idx2 <- 10
        }else if (t_delta == 6){
          idx2 <- daysvec[daysvec %% 7 == (t_delta - 1) | daysvec %% 7 == 0]
          eta_idx <- 11; eta_idx2 <- 12
        }
      }
      nidx <- length(idx)
      idx_first <-  idx[seq(1, nidx, by = 2)]
      nidx_first <- length(idx_first)
      nidx2 <- length(idx2)
      idx_first2 <-  idx2[seq(1, nidx2, by = 2)]
      nidx_first2 <- length(idx_first2)

      if(is_even(nidx)){
        ddelta <- exp(-exp(etaR[2])) * exp(etaR[eta_idx]- exp(etaR[eta_idx])) *
          sum(dOtst[idx_first] * mutvec[idx_first]) + (1-exp(etaR[eta_idx]))
      } else{
        ddelta <- exp(-exp(etaR[2])) * exp(etaR[eta_idx]-exp(etaR[eta_idx])) *
          (sum(dOtst[idx_first[-nidx_first]] * mutvec[idx_first[-nidx_first]]) +
    (Ot[tail(idx_first,1)]/st[tail(idx_first,1)]-1)*mutvec[tail(idx_first,1)])+
          (1-exp(etaR[eta_idx]))
      }

      if(is_even(nidx2)){
    ddelta2 <-  exp(-exp(etaR[2])) * exp(etaR[eta_idx2]- exp(etaR[eta_idx2])) *
          sum(dOtst2[idx_first2] * mutvec[idx_first2]) + (1-exp(etaR[eta_idx2]))
      }else{
       ddelta2 <- exp(-exp(etaR[2])) * exp(etaR[eta_idx2]-exp(etaR[eta_idx2])) *
    (sum(dOtst2[idx_first2[-nidx_first2]] * mutvec[idx_first2[-nidx_first2]]) +
  (Ot[tail(idx_first2,1)]/st[tail(idx_first2,1)]-1)*mutvec[tail(idx_first2,1)])+
      (1-exp(etaR[eta_idx2]))
      }

      ddelta_vec <- c(ddelta,ddelta2)
      return(ddelta_vec)
    }
  }else if(delay_struct=="weekend"){

    dOtst_lag3 <- (-1) * diff(Ot / st, lag = 3)

    deriv_delta <- function(t_delta){

      daysvec <- seq(Ttime)

      if (t_delta == 1 | t_delta == 0) {
        if(t_delta == 1){
          idx <- daysvec[daysvec %% 7 == t_delta |
                           daysvec %% 7 == (t_delta + 1)][-(1:2)]
          eta_idx <- t_delta + 7
          muS <- daysvec[daysvec %% 7 == 6] # Saturdays
        }else if(t_delta == 0){
          idx <- daysvec[daysvec %% 7 == t_delta + 1 |
                           daysvec %% 7 == (t_delta + 2)][-(1:2)]
          eta_idx <- t_delta + 9
          muS <- daysvec[daysvec %% 7 == 0] # Sundays
        }
        nidx <- length(idx)
        idx_first <-  idx[seq(1, nidx, by = 2)]
        nidx_first <- length(idx_first)


        if(is_even(nidx)){
        ddelta <- exp(-exp(etaR[2])) * exp(etaR[eta_idx]- exp(etaR[eta_idx])) *
            sum(dOtst[idx_first] * mutvec[muS]) + (1-exp(etaR[eta_idx]))
        }else{
          ddelta <- exp(-exp(etaR[2])) * exp(etaR[eta_idx]-exp(etaR[eta_idx])) *
            (sum(dOtst[idx_first[-nidx_first]] * mutvec[muS[-nidx_first]]) +
          (Ot[tail(idx_first,1)]/st[tail(idx_first,1)]-1)*mutvec[tail(muS,1)])+
            (1-exp(etaR[eta_idx]))
        }
      }else{
        if(t_delta == 6){
  idx <- daysvec[daysvec %% 7 == (t_delta-1) | daysvec %% 7 == (t_delta-5)][-1]
          eta_idx <- t_delta + 1
          dOtst <- dOtst_lag3
        }else if(t_delta == 2 | t_delta == 3 | t_delta == 4 | t_delta == 5){
      idx <- daysvec[daysvec %% 7 == (t_delta - 1) | daysvec %% 7 == t_delta]
          eta_idx <- t_delta + 1
        }
        nidx <- length(idx)
        idx_first <-  idx[seq(1, nidx, by = 2)]
        nidx_first <- length(idx_first)

        if(is_even(nidx)){
        ddelta <- exp(-exp(etaR[2])) * exp(etaR[eta_idx]- exp(etaR[eta_idx])) *
            sum(dOtst[idx_first] * mutvec[idx_first]) + (1-exp(etaR[eta_idx]))
        }else{
          ddelta <- exp(-exp(etaR[2])) * exp(etaR[eta_idx]-exp(etaR[eta_idx])) *
        (sum(dOtst[idx_first[-nidx_first]] * mutvec[idx_first[-nidx_first]]) +
    (Ot[tail(idx_first,1)]/st[tail(idx_first,1)]-1)*mutvec[tail(idx_first,1)])+
            (1-exp(etaR[eta_idx]))
        }
      }
      return(ddelta)
    }
  }

  if(delay_struct == "twoday"){
    dlogp_result <- c(deriv_v,deriv_rho,deriv_delta(2), deriv_delta(3),
                      deriv_delta(4),deriv_delta(5),deriv_delta(6),
                      deriv_delta(1),deriv_delta(0))
  }else{
  dlogp_result <- c(deriv_v,deriv_rho,sapply(c(seq(2,6),1,0),deriv_delta))
  }

  return(dlogp_result)
}

#---------- Newton-Raphson to find posterior max of hyperparam vector 


# Gradient ascent to determine "wise" initial condition
grad_ascent <- function(eta, dist=3, epsil=0.001){
  
  theta0 <- rep(1,K)  # Initial theta vector
  iter <- 0
  
  while(dist > epsil){
    
    eta_origin <- c(exp(eta[1]),exp(-exp(eta[2:dim_eta])))        
    Laplace <- Laplace_approxR(theta0 = theta0, eta = eta_origin) 
    theta_star <- Laplace$thetastar
    theta0 <- theta_star
    Sigma <- Laplace$Sigmastar
    
    gradp <- dlogpetaR(eta, theta_star)
    step  <- gradp*0.01
    eta_new <- eta+step
    dist <- sqrt(sum((eta_new-eta)^2))
    
    eta <- eta_new
    iter <- iter+1
    return(eta)
  }
  
}


Newton_Raphson <- function(eta0, epsil = 0.01, maxiter = 500){
  
  NRoutput <- matrix(0, nrow = maxiter + 1, ncol = dim_eta + 2)
  if(delay_struct=="oneday" | delay_struct=="weekend"){
    colnames(NRoutput) <- c("Iteration", 
                          c("lambda","rho",paste0(rep("p", 7),seq(1,7))),
                          "Distance")
  }else{
    colnames(NRoutput) <- c("Iteration", 
                            c("lambda","rho",paste0(rep("p", 14),seq(1,14))),
                            "Distance")
  }
  NRoutput[1, 2 : (dim_eta + 1)] <- c(exp(eta0[1]), exp(-exp(eta0[2:dim_eta])))
  NRoutput[1, (dim_eta + 2)] <- NA
  
  iter <- 0
  eta <- eta0
  theta0 <- rep(1,K) 
 
  for (j in 1:maxiter) {
    
    eta_origin <- c(exp(eta[1]),exp(-exp(eta[2:dim_eta])))
    Laplace <- Laplace_approxR(theta0, eta=eta_origin)
    theta_star <- Laplace$thetastar
    Sigma <- Laplace$Sigmastar
    theta0 <- theta_star
    gradf <-dlogpetaR(eta, theta=theta0)        
    Hessf <- hessian(logpetaR,eta,theta=theta0,Sigma=Sigma)  
    Hessf <- (-1)*(adjustPD(-Hessf)$PD)
    d_eta <- (-1) * solve(Hessf, gradf)  
    eta_new <- eta + d_eta            
    dist_temp <- sqrt(sum((eta_new - eta) ^ 2))
    
    if (dist_temp > 100) {
      eta_new <- eta + d_eta / (dist_temp)
    }
    
      step <- 1                                      # Initial step length
      iter.halving <- 1                              # Counter for step halving
      feta <- logpetaR(eta, theta0,Sigma)            # Function value at eta
      feta_new <- logpetaR(eta_new, theta0,Sigma)    # Function value at eta_new
      #Step halving
      while (feta_new <= feta) {
        step <- step * .5 
        eta_new <- eta + (step * d_eta)
        feta_new <- logpetaR(eta_new, theta0,Sigma)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30){break; print("Break")}
      }

    ############################################

    dist <- sqrt(sum((eta_new - eta) ^ 2))         
    dist_prev <- dist
    iter <- iter + 1
    eta <- eta_new
    NRoutput[j + 1, 1] <- iter
    NRoutput[j + 1, 2 : (dim_eta + 1)] <- c(exp(eta_new[1]),
                                            exp(-exp(eta_new[2:dim_eta])))
    NRoutput[j + 1, (dim_eta + 2)] <- dist
    
    
    if(dist < epsil){ # stop algorithm
      listout <- list(eta_optim = c(exp(eta_new[1]),
                                    exp(-exp(eta_new[2:dim_eta]))),
                      info = NRoutput[1 : (iter + 1), ],
                      converge = TRUE,
                      iter = iter,
                      distance = dist,
                      thetastar = theta_star,
                      Hessf = Hessf)
      attr(listout, "class") <- "NewtonRaphson"
      return(listout)
      break
    }
  }
  
  if(dist >= epsil){
    listout <- list(info = NRoutput,
                    converge = FALSE,
                    iter = iter)
    attr(listout, "class") <- "NewtonRaphson"
    return(listout)
  }
  
  
}

# Initial hyperparameter values

lambda_init <- 0.001
rho_init <- 0.7
prob_init <- rep(0.2, (dim_eta - 2))

if(delay_struct=="weekend"){
  lambda_init <- 2
  rho_init <- 0.75
  prob_init <- rep(0.45, (dim_eta - 2))
}

eta0<- c(log(lambda_init), log(-log(rho_init)), log(-log(prob_init)))
eta0 <- grad_ascent(eta0) 

NR <- Newton_Raphson(eta0)
if(NR$converge==FALSE){message("Warning NR did not converge")}

rho_hat <- NR$eta_optim[2]
delta_hat <- NR$eta_optim[3:dim_eta]
theta_hat <- Laplace_approxR(NR$thetastar,NR$eta_optim)$thetastar

# dlogpetaR(c(log(NR$eta_optim[1]),log(-log(NR$eta_optim[2:dim_eta]))),theta_hat)

#--- Credible intervals for Rt

Sighat <- Laplace_approxR(theta_hat, eta=NR$eta_optim)$Sigmastar

qz <- qnorm(alpha * .5, lower.tail=F)
gradg <- function(t){B[t,]}


gtheta_star <- as.numeric(B%*%theta_hat)
sd_postg <- c()
for(j in 1:Ttime){
  sd_postg[j] <- sqrt(as.numeric(t(gradg(j))%*%Sighat%*%(gradg(j))))
}

CIg_theta_low <- gtheta_star - qz * sd_postg
CIg_theta_up <-  gtheta_star + qz * sd_postg


CIlow <- exp(CIg_theta_low)
CIup <- exp(CIg_theta_up)


#------------------------ Estimation with EpiEstim


epiR <- estimate_R(Mtestim,
                   method="non_parametric_si",
                   config=make_config(list(si_distr=c(0,p))))
# plot(epiR,"R")


#------------------------ Estimation with Wallinga-Teunis method

Rwateu <- wallinga_teunis(Mtestim, 
                      method="non_parametric_si",
                      config=list(
                        t_start = seq(epiR$R$t_start[1],tail(epiR$R$t_start,1)), 
                        t_end = seq(epiR$R$t_end[1], tail(epiR$R$t_end,1)),
                                      si_distr = c(0,p),
                                      n_sim = 100))
 

Btdom <- cubicbs(seq(Ttime), lower=tl, upper=tr, K=K)$Bmatrix
Rthat <- exp(as.numeric(Btdom%*%theta_hat))[8:Ttime]

pdf(file="Fig7_Incidence.pdf", width = 12, height = 6)
plot(as.incidence(Ot))
dev.off()

if(delay_struct=="oneday") {
  pdf(file = "Fig7_AUS_oneday.pdf",
      width = 12,
      height = 6)
} else if (delay_struct == "twoday") {
  pdf(file = "Fig7_AUS_twoday.pdf",
      width = 12,
      height = 6)
}

plot(seq(8,Ttime),Rthat,col="black",type="l",xlab="Time",ylab="R",
     xlim=c(0,Ttime), ylim=c(0,2.7), xaxt="n",
     main=paste0("Estimated R ","(",delay_struct,")"),
     cex.main=0.92) #Estimated R with LPS
abline(h=1, lty=3,col="black")
axis(side = 1, at=seq(8,Ttime,by=10),
     labels = as.Date(dataAUS$date)[seq(8,Ttime,by=10)])

polygon(c(seq(8,Ttime),rev(seq(8,Ttime))),c(CIlow[8:Ttime],rev(CIup[8:Ttime])),
        col="azure3",border=NA)
lines(seq(8,Ttime),Rthat,col="blue",type="l",lwd=2)
lines(seq(8,Ttime),epiR$R$`Mean(R)`,col="black",type="l",lty=2, lwd=2)# Epiestim
lines(seq(8,Ttime),Rwateu$R$`Mean(R)`,lty=3,col="brown2",lwd=2)# Wallinga-Teunis
legend("topleft",c("LPS","EpiEstim","Wallinga-Teunis"),
       lty=c(1,2,3), lwd=c(2,2,2),
       col=c("blue","black","brown2"),bty="n",cex=0.8)

dev.off()

#--- (1-alpha) * 100% credible intervals for hyperparameters
Sigma_rhodelta <- round(solve(-NR$Hessf),5)
rhodelta_sd <- sqrt(diag(Sigma_rhodelta)[-1])
rhodelta_hat <- log(-log(c(rho_hat, delta_hat)))

#znorm <- qnorm(alpha * .5, lower = F)
dfS <- dim_eta-5 
tquant <- qt(alpha * .5, df=dfS,lower.tail=F)


CIrhodelta <- matrix(0, nrow = dim_eta - 1, ncol = 2)

#For rho
CIrhodelta[1, 1] <- exp(-exp(rhodelta_hat[1] +
                                 tquant * rhodelta_sd[1] * sqrt(dfS/(dfS-2))))
CIrhodelta[1, 2] <- exp(-exp(rhodelta_hat[1] -
                                  tquant * rhodelta_sd[1] * sqrt(dfS/(dfS-2))))


CIrhodelta[(2:(dim_eta-1)), 1] <- exp(-exp(rhodelta_hat[2:(dim_eta-1)] +
                      tquant * rhodelta_sd[2:(dim_eta-1)] * sqrt(dfS/(dfS-2))))
CIrhodelta[(2:(dim_eta-1)), 2] <- exp(-exp(rhodelta_hat[2:(dim_eta-1)] -
                      tquant * rhodelta_sd[2:(dim_eta-1)] * sqrt(dfS/(dfS-2))))



#--- Stop chronometer for timing
toc <- proc.time()-tic
toc_total <- toc_total + toc[3]

rho_hat







