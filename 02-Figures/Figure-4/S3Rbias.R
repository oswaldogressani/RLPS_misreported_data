#------------------------------------------------------------------------------#
#                                                                              #
#              Laplace-P-splines (LPS) for fast estimation                     #
#           of the reproduction number under misreported data                  #
#     Assessing the effect of the multiplication factor in Scenario 3          #
#                 Copyright, Oswaldo Gressani. All rights reserved.            #
#                                                                              #
#------------------------------------------------------------------------------#


# Load packages
library("blapsr")
library("EpiEstim")
library("numDeriv")
library("crayon")
library("progress")
library("reshape2")
library("ggplot2")
library("latex2exp")
source("simepi.R")


S3Rbias <- function(rho = 0.2, rhotilde = 0.2, S=20, seedval = 1647){

  S <- S                     # Number of replications
  delay_struct <- "weekend"  # Delay structure: "oneday", "twoday" or "weekend" 
  set.seed(seedval)          # Seed
  Ttime <- 30                # Number of days of an epidemic
  rho <- rho                 # True reporting probability
  K <- 10                    # Number of (cubic) B-splines in basis
  dimlat <- K                # Dimension of latent vector
  penorder <- 3              # Order of the penalty
  alpha <- 5/100             # Level for credible intervals
  tl <- 0                    # Lower bound of B-spline basis
  tr <- Ttime                # Upper bound of B-spline basis
  dim_eta <- 9
  pdelay <- c(0.4, 0.5, 0.7, 0.3, 0.4, 0.6, 0.5) # True delay probabilities
  NR_check <- matrix(0, nrow = S, ncol = dim_eta)
  NR_check_raw <-  matrix(0, nrow = S, ncol = dim_eta)
  rho_delta_hat <- matrix(0, nrow = S , ncol = dim_eta - 1)
  coverage_rho_delta <- matrix(0, nrow = S, ncol = dim_eta - 1)
  theta_mat <- matrix(0, nrow = S, ncol = K) 
  BiasRtmat <- matrix(0, nrow = S, ncol = Ttime)
  redraws <- 0
  
  # Influenza like generation interval distribution
  p <- c(0.233,0.359,0.198,0.103,0.053,0.027,0.014,0.007,0.003,0.002,0.001)  
  kmax <- length(p)                        # Max length of generation interval
  
  rhodelta_true <- c(rho,pdelay)           # True reporting proba and delays
  toc_total <- 0                           # Total time required for estimation
  CI_rhodelta_low <- matrix(0, nrow=S, ncol=(dim_eta-1))
  CI_rhodelta_up <- matrix(0, nrow=S, ncol=(dim_eta-1))
  
  #-- Progress bar
  progbar <- progress::progress_bar$new(
    format = crayon::white$yellow("Progressing [:elapsed :spin] [:bar] :percent"),
    total = S,
    clear = FALSE
  )
  
  #-------------------------- Loop starts
  for(s in 1:S){
    
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
    
    #------------------------- Generate the data -----------------------------------
    
    # Smooth function for the reproduction number
    Rsmooth <- function(t) exp(cos(t / 13))
    data <- simul_epidemic(Ttime, rho, p, delay=substr(delay_struct,1,3), 
                           pij=pdelay, Rfunc=Rsmooth)
    Ot <- data[, 3]
    
    while(Ot[1]==0){ # Resimulate data if 0 obs in first day
      data <- simul_epidemic(Ttime, rho, p, delay=substr(delay_struct,1,3), 
                             pij=pdelay, Rfunc=Rsmooth)
      Ot <- data[, 3]
      redraws <- redraws + 1
    }
    
    #-------------------------- LPS for fast estimation of Rt ----------------------
    
    tic <- proc.time()                               # Start tic-tac of the clock!
    B <- cubicbs(seq(Ttime), tl, tr, K = K)$Bmatrix  # B-spline basis 
    
    #--- Penalty matrix (P) and precision matrix (Q) of spline parameters
    
    D <- diag(K)
    for(j in 1:penorder){D <- diff(D)}
    P <- t(D) %*% D
    P <- P + diag(1e-05, K)                               # Penalty matrix
    Q <- function(lambda) lambda * P                      # Precision matrix
    a_lambda <- 1e-05; b_lambda <- 1e-05                  # Hyperpriors
    
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
    
    # Correcting for 0 obs on weekends
    Ott <- Ot
    Ott[6]  <- mean(Ot[1:5])
    Ott[7]  <- mean(Ot[1:6])
    Ott[13] <- mean(Ot[8:12])
    Ott[14] <- mean(Ot[8:13])
    Ott[20] <- mean(Ot[15:19])
    Ott[21] <- mean(Ot[15:20])
    Ott[27] <- mean(Ot[22:26])
    Ott[28] <- mean(Ot[22:27])
    Mtestim <- round((1 / rhotilde) * Ott)
    
    #----- Estimation of Rt and hyperparameter vector (lambda, rho, delta)
    
    sR <- function(theta, eta){
      
      Mpvec <- c()
      for (t in 1:Ttime) {
        if (t == 1) {
          Mpvec[t] <- Mtestim[1]
        } else if (t >= 2 && t <= kmax) {
          Mpvec[t] <- sum(rev(Mtestim[1:(kmax - 1)][1:(t - 1)]) * p[1:(kmax - 1)][1:(t - 1)])
        } else if (t > kmax && t <= Ttime) {
          Mpvec[t] <- sum(rev(Mtestim[(t - kmax):(t - 1)]) * p)
        }
      }
      
      st_theta <- c()
      Rtvec <- as.numeric(exp(B %*% theta))
      mutvec <- Rtvec * Mpvec 
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
      
      Mpvec <- c()
      for (t in 1:Ttime) {
        if (t == 1) {
          Mpvec[t] <- Mtestim[1]
        } else if (t >= 2 && t <= kmax) {
          Mpvec[t] <- sum(rev(Mtestim[1:(kmax - 1)][1:(t - 1)]) * p[1:(kmax - 1)][1:(t - 1)])
        } else if (t > kmax && t <= Ttime) {
          Mpvec[t] <- sum(rev(Mtestim[(t - kmax):(t - 1)]) * p)
        }
      }
      
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
      
      Mpvec <- c()
      for (t in 1:Ttime) {
        if (t == 1) {
          Mpvec[t] <- Mtestim[1]
        } else if (t >= 2 && t <= kmax) {
          Mpvec[t] <- sum(rev(Mtestim[1:(kmax - 1)][1:(t - 1)]) * p[1:(kmax - 1)][1:(t - 1)])
        } else if (t > kmax && t <= Ttime) {
          Mpvec[t] <- sum(rev(Mtestim[(t - kmax):(t - 1)]) * p)
        }
      }
      
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
    
    #----- Laplace approximations
    Laplace_approxR <- function(theta0,eta){
      
      dist <- 3
      epsilon <- 0.01
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
      }
      
      return(list(thetastar=theta0,
                  Sigmastar=Sigma_new))
      
    }
    
    #----------------------- Hyperparameter posterior -------------------------
    
    # Function that back transforms the hyperparameters in original scale
    invtransform <- function(x){
      backtransform <- c(exp(x[1]),exp(-exp(x[2:length(x)])))
      return(backtransform)
    }
    
    # Eta is transformed to etaR where all parameters live on real line
    
    sR2 <- function(etaR, theta){
      
      Mpvec <- c()
      for (t in 1:Ttime) {
        if (t == 1) {
          Mpvec[t] <- Mtestim[1]
        } else if (t >= 2 && t <= kmax) {
          Mpvec[t] <- sum(rev(Mtestim[1:(kmax - 1)][1:(t - 1)]) * p[1:(kmax - 1)][1:(t - 1)])
        } else if (t > kmax && t <= Ttime) {
          Mpvec[t] <- sum(rev(Mtestim[(t - kmax):(t - 1)]) * p)
        }
      }
      
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
      st_theta[which(st_theta <= 0)] <- 1e-8 
      val <- exp(-exp(etaR[2])) * st_theta
      return(val)
    }
    
    logpetaR <- function(etaR, theta, Sigma) {
      st <- sR2(etaR, theta)
      val <- sum(Ot * log(st) - st) -
        as.numeric(0.5 * t(theta) %*% (exp(etaR[1]) * P) %*% theta) +
        0.5 * log(det(Sigma)) +
        (0.5 * K + a_lambda) * etaR[1] - b_lambda * exp(etaR[1]) +
        sum(etaR[2:dim_eta] - exp(etaR[2:dim_eta]))
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
      
      Mpvec <- c()
      for (t in 1:Ttime) {
        if (t == 1) {
          Mpvec[t] <- Mtestim[1]
        } else if (t >= 2 && t <= kmax) {
          Mpvec[t] <- sum(rev(Mtestim[1:(kmax - 1)][1:(t - 1)]) * p[1:(kmax - 1)][1:(t - 1)])
        } else if (t > kmax && t <= Ttime) {
          Mpvec[t] <- sum(rev(Mtestim[(t - kmax):(t - 1)]) * p)
        }
      }
      
      Rtvec <- as.numeric(exp(B %*% theta))
      mutvec <- Rtvec * Mpvec 
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
    grad_ascent <- function(eta, dist=3, epsil=0.01){
      
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
          step <- step * .5 # Decrease step size
          eta_new <- eta + (step * d_eta)
          feta_new <- logpetaR(eta_new, theta0,Sigma)
          iter.halving <- iter.halving + 1
          if (iter.halving > 30){break; print("Break")}
        }
        dist <- sqrt(sum((eta_new - eta) ^ 2))  # Distance          
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
    prob_init <- rep(0.3, (dim_eta - 2))
    lambda_init <- 100
    rho_init <- rho + runif(1, min = 0, max = 0.2)
    if(rho == 1) {
      rho_init <- rho - runif(1, min = 0, max = 0.2)
    }
    if(rho == 1 && rhotilde == 0.8){
      rho_init <- 0.92
    }
    if(rho == 0.8 && rhotilde == 0.6){
      rho_init <- 0.78; prob_init <- rep(0.35, (dim_eta - 2))
    }
    if(rho == 0.8 && rhotilde == 0.4){
      rho_init <- 0.78; prob_init <- rep(0.35, (dim_eta - 2))
    }
    if(rho == 0.8 && rhotilde == 0.2){
      rho_init <- 0.78; prob_init <- rep(0.35, (dim_eta - 2))
    }
    
    eta0 <- c(log(lambda_init), log(-log(rho_init)), log(-log(prob_init)))
    eta0 <- grad_ascent(eta0) 
    
    NR <- Newton_Raphson(eta0)
    if(NR$converge==FALSE){message("Warning NR did not converge")}
    
    rho_hat <- NR$eta_optim[2]
    delta_hat <- NR$eta_optim[3:dim_eta]
    theta_hat <- Laplace_approxR(NR$thetastar,NR$eta_optim)$thetastar
    
    etaR_check <- c(log(NR$eta_optim[1]), log(-log(NR$eta_optim[2])),
                    log(-log(NR$eta_optim[3:dim_eta])))
    
    etaR_mat_check <- matrix(rep(etaR_check, 40), nrow = 40, byrow = T)
    etaR_mat_check[,2] <- seq(log(-log(NR$eta_optim[2]-0.0015)),
                              log(-log(NR$eta_optim[2]+0.0015)),
                              length=40)
    
    NR_checked <- matrix(0, nrow = 40, ncol = dim_eta)
    for(j in 1:40){
      NR_checked[j,] <- dlogpetaR(etaR_mat_check[j,],theta=theta_hat)
    }
    
    min_row <- which(abs(NR_checked[, 2]) == min(abs(NR_checked[, 2])))
    
    NR_check[s,] <- NR_checked[min_row,]
    NR_check_raw[s, ] <- dlogpetaR(etaR_check, theta_hat)
    
    #--- Stop chronometer for timing
    toc <- proc.time()-tic
    toc_total <- toc_total + toc[3]

    #--- Compute bias of Rt over days t=1,...,T.
    BiasRtmat[s, ] <- exp(as.numeric(B%*%theta_hat))-sapply(seq(Ttime),Rsmooth)
    progbar$tick()
    }
    BiasRt <- round(mean(colMeans(BiasRtmat)),3) #/ round(mean(rowMeans(BiasRtmat)),3)
    outlist <- list(BiasRt = BiasRt)
}

# ------------- Compute R biases under rho (mis)match --------------------------

# True rho = 0.2
S3_0202 <- S3Rbias(rho = 0.2, rhotilde = 0.2, S = 20)
S3_0204 <- S3Rbias(rho = 0.2, rhotilde = 0.4, S = 20)
S3_0206 <- S3Rbias(rho = 0.2, rhotilde = 0.6, S = 20)
S3_0208 <- S3Rbias(rho = 0.2, rhotilde = 0.8, S = 20)
S3_021  <- S3Rbias(rho = 0.2, rhotilde = 1.0, S = 20)

# True rho = 0.4
S3_0402 <- S3Rbias(rho = 0.4, rhotilde = 0.2, S = 20)
S3_0404 <- S3Rbias(rho = 0.4, rhotilde = 0.4, S = 20)
S3_0406 <- S3Rbias(rho = 0.4, rhotilde = 0.6, S = 20)
S3_0408 <- S3Rbias(rho = 0.4, rhotilde = 0.8, S = 20)
S3_041  <- S3Rbias(rho = 0.4, rhotilde = 1.0, S = 20)

# True rho = 0.6
S3_0602 <- S3Rbias(rho = 0.6, rhotilde = 0.2, S = 20)
S3_0604 <- S3Rbias(rho = 0.6, rhotilde = 0.4, S = 20)
S3_0606 <- S3Rbias(rho = 0.6, rhotilde = 0.6, S = 20)
S3_0608 <- S3Rbias(rho = 0.6, rhotilde = 0.8, S = 20)
S3_061  <- S3Rbias(rho = 0.6, rhotilde = 1.0, S = 20)

# True rho = 0.8
S3_0802 <- S3Rbias(rho = 0.8, rhotilde = 0.2, S = 20) #
S3_0804 <- S3Rbias(rho = 0.8, rhotilde = 0.4, S = 20) # ok
S3_0806 <- S3Rbias(rho = 0.8, rhotilde = 0.6, S = 20) # ok
S3_0808 <- S3Rbias(rho = 0.8, rhotilde = 0.8, S = 20) # ok
S3_081  <- S3Rbias(rho = 0.8, rhotilde = 1.0, S = 20) # ok

# Note: seedval needs to be the same for one block above (i.e. for 
# simulations S3_0802-S3_081). As such we ensure that the bias of Rt
# has really its source in the discrepancy between rho and rhotilde and
# not in the variability of the sampled values in the different scenarios
# within a block. This constraint is not required between blocks since
# there is variability among the sampled values anyway as different rho
# vlaues generate different data.

# True rho = 1
S3_102 <- S3Rbias(rho = 1, rhotilde = 0.2, S = 20)
S3_104 <- S3Rbias(rho = 1, rhotilde = 0.4, S = 20)
S3_106 <- S3Rbias(rho = 1, rhotilde = 0.6, S = 20)
S3_108 <- S3Rbias(rho = 1, rhotilde = 0.8, S = 20)
S3_11  <- S3Rbias(rho = 1, rhotilde = 1.0, S = 20)

# ------------- Compute Bias matrix --------------------------

RBiasMatrix <- matrix(0, nrow = 5, ncol = 5) # 5x5 matrix
rownames(RBiasMatrix) <- c("Rho=0.2","Rho=0.4","Rho=0.6","Rho=0.8","Rho=1")
colnames(RBiasMatrix) <- c("MF=0.2","MF=0.4","MF=0.6","MF=0.8","MF=1")

#-- Populate row 1
RBiasMatrix[1,1] <- S3_0202$BiasRt
RBiasMatrix[1,2] <- S3_0204$BiasRt
RBiasMatrix[1,3] <- S3_0206$BiasRt
RBiasMatrix[1,4] <- S3_0208$BiasRt
RBiasMatrix[1,5] <- S3_021$BiasRt

#-- Populate row 2
RBiasMatrix[2,1] <- S3_0402$BiasRt
RBiasMatrix[2,2] <- S3_0404$BiasRt
RBiasMatrix[2,3] <- S3_0406$BiasRt
RBiasMatrix[2,4] <- S3_0408$BiasRt
RBiasMatrix[2,5] <- S3_041$BiasRt

#-- Populate row 3
RBiasMatrix[3,1] <- S3_0602$BiasRt
RBiasMatrix[3,2] <- S3_0604$BiasRt
RBiasMatrix[3,3] <- S3_0606$BiasRt
RBiasMatrix[3,4] <- S3_0608$BiasRt
RBiasMatrix[3,5] <- S3_061$BiasRt

#-- Populate row 4
RBiasMatrix[4,1] <- S3_0802$BiasRt
RBiasMatrix[4,2] <- S3_0804$BiasRt
RBiasMatrix[4,3] <- S3_0806$BiasRt
RBiasMatrix[4,4] <- S3_0808$BiasRt
RBiasMatrix[4,5] <- S3_081$BiasRt

#-- Populate row 5
RBiasMatrix[5,1] <- S3_102$BiasRt
RBiasMatrix[5,2] <- S3_104$BiasRt
RBiasMatrix[5,3] <- S3_106$BiasRt
RBiasMatrix[5,4] <- S3_108$BiasRt
RBiasMatrix[5,5] <- S3_11$BiasRt

# ------------- Heatmap of Bias matrix --------------------------
RBiasMatrix2 <- (abs(RBiasMatrix))/max(RBiasMatrix)
RBiasMatrix2 <- round(RBiasMatrix2,3)


df <- melt(RBiasMatrix2)
colnames(df) <- c("x","y","value")
plotheat <- ggplot(df, aes(x=x,y=y,fill=value))+
  geom_tile(color = "white", lwd = 1.5) +
  scale_fill_gradient(low = "#55E268", high = "blue",
                      limits = c(0,1),
                      breaks = c(0,0.25,0.50,0.75,1),
                      labels = c(0,0.25,0.50,0.75,1)) +
  geom_text(aes(label = value), color = "white", size = 4) +
  coord_fixed() +
  ggtitle("Scenario 3") + 
  theme(axis.text.x = element_text(face="bold", color="black", size=11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "left",
        plot.title = element_text(size = 15, hjust = 0.5),
        legend.key.height = unit(1.5, "cm"),
        axis.text.y = element_text(face="bold", color="black", size=11)) +
  scale_x_discrete(labels=c("Rho=0.2" = TeX("$\\rho=0.2$"), 
                            "Rho=0.4" = TeX("$\\rho=0.4$"),
                            "Rho=0.6" = TeX("$\\rho=0.6$"),
                            "Rho=0.8" = TeX("$\\rho=0.8$"),
                            "Rho=1" = TeX("$\\rho=1$"))) + 
  scale_y_discrete(labels=c("MF=0.2" = TeX("$\\tilde{rho}=0.2$"), 
                            "MF=0.4" = TeX("$\\tilde{rho}=0.4$"),
                            "MF=0.6" = TeX("$\\tilde{rho}=0.6$"),
                            "MF=0.8" = TeX("$\\tilde{rho}=0.8$"),
                            "MF=1" = TeX("$\\tilde{rho}=1$")))

#--- Save as pdf
pdf(file="S3RBias_rho.pdf", width = 7, height = 7)
plotheat
dev.off()

#--- Save as png
# png(file="S3RBias_rho.png", width = 700, height = 700)
# plotheat
# dev.off()

#--- Save workspace
# save.image("S3Rbias.RData")



















