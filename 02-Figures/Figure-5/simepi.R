#------------------------------------------------------------------------------#
#                                                                              #
#             Simulation of epidemic data for a model accounting for           #
#               underreporting and delay (following Azmon et al. 2013)         #    
#                 Copyright, Oswaldo Gressani. All rights reserved.            #
#                                                                              #
#------------------------------------------------------------------------------#

# The function simul_epidemic takes as input the number of days of an epidemic,
# the reporting probability (rho), the generation interval distribution (pgen),
# the delay structure (delay) and the true delay probabilities (pij). 

simul_epidemic <- function(Days, rho, pgen, delay=c("one","two","weekend"), 
                           pij, Rfunc){
  
# Days  --> Number of days (T) of an epidemic (an integer >= 8)
# rho   --> Reporting probability in (0,1)
# pgen  --> Generation interval distribution (a vector of proba)
# delay --> Delay structure "one", "two", "weekend"
# pij   --> Vector of true delay probabilities
# Rfunc --> True function of the reproduction number 

k <- length(pgen)                     # Maximum length of generation interval
Rsmooth <- Rfunc                      # Reproduction number as smooth function

#---- Generation of new (unobserved) disease counts (M) 
mut <- c()          # Mean of y
M <- c()            # Latent disease counts

for (t in 1:Days) {
  if (t == 1) {
    mut[t] <- 10   
    M[t] <- 10 # start with 10 index cases
  } else if (t >= 2 && t <= k) {
    mut[t] <- Rsmooth(t) *
      sum(rev(M[1:(k - 1)][1:(t - 1)]) *
            pgen[1:(k - 1)][1:(t - 1)])
    M[t] <- rpois(n = 1, lambda = mut[t])
  } else if (t > k && t <= Days) {
    mut[t] <- Rsmooth(t) * sum(rev(M[(t - k):(t - 1)]) * pgen)
    M[t] <- rpois(n = 1, lambda = mut[t])
  }
}

#---- Generation of reported disease counts (N)

N <- rpois(Days, rho * mut)

#---- Definition of composition matrix (C)

delay_structure <- match.arg(delay)

if(delay_structure == "one"){

dMTu <- pij[1]; dTuW <- pij[2]; dWTh <- pij[3] ; dThF <- pij[4]; dFSa <- pij[5]
dSaSu <- pij[6]; dSuM <- pij[7]
  C <- matrix(0, nrow = 7, ncol = 7)
    C[1,1] <- 1-dMTu; C[1,2] <- dMTu
      C[2,2] <- 1-dTuW; C[2,3] <- dTuW
        C[3,3] <- 1-dWTh; C[3,4] <- dWTh
          C[4,4] <- 1-dThF; C[4,5] <- dThF
            C[5,5] <- 1-dFSa; C[5,6] <- dFSa
              C[6,6] <- 1-dSaSu; C[6,7] <- dSaSu
                C[7,1] <- dSuM; C[7,7] <- 1-dSuM

}

if(delay_structure == "two"){
  dMTu <- pij[1]; dMW <- pij[2];     # Monday 
  dTuW <- pij[3]; dTuTh <- pij[4]    # Tuesday
  dWTh <- pij[5]; dWF <- pij[6]      # Wednesday
  dThF <- pij[7]; dThS <- pij[8]     # Thursday
  dFSa <- pij[9]; dFSu <- pij[10]    # Friday
  dSaSu <- pij[11]; dSaM <- pij[12]  # Saturday
  dSuM <- pij[13]; dSuT <- pij[14]   # Sunday
    C <- matrix(0, nrow = 7, ncol = 7)
    C[1,1] <- 1-dMTu-dMW; C[1,2] <- dMTu; C[1,3] <- dMW
      C[2,2] <- 1-dTuW-dTuTh; C[2,3] <- dTuW; C[2,4] <- dTuTh
       C[3,3] <- 1-dWTh-dWF; C[3,4] <- dWTh; C[3,5] <- dWF
        C[4,4] <- 1-dThF-dThS; C[4,5] <- dThF; C[4,6] <- dThS
         C[5,5] <- 1-dFSa-dFSu; C[5,6] <- dFSa; C[5,7] <- dFSu
           C[6,1] <- dSaM; C[6,6] <- 1-dSaSu-dSaM; C[6,7] <- dSaSu
            C[7,1] <- dSuM; C[7,2] <- dSuT; C[7,7] <- 1-dSuM-dSuT
 }

if(delay_structure == "weekend"){
  
dMTu <- pij[1]; dTuW <- pij[2]; dWTh <- pij[3]; dThF <- pij[4]; dFM <- pij[5]
dSaM <- pij[6]; dSuM <- pij[7]

C <- matrix(0, nrow = 7, ncol = 7)
  C[1,1] <- 1-dMTu; C[1,2] <- dMTu
    C[2,2] <- 1-dTuW; C[2,3] <- dTuW
     C[3,3] <- 1-dWTh; C[3,4] <- dWTh
      C[4,4] <- 1-dThF; C[4,5] <- dThF
       C[5,1] <- dFM; C[5,5] <- 1- dFM
        C[6,1] <- 1-dSaM; C[6,2] <- dSaM
         C[7,1] <- 1-dSuM; C[7,2] <- dSuM
  
}
  
#---- Generation of disease counts subject to underreporting and delay (O)

mutd <- rep(0, Days)
mutd[1:7] <- diag(apply(C * mut[1:7], 2, cumsum))

for(t in 8:Days){
  if(t%%7 == 1){
    mutd[t] <- as.numeric(C[, 1]%*%c(mut[t],mut[t-6],mut[t-5],mut[t-4],mut[t-3],
                                     mut[t-2],mut[t-1]))
  }
  if(t%%7 == 2){
    mutd[t] <- as.numeric(C[, 2]%*%c(mut[t-1],mut[t],mut[t-6],mut[t-5],mut[t-4],
                                     mut[t-3],mut[t-2]))
  }
  if(t%%7 == 3){
    mutd[t] <- as.numeric(C[, 3]%*%c(mut[t-2],mut[t-1],mut[t],mut[t-6],mut[t-5],
                                     mut[t-4],mut[t-3]))
  }
  if(t%%7 == 4){
    mutd[t] <- as.numeric(C[, 4]%*%c(mut[t-3],mut[t-2],mut[t-1],mut[t],mut[t-6],
                                     mut[t-5],mut[t-4]))
  }
  if(t%%7 == 5){
    mutd[t] <- as.numeric(C[, 5]%*%c(mut[t-4],mut[t-3],mut[t-2],mut[t-1],mut[t],
                                     mut[t-6],mut[t-5]))
  }
  if(t%%7 == 6){
    mutd[t] <- as.numeric(C[, 6]%*%c(mut[t-5],mut[t-4],mut[t-3],mut[t-2],mut[t-1],
                                     mut[t],mut[t-6]))
  }
  if(t%%7 == 0){
    mutd[t] <- as.numeric(C[, 7]%*%c(mut[t-6],mut[t-5],mut[t-4],mut[t-3],mut[t-2],
                                     mut[t-1],mut[t]))
  }
}  

O <- rpois(Days,(rho*mutd))


dataset <- cbind(M,N,O) 
return(dataset)

}





