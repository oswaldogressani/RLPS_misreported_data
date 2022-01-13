#------------------------------------------------------------------------------#
#                                                                              #
#                 Laplace-P-splines (LPS) for fast estimation                  #
#                        of the reproduction number                            #
#                         Plots of incidence data                              #
#                        © Oswaldo Gressani, 2022                              #
#                                                                              #
#------------------------------------------------------------------------------#

library("incidence")
library("ggplot2")
library("gridExtra")

set.seed(7589)

Ttime <- 30                              # Number of days
rho <- 0.2                               # True reporting rate
pdelay <- c(0.4,0.5,0.7,0.3,0.4,0.6,0.5) # True delay probabilities
p <- c(0.2,0.5,0.3)                      # Generation interval distribution
Rsmooth <- function(t) exp(cos(t / 13))  # True Rt

## One-day ---------------------------
source("simul_epidemic.R")
data <- simul_epidemic(Ttime, rho, p, delay="one", pij=pdelay, Rfunc=Rsmooth)

Mt <- data[,1] # Unobserved (true) number of cases
Ot <- data[,3] # Observed number of cases

incid_data <- as.incidence(as.data.frame(cbind(Ot,Mt)),dates = seq(1,Ttime))

plot1 <- plot(incid_data,stack = TRUE, border="grey",
              color=c("red","orange"))+ggtitle("One-day")+
                   theme(plot.title=element_text(hjust=0.5),
                         legend.title = element_blank())

## Two-day ---------------------------
pdelay <- c(0.3,0.4,0.4,0.3,0.3,0.4,0.5,0.3,0.3,0.4,0.5,0.3,0.3,0.6)
source("simul_epidemic.R")
data <- simul_epidemic(Ttime, rho, p, delay="two", pij=pdelay, Rfunc=Rsmooth)

Mt <- data[,1] # Unobserved (true) number of cases
Ot <- data[,3] # Observed number of cases

incid_data <- as.incidence(as.data.frame(cbind(Ot,Mt)),dates = seq(1,Ttime))

plot2 <- plot(incid_data,stack = TRUE, border="grey",
              color=c("red","orange"))+ggtitle("Two-day")+
                 theme(plot.title=element_text(hjust=0.5),
                       legend.title = element_blank())

## Weekend ---------------------------
pdelay <- c(0.4,0.5,0.7,0.3,0.4,0.6,0.5) # True delay probabilities
source("simul_epidemic.R")
data <- simul_epidemic(Ttime, rho, p, delay="weekend", pij=pdelay, 
                       Rfunc=Rsmooth)

Mt <- data[,1] # Unobserved (true) number of cases
Ot <- data[,3] # Observed number of cases

incid_data <- as.incidence(as.data.frame(cbind(Ot,Mt)),dates = seq(1,Ttime))

plot3 <- plot(incid_data,stack = TRUE, border="grey",
              color=c("red","orange"))+ggtitle("Weekend")+
                theme(plot.title=element_text(hjust=0.5),
                      legend.title = element_blank())

pdf(file="Figure1_Incidence_Scenario_1-2-3.pdf", width = 14, height = 5)

grid.arrange(plot1,plot2,plot3,ncol=3)

dev.off()



