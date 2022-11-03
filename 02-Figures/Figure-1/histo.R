#------------------------------------------------------------------------------#
#                                                                              #
#                 Laplace-P-splines (LPS) for fast estimation                  #
#                        of the reproduction number                            #
#                   Plot of incidence data (Scenarios 1-3)                     #
#                 Copyright, Oswaldo Gressani. All rights reserved.            #
#                                                                              #
#------------------------------------------------------------------------------#

library("incidence")
library("ggplot2")
library("gridExtra")

set.seed(7589)

Ttime <- 30                              # Number of days
rho <- 0.2                               # True reporting rate

# Generation interval distribution (Influenza like generation interval)
p <- c(0.233,0.359,0.198,0.103,0.053,0.027,0.014,0.007,0.003,0.002,0.001) 
# Generation interval approximated with discr_si function of EpiEstim package
# si <- round(EpiEstim::discr_si(k = seq(1,11), mu = 2.6, sigma=1.5),3)
# si
# sum(si)
pdelay <- c(0.4,0.5,0.7,0.3,0.4,0.6,0.5) # True delay probabilities
Rsmooth <- function(t) exp(cos(t / 13))  # True Rt

## One-day ---------------------------
source("simepi.R")
data <- simul_epidemic(Ttime, rho, p, delay="one", pij=pdelay, Rfunc=Rsmooth)
while(data[,3][1]==0){ # Resimulate data if 0 obs in first day
  data <- simul_epidemic(Ttime, rho, p, delay="one", pij=pdelay, Rfunc=Rsmooth)
}

Mt <- data[,1] # Unobserved (true) number of cases
Ot <- data[,3] # Observed number of cases


incid_data <- as.incidence(as.data.frame(cbind(Ot,Mt)),dates = seq(1,Ttime))

plot1 <- plot(incid_data,stack = TRUE, border="grey",
              color=c("red","orange"))+ggtitle("One-day (Scenario 1) ")+
                   labs(x = "Time") +
                   theme(plot.title=element_text(hjust=0.5),
                         axis.title.x = element_text(size = 14),
                         axis.title.y = element_text(size = 14),
                         axis.text.x = element_text(size = 13),
                         axis.text.y = element_text(size = 14),
                         legend.text = element_text(size = 14),
                         legend.title = element_blank())

## Two-day ---------------------------
pdelay <- c(0.3,0.4,0.4,0.3,0.3,0.4,0.5,0.3,0.3,0.4,0.5,0.3,0.3,0.6)
source("simepi.R")
data <- simul_epidemic(Ttime, rho, p, delay="two", pij=pdelay, Rfunc=Rsmooth)
while(data[,3][1]==0){ # Resimulate data if 0 obs in first day
  data <- simul_epidemic(Ttime, rho, p, delay="two", pij=pdelay, Rfunc=Rsmooth)
}

Mt <- data[,1] # Unobserved (true) number of cases
Ot <- data[,3] # Observed number of cases

incid_data <- as.incidence(as.data.frame(cbind(Ot,Mt)),dates = seq(1,Ttime))

plot2 <- plot(incid_data,stack = TRUE, border="grey",
              color=c("red","orange"))+ggtitle("Two-day (Scenario 2)") +
              labs(x = "Time") +
                 theme(plot.title=element_text(hjust=0.5),
                       axis.title.x = element_text(size = 14),
                       axis.title.y = element_text(size = 14),
                       axis.text.x = element_text(size = 13),
                       axis.text.y = element_text(size = 14),
                       legend.text = element_text(size = 14),
                       legend.title = element_blank())

## Weekend ---------------------------
pdelay <- c(0.4,0.5,0.7,0.3,0.4,0.6,0.5) # True delay probabilities
source("simepi.R")
data <- simul_epidemic(Ttime, rho, p, delay="weekend", pij=pdelay, 
                       Rfunc=Rsmooth)
while(data[,3][1]==0){ # Resimulate data if 0 obs in first day
  data <- simul_epidemic(Ttime, rho, p, delay="weekend", pij=pdelay, 
                         Rfunc=Rsmooth)
}

Mt <- data[,1] # Unobserved (true) number of cases
Ot <- data[,3] # Observed number of cases

incid_data <- as.incidence(as.data.frame(cbind(Ot,Mt)),dates = seq(1,Ttime))

plot3 <- plot(incid_data,stack = TRUE, border="grey",
              color=c("red","orange"))+ggtitle("Weekend (Scenario 3)") +
              labs(x = "Time") +
                theme(plot.title=element_text(hjust=0.5),
                      axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14),
                      axis.text.x = element_text(size = 13),
                      axis.text.y = element_text(size = 14),
                      legend.text = element_text(size = 14),
                      legend.title = element_blank())

pdf(file="Fig1.pdf", width = 14, height = 5)

grid.arrange(plot1,plot2,plot3,ncol=3)

dev.off()



