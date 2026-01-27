#########################################
# Breast cancer - MNBcr/ MPIGcr - cure rate
#########################################
rm(list=ls(all=TRUE))
library(compiler);enableJIT(13)
library(pracma);library(asaur);library(PScr);library(readxl);library(devtools)
library(tidyverse);library(foreign);library(ggplot2);library(ggfortify);
library(survival);library(survminer);library(ggsurvfit);library(gt)
require(Bessel);require(VGAM)

setwd("G:/Meu Drive/MÁRCIA/Artigos/A new Poisson mixture cure rate model/programas/aplicacaoNOVA")

source("G:/Meu Drive/MÁRCIA/Artigos/A new Poisson mixture cure rate model/programas/aplicacaoNOVA/EM-MIG-MGA.R")


#===============================================================================
passoE<-function(psi,Data,dist,model)
#-------------------------------------------------------------------------------
{
  t     <- Data$t                  # tempo
  delta <- Data$delta              # censoring indicator
  x     <- matrix(Data$x,nrow=r)   # explanatory variables
  beta  <- matrix(psi[1:r],ncol=1) # vector of regression parameters 
  alpha <- psi[r+1]                # parameter for W
  nu    <- psi[r+2]                # parameter for W
  m     <- psi[r+3]                # parameter from  
  theta <- exp(t(x)%*%beta)        # regression estructure - parameter from Poisson (on mixture)
  #-----------------------------------------------------------------------------
  # survival function for W (dist)
  #-----------------------------------------------------------------------------
  if (dist == 1) { # Weibull
    S    <- exp(-exp(alpha) * t^nu)
  } else if (dist == 2) { # Log-normal
    S    <- plnorm(t, meanlog = alpha, sdlog = nu, lower.tail = FALSE)
  } else if (dist == 3) { # Gamma
    S    <- pgamma(t, shape = alpha, rate = nu, lower.tail = FALSE)
  } else if (dist == 4) { # Birnbaum-Saunders
    #S    <- pnorm((t - alpha) / (nu * sqrt(t)), lower.tail = FALSE)
    S    <-  pbisa(t, shape = alpha, scale = nu,  lower.tail = FALSE)
  } else {
    stop("Distribuição invalida. Escolha 1 (Weibull), 2 (Log-normal), 3 (Gamma) ou 4 (BS).")
  }
  #-----------------------------------------------------------------------------
  A<- B<- G <- H <- NULL
  #-MIG CASE (mix=1) -----------------------------------------------------------
  if(model==1){ 
    #--auxiliar values for conditional expectations-----------------------------
    p = delta - 0.5
    ay0 = 2*theta*(1-S) +  (1-m)/m
    ay1 = 2*theta*(1-S) +  (m+1)*(1-m)/(m^2)
    by0 = m*(1-m)
    by1 = (((m+1)^3)*(1-m))/m^2
    sab0    = sqrt(ay0*by0)
    sab1    = sqrt(ay1*by1)
    aux0_omega = (m^(5/2))*exp(1-m)*(besselK(sab0,p)/((ay0/by0)^(p/2)))
    aux1_omega = ((1-m)*exp(((m+1)*(1-m^2))/(m^2))*((m+1)^(3/2))*(besselK(sab1,p)/((ay1/by1)^(p/2))))
    #---------------------------------------------------------------------------
    omega = aux1_omega/(aux0_omega + aux1_omega)
    #--Conditional expectations-------------------------------------------------
    Y     = omega
    Z     = ((sqrt(by0)*besselK(sab0,p+1))/(sqrt(ay0)*besselK(sab0,p)))*(1-omega) + 
      ((sqrt(by1)*besselK(sab1,p+1))/(sqrt(ay1)*besselK(sab1,p)))*omega
    N     = delta + theta*S*Z 
    Z.0 =  (sqrt(by0)*besselK(sab0,p+1))/(sqrt(ay0)*besselK(sab0,p))
    Z.1 =  (((sqrt(by1)*besselK(sab1,p+1))/(sqrt(ay1)*besselK(sab1,p))))
    kappa.0 =  ((sqrt(ay0)*besselK(sab0,p+1))/(sqrt(by0)*besselK(sab0,p)) - ((2*p)/by0))
    kappa.1 =  ((sqrt(ay1)*besselK(sab1,p+1))/(sqrt(by1)*besselK(sab1,p)) - ((2*p)/by1))
    B = cbind(omega, Z.0, Z.1 , kappa.0, kappa.1)
    colnames(B)<-c("omega","Z.0","Z.1","kappa.0","kappa.1")
  }
  #-MGA CASE (mix=2)------------------------------------------------------------
  if(model==2){
    #--auxiliar values for conditional expectations-----------------------------
    aux0_nu = m*( (((1-m)^delta)*(((1-m)/m)^(1-m))) /( ((1-m)/m + theta*(1-S))^(1-m+delta) ) )
    aux1_nu = (1-m)*(( ((m+1)^delta) * ( ((1-m^2)/(m^2)) ^ (((m+1)*(1-m^2))/m^2 + delta) ))/
                       (( (1-m^2)/m^2 + theta*(1-S))^(((m+1)*(1-m^2))/m^2 + delta ) )  ) 
    #---------------------------------------------------------------------------
    nu =  aux1_nu/( aux0_nu +  aux1_nu)
    #--Conditional expectations-------------------------------------------------
    Y    = nu
    Z    = ((m*(1-m+delta))/(1+m*(theta*(1-S)-1)))*(1-nu) + ((((m+1)^2)*(1-m) + 
             (m^2)*delta)/(1+(m^2)*(theta*(1-S)-1)))*nu
    N    = delta + theta*S*Z 
    a.Z0=1-m+delta
    a.Z1=(1+m)^2*(1-m)/m^2+delta
    b.Z0=(1-m)/m+theta*(1-S)
    b.Z1=(1-m^2)/m^2+theta*(1-S)
    A<-cbind(nu,a.Z0,a.Z1,b.Z0,b.Z1)
    colnames(A)<-c("nu","a.Z0","a.Z1","b.Z0","b.Z1")
  }
  #--output-latent variables----------------------------------------------------
  return(list(Y = c(Y), N = c(N), Z = c(Z), A=A ,B=B) )
  #-------------------------------------------------------------------------------
} # end E-Step
#-------------------------------------------------------------------------------


#===============================================================================
#   Q-functions = Q1+Q2+Q3
#-------------------------------------------------------------------------------
Q1<-function(lambda,N,Data,dist)
  #-------------------------------------------------------------------------------
{
  t     <- Data$t
  delta <- Data$delta
  alpha <- lambda[1]
  nu    <- lambda[2]
  #-----------------------------------------------------------------------------
  # survival function and log density for W (dist)
  #-----------------------------------------------------------------------------
  if (dist == 1) { # Weibull
    S    <- exp(-exp(alpha) * t^nu)
    logf <- log(nu) + (nu - 1) * log(t) + alpha - exp(alpha) * t^nu
  } else if (dist == 2) { # Log-normal
    S    <- plnorm(t, meanlog = alpha, sdlog = nu, lower.tail = FALSE)
    logf <- dlnorm(t, meanlog = alpha, sdlog = nu, log = TRUE)
  } else if (dist == 3) { # Gamma
    S    <- pgamma(t, shape = alpha, rate = nu, lower.tail = FALSE)
    logf <- dgamma(t, shape = alpha, rate = nu, log = TRUE)
  } else {
    stop("Distribuição invalida. Escolha 1 (Weibull), 2 (Log-normal), 3 (Gamma) ou 4 (BS).")
  }
  #-output---------------------------------------------------------------------
  q1_function <- -sum((N - delta) * log(S) + delta * logf)
  return(q1_function)
  #------------------------
} # end Q1
#-------------------------------------------------------------------------------
Q2<-function(beta,N,Z,Data)
  #-------------------------------------------------------------------------------
{
  x     <- matrix(Data$x,nrow=r)    
  beta  <- matrix(beta,ncol=1)  
  #-----------------------------------------------------------------------------
  theta     <- exp(t(x)%*%beta)
  log.theta <- log(theta)
  #-output--Equation 8----------------------------------------------------------
  q2_function <- -sum( N*log.theta - Z*theta)
  return(q2_function)
  #------------------------
} # end Q2 
#-------------------------------------------------------------------------------
Q3 <- function(m, Y, A, B,model) {
  #-------------------------------------------------------------------------------
  n <- length(Y)
  #-MIG case---------------------------------------------------------------------- 
  if (model == 1) {
    
    q3_function <- (3/2)*(log(m)*(1-B[, "omega"]) +log(m+1)*B[, "omega"] ) +
      (1/2)*log((1-m)/m^2) - ((1-m)/(2*m^2))*(m*B[,"Z.0"]*(1-B[, "omega"]) + (m+1)*B[,"Z.1"]*B[, "omega"] -
        2*((m^2)*(1-B[, "omega"]) + ((m+1)^2)*B[, "omega"])  +
         (m^3)*B[,"kappa.0"]*(1-B[, "omega"]) +((m+1)^3)*B[,"kappa.1"]*B[, "omega"]     ) +
      Y*log(1-m)+(1-Y)*log(m)
    q3_function=-sum(q3_function)
    return(q3_function)
    #-MGA Case----------------------------------------------------------------------
  } else if (model == 2) {
    q3_function <- ((1 - m) / m^2) * (log((1 - m) / m^2) *(m^2*(1-A[, "nu"])+(m+1)^2*A[, "nu"])+
                    m^2*log(m)*(1-A[, "nu"])+(m+1)^2*log(m+1)*A[, "nu"]+
                     m^2*(digamma(A[, "a.Z0"])-log(A[, "b.Z0"]))*(1-A[, "nu"]) +
                      (m+1)^2*(digamma(A[, "a.Z1"])-log(A[, "b.Z1"]))*A[, "nu"] -
                      m*(A[, "a.Z0"]/A[, "b.Z0"])*(1-A[, "nu"])-(m+1)*(A[, "a.Z1"]/A[, "b.Z1"])*A[, "nu"])-
      lgamma(1 - m) * (1-A[, "nu"]) - lgamma((m+1)^2*(1-m)/m^2) * (A[, "nu"])+
      Y*log(1-m)+(1-Y)*log(m)
    q3_function=-sum(q3_function)
    return(q3_function)
  } else {
    stop("Modelo invalido: use 1 (MIG) ou 2 (MGA).")
  }
  #------------------------
} # end  Q3




#===============================================================================
llikeobserved<-function(psi,Data,dist,model,real=FALSE)
#-------------------------------------------------------------------------------
{
  t     <- Data$t
  delta <- Data$delta
  x     <- matrix(Data$x,nrow=r)
  beta  <- matrix(psi[1:r],ncol=1)
  alpha <- psi[r+1]
  nu    <- psi[r+2]
  m     <- psi[r+3]
  theta <- exp(t(x)%*%beta)  #parameter from Poisson (on mixture)
  #-----------------------------------------------------------------------------
  # survival function and log density for W (dist)
  #-----------------------------------------------------------------------------
  if (dist == 1) { # Weibull
    S    <- exp(-exp(alpha) * t^nu)
    logf <- log(nu) + (nu - 1) * log(t) + alpha - exp(alpha) * t^nu
  } else if (dist == 2) { # Log-normal
    S    <- plnorm(t, meanlog = alpha, sdlog = nu, lower.tail = FALSE)
    logf <- dlnorm(t, meanlog = alpha, sdlog = nu, log = TRUE)
  } else if (dist == 3) { # Gamma
    S    <- pgamma(t, shape = alpha, rate = nu, lower.tail = FALSE)
    logf <- dgamma(t, shape = alpha, rate = nu, log = TRUE)
  }else {
    stop("Distribuição invalida. Escolha 1 (Weibull), 2 (Log-normal), 3 (Gamma) ou 4 (BS).")
  }
  #--MIG case-------------------------------------------------------------------
  if (model == 1) {
    aux1_spop = m * exp((1 - m) * (1 - sqrt(1 - (2 * m * theta*(S-1)) / (1 - m))))
    aux2_spop = (1 - m) * exp(((1 + m) * (1 - m^2)) / m^2 * (1 - sqrt(1 - (2 * m^2 * theta*(S-1)) / (1 - m^2))))
    #---------------------------------------------------------------------------
    aux1_fpop = sqrt(1 + (2*m*theta*(1-S))/(1-m))
    aux2_fpop = sqrt(1 + (2*(m^2)*theta*(1-S))/(1-m^2))
    #---------------------------------------------------------------------------
    log.Spop = log(aux1_spop + aux2_spop)
    log.fpop = log(theta) + logf + log( ( (m^2)/aux1_fpop ) * exp( (1-m)*(1-aux1_fpop) )  +  
                                          ( (1-m^2)/aux2_fpop ) * exp( (((m+1)*(1-m^2))/(m^2))*(1-aux2_fpop) ) )
    #-output--------------------------------------------------------------------
    loglik_obs = - sum(delta*log.fpop +(1-delta)*log.Spop)
    return(loglik_obs)
    #-MGA case--------------------------------------------------------------------
  } else if (model == 2) {
    aux1_spop = 1- (m*theta*(S-1)) / (1-m) 
    aux2_spop = 1- (m^2*theta*(S-1)) / (1-m^2)
    #---------------------------------------------------------------------------
    aux1_fpop =  1 + (m*theta*(1-S)) / (1-m) 
    aux2_fpop =  1 + (m^2*theta*(1-S)) / (1-m^2)
    #---------------------------------------------------------------------------
    log.Spop = log( m*(aux1_spop)^(-(1-m)) + (1-m)*(aux2_spop)^(-((m+1)*(1-m^2))/m^2))
    log.fpop = log(theta)+logf+ log( (m^2)*(aux1_spop)^(-(1-m)-1) + 
                                       (1-m^2)*(aux2_spop)^(-((m+1)*(1-m^2))/m^2 -1))
    #-output----------------------------------------------------------------------
    loglik_obs =  - sum(delta*log.fpop +(1-delta)*log.Spop)
    return(loglik_obs)
  } else {
    stop("Modelo invalido: use 1 (MIG) ou 2 (MGA).")
  }
  #------------------------
} # end llikeobserved








#database-----------------------------------------
x          <- read.dbf("pacigeral.dbf",as.is=FALSE)
#------------------------------------------------
## C50 --> cancer de mama
a2   <- x %>% 
  filter(TOPOGRUP=="C50") %>% 
  filter(DTDIAG >= "2008-12-31" & DTDIAG <="2016-12-31") %>% 
  #filter(DTULTINFO <= "2016-12-31") %>% 
  mutate(tempo=as.numeric(difftime(DTULTINFO,DTDIAG,units="days"))/365) %>% 
  filter(tempo > 0) %>%
  filter(IDADE !=0) %>%
  filter(IDADE !=1) %>%
  filter(SEXO != 1)%>%
  mutate(status=ifelse(ULTINFO==3,1,0)) %>% 
  mutate(ECGRUP  = as.character(ECGRUP)) %>% 
  filter(ECGRUP != "X") %>%
  filter(ECGRUP != "Y") %>% 
  filter(ECGRUP != "0") %>% 
  mutate(ECGRUP  = as.factor(ECGRUP))

#-------------------------------------------
# MODELING
#-------------------------------------------
# model matrix 
#-------------------------------------------
# Stage of disease + Individualized Treatments + Age 
x  <- t(model.matrix(~ as.factor(a2$ECGRUP)+as.factor(a2$CIRURGIA)+ 
                       as.factor(a2$RADIO)+as.factor(a2$QUIMIO)+ a2$IDADE   ))
#-------------------------------------------
#Data
#-------------------------------------------
t       <- a2$tempo
delta   <- a2$status
n       <- length(t) # Sample size
r <- nrow(x)
Data <- list(t=t,delta=delta,x=x)




#------------------------------
maximo.logMIG = matrix(NA,3,1) 
AIC_MIG       = matrix(NA,3,1)  
BIC_MIG       = matrix(NA,3,1)  
maximo.logMGA = matrix(NA,3,1)  
AIC_MGA       = matrix(NA,3,1)  
BIC_MGA       = matrix(NA,3,1) 
#------------------------------
# distribution for lifetime 
#------------------------------
# dist=1 #Weibull
# dist=2 #log-normal
# dist=3 #gamma
# dist=4 #BS
#------------------------------------------------------------
for(dist in c(1,2,3)){

#------------------------------
# MIXTURE MODEL - MIG
#------------------------------
model = 1
i <- 1;  dif <- 10;  lower <- c(-Inf,0)
#if(dist==3 || dist==4){lower=c(0,0)}
if(dist==1 || dist==2){alpha.ini   <- -log(mean(t))}
if(dist==3 || dist==4){alpha.ini   <- log(mean(t))} 
if(dist==3 || dist==4){lower=c(1e-5, 1e-5)}

#-------------------------------------------
#-- for initial values of beta0
KM = survfit(Surv(t, delta) ~ 1)$surv
beta0.ini = KM[length(KM)]
#-------------------------------------------
#initial values for parameter estimation
beta.ini    <- c(beta0.ini,rep(0,r-1));beta.ini
#alpha.ini   <- log(mean(t))
nu.ini      <- 1
m.ini       <- 0.5
# inicial values for vector of parameters
psi.ini     <- c(beta.ini,alpha.ini,nu.ini,m.ini) 
psi         <- psi.ini

#-------------------------------------------
#  set.seed(3072022)
#------------------------------------------- 
# EM algorithm
#-------------------------------------------
while(i<=10000 && dif>1e-4)
{
  #-Passo E-----------------------------------
  latentes  = passoE(psi,Data,dist,model)
  #-------------------------------------------
  Y         = latentes$Y
  Z         = latentes$Z
  N         = latentes$N
  A         = latentes$A
  B         = latentes$B
  #-PASSO M----------------------------------- 
  beta      = psi[1:r]
  lambda    = psi[(r+1):(r+2)]
  m         = psi[r+3]
  #-------------------------------------------
  maximo1   = optim(lambda,Q1,method = "L-BFGS-B",lower=lower,N=N,Data=Data,dist=dist,hessian=F)
  lambda    = maximo1$par  # estimando lambda
  
  #-------------------------------------------
  maximo2   = optim(beta,Q2,method = "BFGS",N=N,Z=Z,Data=Data,hessian=F)
  beta      = maximo2$par #estimando beta
  
  #------------------------------------------- 
  maximo3   = optim(m,Q3,method = "Brent",lower=0.0000001,upper=0.9999999,Y=Y,A=A,B=B,model=model,hessian=F)
  m         = maximo3$par  #estimando m
  
  #--guardadndo estimativas num vetor auxiliar----- 
  psi.aux   = c(beta,lambda, m) 
  #----criterio EM---------------------------------
  dif       = max(abs(psi-psi.aux))  
  #--atualizacao psi-------------------------------
  psi       = psi.aux #atualiza psi
  #--próxima iteracao-----------------------------
  i=i+1
  cat("iteracao=", i, "\n")
}# end EM

#-------------------------------------------
psi.MIG       = matrix(psi,nrow=1);colnames(psi.MIG) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");round(psi.MIG,7)

#------------------------------------
# Criterios AIC e BIC p/ o Mix MIG
maximo.logMIG[dist,]  = llikeobserved(psi.MIG,Data=Data,dist =dist,model=1);maximo.logMIG 
AIC_MIG[dist,]        = 2*(maximo.logMIG[dist,]) + 2*length(psi.MIG);AIC_MIG
BIC_MIG[dist,]        = 2*(maximo.logMIG[dist,]) + log(length(t))*length(psi.MIG);BIC_MIG

#==========================================================================================
# MODELO MGA
#------------------------------------
model = 2
i <- 1;  dif <- 10;  lower <- c(-Inf,0)
if(dist==1 || dist==2){alpha.ini   <- -log(mean(t))}
if(dist==3 || dist==4){alpha.ini   <- log(mean(t))} 
if(dist==3 || dist==4){lower=c(1e-5, 1e-5)}

#-------------------------------------------
#-- for initial values of beta0
KM = survfit(Surv(t, delta) ~ 1)$surv
beta0.ini = KM[length(KM)]
#-------------------------------------------
#initial values for parameter estimation
beta.ini    <- c(beta0.ini,rep(0,r-1));beta.ini
nu.ini      <- 1
m.ini       <- 0.5
# inicial values for vector of parameters
psi.ini     <- c(beta.ini,alpha.ini,nu.ini,m.ini) 
psi         <- psi.ini
#-------------------------------------------
#  set.seed(3072022)
#------------------------------------------- 
# EM algorithm
#-------------------------------------------
while(i<=10000 && dif>1e-4)
{
  #-Passo E-----------------------------------
  latentes  = passoE(psi,Data,dist,model)
  #-------------------------------------------
  Y         = latentes$Y
  Z         = latentes$Z
  N         = latentes$N
  A         = latentes$A
  B         = latentes$B
  #-PASSO M----------------------------------- 
  beta      = psi[1:r]
  lambda    = psi[(r+1):(r+2)]
  m         = psi[r+3]
  #-------------------------------------------
  maximo1   = optim(lambda,Q1,method = "L-BFGS-B",lower=lower,N=N,Data=Data,dist=dist,hessian=F)
  lambda    = maximo1$par  # estimando lambda
  
  #------------------------------------------- 
  maximo2   = optim(beta,Q2,method = "BFGS",N=N,Z=Z,Data=Data,hessian=F)
  beta      = maximo2$par #estimando beta
  
  #------------------------------------------- 
  maximo3   = optim(m,Q3,method = "Brent",lower=0.0000001,upper=0.9999999,Y=Y,A=A,B=B,model=model,hessian=F)
  m         = maximo3$par  #estimando m
  
  #--guardadndo estimativas num vetor auxiliar----- 
  psi.aux   = c(beta,lambda, m) 
  #----criterio EM---------------------------------
  dif       = max(abs(psi-psi.aux))  
  #--atualizacao psi-------------------------------
  psi       = psi.aux #atualiza psi
  #--próxima iteracao-----------------------------
  i=i+1
  cat("iteracao=", i, "\n")
  #cat("replica=", jj, "\n")
}# end EM

#-------------------------------------------
psi.MGA      = matrix(psi,nrow=1);colnames(psi.MGA) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");round(psi,7)
#------------------------------------
# Criterios AIC e BIC p/ o Mix MGA
maximo.logMGA[dist,]  = llikeobserved(psi.MGA,Data=Data,dist =dist,model=model);maximo.logMGA 
AIC_MGA[dist,]        = 2*(maximo.logMGA[dist,]) + 2*length(psi.MGA);AIC_MGA
BIC_MGA[dist,]        = 2*(maximo.logMGA[dist,]) + log(length(t))*length(psi.MGA);BIC_MGA
 
}

#------------------------------------------------------------




#===================================================================== 
#COMPARING MODEL - BIN NEG, POISSON e BERNOULLI - Package (PScr)  
library(PScr)
#===================================================================== 
# model = distribution to be used for the concurrent causes: 
# 1 for Poisson,
# 3 for negative binomial, 
# 4 for bernoulli 
#------------------------------------------------------------
# Distribution for lifetime
# dist= distribution to be used for the time-to-event:
# 2 for Weibull, 
# 3 for gamma and 
# 5 log-normal
#----------
llikeEstimated_Bern <- matrix(NA, nrow = 3, ncol = 1)
AIC_Bern             <- matrix(NA, nrow = 3, ncol = 1)
BIC_Bern             <- matrix(NA, nrow = 3, ncol = 1)
#----------
llikeEstimated_PO    <- matrix(NA, nrow = 3, ncol = 1)
AIC_PO               <- matrix(NA, nrow = 3, ncol = 1)
BIC_PO               <- matrix(NA, nrow = 3, ncol = 1)
#----------
llikeEstimated_BN    <- matrix(NA, nrow = 3, ncol = 1)
AIC_BN               <- matrix(NA, nrow = 3, ncol = 1)
BIC_BN               <- matrix(NA, nrow = 3, ncol = 1)
#---------- 
dist_values <- c(2, 3, 5)
for (j in 1:3) {
  dist <- dist_values[j]
  
  # ----------------- Bernoulli -----------------
  EM_Bern  <- EM.PScr(t, delta, x, model = 4, dist = dist, max.iter = 10000)
  llikeEstimated_Bern[j,] <- EM_Bern$loglike
  AIC_Bern[j,]             <- EM_Bern$AIC
  BIC_Bern[j,]             <- EM_Bern$BIC
  print(paste("model-bern, dist =", dist))
  
  # ----------------- Poisson -------------------
  EM_PO <- EM.PScr(t, delta, x, model = 1, dist = dist, max.iter = 10000)
  llikeEstimated_PO[j,] <- EM_PO$loglike
  AIC_PO[j,]            <- EM_PO$AIC
  BIC_PO[j,]            <- EM_PO$BIC
  print(paste("model-poisson, dist =", dist))
  
  # ------------- Negative Binomial -------------
  EM_BN <- EM.PScr(t, delta, x, model = 3, dist = dist, max.iter = 10000)
  llikeEstimated_BN[j,] <- EM_BN$loglike
  AIC_BN[j,]            <- EM_BN$AIC
  BIC_BN[j,]            <- EM_BN$BIC
  print(paste("model-BN, dist =", dist))
}
#==========================================================================================
# Selection Criteria for MPIGcr, MNBcr, NBcr, PTcr and BERcr
#------------------------------------
#dist Weibull
AIC_1 = round(c(AIC_MIG[1], AIC_MGA[1], AIC_BN[1], AIC_PO[1], AIC_Bern[1]),2)
#dist Log-normal
AIC_2 = round(c(AIC_MIG[2],AIC_MGA[2],AIC_BN[3],AIC_PO[3],AIC_Bern[3]),2)
#dist Gamma
AIC_3 = round(c(AIC_MIG[3],AIC_MGA[3],AIC_BN[2],AIC_PO[2],AIC_Bern[2]),2)

#------------------------------------
#dist Weibull
BIC_1 = round(c(BIC_MIG[1],BIC_MGA[1],BIC_BN[1],BIC_PO[1],BIC_Bern[1]),2)
#dist Log-normal
BIC_2 = round(c(BIC_MIG[2],BIC_MGA[2],BIC_BN[3],BIC_PO[3],BIC_Bern[3]),2)
#dist Gamma
BIC_3 = round(c(BIC_MIG[3],BIC_MGA[3],BIC_BN[2],BIC_PO[2],BIC_Bern[2]),2)

#------------------------------------
#dist Weibull
llike_1 = round(c(maximo.logMIG[1],maximo.logMGA[1],llikeEstimated_BN[1],llikeEstimated_PO[1],llikeEstimated_Bern[1]),2)
#dist Log-normal
llike_2 = round(c(maximo.logMIG[2],maximo.logMGA[2],llikeEstimated_BN[3],llikeEstimated_PO[3],llikeEstimated_Bern[3]),2)
#dist Gamma
llike_3 = round(c(maximo.logMIG[3],maximo.logMGA[3],llikeEstimated_BN[2],llikeEstimated_PO[2],llikeEstimated_Bern[2]),2)





#===================================================
# Breast cancer  
# Table 4: AIC and BIC obtained by fitting the MPIGcr, MNBcr, Negative Binomial 
# cure rate (NBcr), PTcr, and Bernoulli cure rate (BERcr) models to the breast cancer dataset.
#===================================================

AIC =  matrix(c(AIC_1,AIC_2,AIC_3), 3, 5, byrow = TRUE) 
BIC =  matrix(c(BIC_1,BIC_2,BIC_3), 3, 5, byrow = TRUE) 
llikeEstimated = matrix(c(llike_1,llike_2,llike_3), 3, 5, byrow = TRUE)

criterios = rbind(AIC,BIC,llikeEstimated )
 criterios_df = data.frame(
  Criteria = rep(c("AIC", "BIC", "log-like"), each = 3),
`Time distribution` = rep(c("Weibull", "Log-normal", "Gamma"), times = 3),
  MIG = criterios[, 1],
  MGA = criterios[, 2],
  BN  = criterios[, 3],
  POI = criterios[, 4],
  BER = criterios[, 5]
)
 
criterios_df %>%
  gt() %>%
  tab_spanner(
    label = "Model",
    columns = c(MIG, MGA, BN, POI, BER)
  ) %>%
  tab_header(
    title = "Criteria, Time distribution and Models"
  )


name1 = paste("Table-Selection_criteria")
write.table(criterios_df,file=paste(name1,".txt"))
write.csv(criterios_df,file=paste(name1,".csv"))



#===================================================
# saving external file to be used 
save(list = ls(),file = "Aplication-Table-Selection_criteria.RData")
#===================================================
