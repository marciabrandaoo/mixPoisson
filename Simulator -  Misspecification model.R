#===================================================
# Simulation Study
# Table 3: Differences in AIC and BIC between models fitted with MNBcr, MPIGcr and 
# PTcr competing causes, and  empirical rejection rate under H0 of the likelihood ratio test, 
# for data generated with PTcr competing causes under different sample sizes
#===================================================
rm(list=ls(all=TRUE))
library(compiler);enableJIT(13)
library(survival)
library(pracma)
library(PScr)



setwd("G:/Meu Drive/MÁRCIA/Artigos/A new Poisson mixture cure rate model/programas/programas/sim_misspecification")
set.seed(104729)

#------------------------
# sample size

sample_size = c(200,300,500,1000)
num_n = length(sample_size)
 
media_mga    = matrix(NA,num_n,11)
mediana_mga  = matrix(NA,num_n,11)
SD_mga       = matrix(NA,num_n,11)
bias_mga     = matrix(NA,num_n,11)
rel_bias_mga = matrix(NA,num_n,11)
SE_mga       = matrix(NA,num_n,11)
rmse_mga     = matrix(NA,num_n,11)
cover_mga    = matrix(NA,num_n,11)

media_poi    = matrix(NA,num_n,11)
mediana_poi  = matrix(NA,num_n,11)
SD_poi       = matrix(NA,num_n,11)
bias_poi     = matrix(NA,num_n,11)
rel_bias_poi = matrix(NA,num_n,11)
SE_poi       = matrix(NA,num_n,11)
rmse_poi     = matrix(NA,num_n,11)
cover_poi    = matrix(NA,num_n,11)
#----------------------------------------
mean_test_statistic_MIG = matrix(NA,num_n,1)
mean_difAIC_MIG   = matrix(NA,num_n,1)
mean_difBIC_MIG   = matrix(NA,num_n,1)
sd_difAIC_MIG     = matrix(NA,num_n,1)
sd_difBIC_MIG     = matrix(NA,num_n,1)
perc_difAIC_MIG   = matrix(NA,num_n,1)
perc_difBIC_MIG   = matrix(NA,num_n,1)
perc_TR_MIG       = matrix(NA,num_n,1) 
#------------------------------------------
mean_test_statistic_MGA = matrix(NA,num_n,1)
mean_difAIC_MGA   = matrix(NA,num_n,1)
mean_difBIC_MGA   = matrix(NA,num_n,1)
sd_difAIC_MGA     = matrix(NA,num_n,1)
sd_difBIC_MGA     = matrix(NA,num_n,1)
perc_difAIC_MGA   = matrix(NA,num_n,1)
perc_difBIC_MGA   = matrix(NA,num_n,1)
perc_TR_MGA       = matrix(NA,num_n,1) 
#-------------------------------------------- 

# for each n 
for (w in 1:num_n){
  
  n = sample_size[w]
#-------------------------------------------
## Global variables
#-------------------------------------------
NREP    <- 1000  ## Monte Carlo trials
ALPHANS <- 0.05 # confidence level 
#-------------------------------------------
#parametros reais
beta.real   <- c(-1.4313,1.1763, 2.5124, 4.1351,-0.7945,-0.3415,0.0415, 0.0066)
alpha.real  <- -3.68
nu.real     <- 1.42
m.real      <- 0.66
lambda.real <- c(alpha.real,nu.real);lambda.real

max.censura  <- 20
r            <- length(beta.real)
psi.real     <- c(beta.real,alpha.real,nu.real,m.real);psi.real
psi          <- psi.real
#-------------------------------------------
# distribution for lifetime 
# dist=1 #Weibull
dist   = 1  
#-------------------------------------------
## Sample loop
AIC_MIG_mle <- matrix(NA,NREP,length(n)) 
AIC_MGA_mle <- matrix(NA,NREP,length(n))
AIC_POI_mle  <- matrix(NA,NREP,length(n))
#-------
BIC_MIG_mle <- matrix(NA,NREP,length(n)) 
BIC_MGA_mle <- matrix(NA,NREP,length(n))
BIC_POI_mle  <- matrix(NA,NREP,length(n)) 
#-------
maximo.logMIG <- matrix(NA,NREP,length(n))
maximo.logMGA <- matrix(NA,NREP,length(n))
maximo.logPOI <- matrix(NA,NREP,length(n))
#-------
difAIC_MIG_mle <- matrix(NA,NREP,length(n))
difBIC_MIG_mle <- matrix(NA,NREP,length(n))
difAIC_MGA_mle <- matrix(NA,NREP,length(n))
difBIC_MGA_mle <- matrix(NA,NREP,length(n))
#-------
#-------
test_statistic_MIG <- matrix(NA,NREP,length(n))
test_statistic_MGA<- matrix(NA,NREP,length(n))


count_mle_mig<- 0L
count_mle_mga<- 0L
#------------------
#-------------------------------------------
# MC Study
jj=1
while(jj<=NREP){

cat("replicate=", jj, "\n")
cat("sample size=", w, "\n")
#------------------------------------------- 
# genarating Data from Poisson distribution 
#-------------------------------------------  
  Data   <- data.generator_PTcr(n,beta.real,alpha.real,nu.real,max.censura)
  t      <- Data$t
  delta  <- Data$delta
  x      <- Data$x
#---------------------------------------------------------------------------------------------------------------------------
# testing MIG CASE
#---------------------------------------------------------------------------------------------------------------------------
  i     = 1
  dif   = 10
  lower = c(-Inf,0)
  if(dist==3 || dist==4){lower=c(0,0)}
  #-------------------------------------------
  # mixture model 
  # model=1: MIG; model=2: MGA
  model = 1
#-------------------------------------------
  psi         = psi.real
#------------------------------------------- 
# EM algorithm
#-------------------------------------------
  while(i<=10000 && dif>1e-4)
  {
    latentes  = passoE1(psi,Data,dist,model)
    Y         = latentes$Y
    Z         = latentes$Z
    N         = latentes$N
    B         = latentes$B
    beta      = psi[1:r]
    lambda    = psi[(r+1):(r+2)]
    m         = psi[r+3]
    maximo1   = optim(lambda,Q1,method = "L-BFGS-B",lower=lower,N=N,Data=Data,dist=dist,hessian=F)
    lambda    = maximo1$par
    maximo2   = optim(beta,Q2,method = "BFGS",N=N,Z=Z,Data=Data,hessian=F)
    beta      = maximo2$par  
    maximo3   = optim(m,Q3,method = "Brent",lower=0.0000001,upper=0.9999999,Y=Y ,B=B,model=model,hessian=F)
    m         = maximo3$par  
    psi.aux   = c(beta,lambda, m) 
    dif       = max(abs(psi-psi.aux))  
    psi       = psi.aux 
    i=i+1
  }

#-------------------------------------------
psi_mig      = matrix(psi,nrow=1);colnames(psi_mig) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");round(psi_mig,4)
  
#---------------------------------------------------------------------------------------------------------------------------
#  testing MGA CASE
#---------------------------------------------------------------------------------------------------------------------------
t      <- Data$t
delta  <- Data$delta
x      <- Data$x
#-------------------------------------------
# mixture model 
#  MGA
model = 2
#-------------------------------------------  
i     = 1
dif   = 10
lower = c(-Inf,0)
if(dist==3 || dist==4){lower=c(0,0)}
#-------------------------------------------
psi         = psi.real
#------------------------------------------- 
# EM algorithm
#-------------------------------------------
while(i<=10000 && dif>1e-4)
{
  latentes2  = passoE1(psi,Data,dist,model=2)
  Y         = latentes2$Y
  Z         = latentes2$Z
  N         = latentes2$N
  A         = latentes2$A
  beta      = psi[1:r]
  lambda    = psi[(r+1):(r+2)]
  m         = psi[r+3]
  maximo1   = optim(lambda,Q1,method = "L-BFGS-B",lower=lower,N=N,Data=Data,dist=dist,hessian=F)
  lambda    = maximo1$par  
  maximo2   = optim(beta,Q2,method = "BFGS",N=N,Z=Z,Data=Data,hessian=F)
  beta      = maximo2$par 
  maximo3   = optim(m,Q3,method = "Brent",lower=0.0000001,upper=0.9999999,Y=Y,A=A ,model=2,hessian=F)
  m         = maximo3$par 
  psi.aux   = c(beta,lambda, m) 
  dif       = max(abs(psi-psi.aux))  
  psi       = psi.aux #atualiza psi
  #----
  i=i+1
}# end EM

#-------------------------------------------
psi_mga      = matrix(psi,nrow=1);colnames(psi_mga) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");round(psi_mga,4)
 
#---------------------------------------------------------------------------------------------------------------------------
#  testing POISSON CASE
#---------------------------------------------------------------------------------------------------------------------------
# M~ Poisson (1); W~Weibul (2)   - Package(PScr)
#-------------------------------------------------------------
EM_PO  = EM.PScr(t, delta, x, model=1, dist=2,max.iter = 10000)
# estimated vector of parameters 
estimates_poi = EM_PO$estimate
psi_poi   =  matrix(estimates_poi[,1],1);colnames(psi_poi) = c(paste("beta",0:(r-1),sep=""),"alpha","sigma");round(psi_poi,4)
 
  
    #-------------------------
    #  AIC e BIC for MPIGcr
    maximo.logMIG[jj] = -llikeobserved(psi_mig,Data=Data,dist =dist,model=1) 
    AIC_MIG_mle[jj]       = -2*(maximo.logMIG[jj]) + 2*length(psi_mig) 
    BIC_MIG_mle[jj]       = -2*(maximo.logMIG[jj]) + log(length(t))*length(psi_mig) 
    #------------------------------------
    #  AIC e BIC for MBNcr
    maximo.logMGA[jj] = -llikeobserved(psi_mga,Data=Data,dist =dist,model=2) 
    AIC_MGA_mle[jj]       = -2*(maximo.logMGA[jj]) + 2*length(psi_mga) 
    BIC_MGA_mle[jj]       = -2*(maximo.logMGA[jj]) + log(length(t))*length(psi_mga)  
    #------------------------------------
    #AIC e BIC  for  Poisson
    maximo.logPOI[jj] <- -EM_PO$loglike; 
    AIC_POI_mle[jj] = EM_PO$AIC 
    BIC_POI_mle[jj] = EM_PO$BIC 
    
    #--- diferenças------------------------------------------
    # MIG x poisson
    difAIC_MIG_mle[jj]  <- AIC_MIG_mle[jj] - AIC_POI_mle[jj]
    difBIC_MIG_mle[jj]  <- BIC_MIG_mle[jj] - BIC_POI_mle[jj]
    # MGA x poisson
    difAIC_MGA_mle[jj]  <- AIC_MGA_mle[jj] - AIC_POI_mle[jj]
    difBIC_MGA_mle[jj]  <- BIC_MGA_mle[jj] - BIC_POI_mle[jj]
    # A: Modelo POISSON e  B: MIG ou MGA 
    test_statistic_MIG[jj] <-   -2 * ( maximo.logPOI[jj] - (  maximo.logMIG[jj]))
    test_statistic_MGA[jj] <-   -2 * ( maximo.logPOI[jj] - (  maximo.logMGA[jj]))
    ALPHA        <- 0.05
    
    # Comparando a estatística do teste com o valor crítico
    if (test_statistic_MIG[jj] > qchisq(1 - 2*ALPHA, df = 1)) {
      count_mle_mig <- count_mle_mig + 1L  # rejeição  
    }
    if (test_statistic_MGA[jj] > qchisq(1 - 2*ALPHA, df = 1)) {
      count_mle_mga <- count_mle_mga + 1L  # rejeição 
    } 
    
 
  #-------------------------
  jj<-jj+1  
  #-------------------------
  
}


##-------------------------------------------------------------------------------
## Numerical Results  - MIG
#-------------------------------------------------------------------------------
# diferences MIG MLE
#-------------------------------------------------------------------------------
mean_difAIC_MIG_mle  <- mean(difAIC_MIG_mle) ; mean_difAIC_MIG_mle  
mean_difBIC_MIG_mle  <- mean(difBIC_MIG_mle) ; mean_difBIC_MIG_mle  
sd_difAIC_MIG_mle    <- sd(difAIC_MIG_mle) ; sd_difAIC_MIG_mle 
sd_difBIC_MIG_mle    <- sd(difBIC_MIG_mle) ; sd_difBIC_MIG_mle 
perc_difAIC_MIG_mle  <- mean(difAIC_MIG_mle<0);perc_difAIC_MIG_mle 
perc_difBIC_MIG_mle  <- mean(difBIC_MIG_mle<0);perc_difBIC_MIG_mle 
perc_TR_MIG_mle      <- count_mle_mig/NREP ;perc_TR_MIG_mle 
mean_test_statistic_MIG_mle <- mean(test_statistic_MIG)


#------------------------------------------------------------------------------- 
# diferences MGA MLE
#-------------------------------------------------------------------------------
mean_difAIC_MGA_mle  <- mean(difAIC_MGA_mle) ; mean_difAIC_MGA_mle  
mean_difBIC_MGA_mle  <- mean(difBIC_MGA_mle) ; mean_difBIC_MGA_mle  
sd_difAIC_MGA_mle    <- sd(difAIC_MGA_mle) ; sd_difAIC_MGA_mle 
sd_difBIC_MGA_mle    <- sd(difBIC_MGA_mle) ; sd_difBIC_MGA_mle 
perc_difAIC_MGA_mle  <- mean(difAIC_MGA_mle<0) ; perc_difAIC_MGA_mle
perc_difBIC_MGA_mle  <- mean(difBIC_MGA_mle<0) ; perc_difBIC_MGA_mle
perc_TR_MGA_mle      <- count_mle_mga/NREP  ; perc_TR_MGA_mle 
mean_test_statistic_MGA <- mean(test_statistic_MGA)

#---------------------------------------------------
 
#---------------------------------------------------
# diferencas MIG
mean_difAIC_MIG[w] = mean_difAIC_MIG_mle  
mean_difBIC_MIG[w] = mean_difBIC_MIG_mle
sd_difAIC_MIG[w]   = sd_difAIC_MIG_mle      
sd_difBIC_MIG[w]   = sd_difBIC_MIG_mle  
perc_difAIC_MIG[w] =  perc_difAIC_MIG_mle 
perc_difBIC_MIG[w] = perc_difBIC_MIG_mle  
perc_TR_MIG[w] = perc_TR_MIG_mle
mean_test_statistic_MIG[w] = mean_test_statistic_MIG_mle  
#---------------------------------------------------
# diferencas MGA
mean_difAIC_MGA[w] = mean_difAIC_MGA_mle  
mean_difBIC_MGA[w] = mean_difBIC_MGA_mle
sd_difAIC_MGA[w]   = sd_difAIC_MGA_mle      
sd_difBIC_MGA[w]   = sd_difBIC_MGA_mle  
perc_difAIC_MGA[w] =  perc_difAIC_MGA_mle 
perc_difBIC_MGA[w] = perc_difBIC_MGA_mle  
perc_TR_MGA[w] = perc_TR_MIG_mle
mean_test_statistic_MGA[w] = mean_test_statistic_MIG_mle 


}


#---------------------------------------------------
# differences MPIGcr x PTcr
#---------------------------------------------------
resumo_mle_MIG <- data.frame(
  sample_size      = sample_size,    
  mean_difAIC_MIG  = round(mean_difAIC_MIG,4),
  mean_difBIC_MIG  = round(mean_difBIC_MIG,4),
  sd_difAIC_BIC_MIG    = round(sd_difAIC_MIG,4),
  perc_difAIC_MIG  = round(perc_difAIC_MIG,4),
  perc_difBIC_MIG  = round(perc_difBIC_MIG,4),
  perc_TR_MIG      = round(perc_TR_MIG,4)
) 

#---------------------------------------------------
# differences MNBcr x PTcr
#---------------------------------------------------
resumo_mle_MGA <- data.frame(
  sample_size      = sample_size,    
  mean_difAIC_MGA  = round(mean_difAIC_MGA,4),
  mean_difBIC_MGA  = round(mean_difBIC_MGA,4),
  sd_difAIC_BIC_MGA    = round(sd_difAIC_MGA,4),
  perc_difAIC_MGA  = round(perc_difAIC_MGA,4),
  perc_difBIC_MGA  = round(perc_difBIC_MGA,4),
  perc_TR_MGA      = round(perc_TR_MGA,4)
) 


resumo_mle_MIG2 <- resumo_mle_MIG
resumo_mle_MIG2$model <- "MPIGcr"
names(resumo_mle_MIG2) <- sub("_MIG$", "", names(resumo_mle_MIG2))

resumo_mle_MGA2 <- resumo_mle_MGA
resumo_mle_MGA2$model <- "MNBcr"
names(resumo_mle_MGA2) <- sub("_MGA$", "", names(resumo_mle_MGA2))

resumo_all <- rbind(resumo_mle_MGA2, resumo_mle_MIG2)
resumo_all <- resumo_all[, c("model","sample_size","mean_difAIC","mean_difBIC",
                             "sd_difAIC_BIC","perc_difAIC","perc_difBIC","perc_TR")]
resumo_all





#===================================================
# Simulation Study
# Table: Table 3: Differences in AIC and BIC between models fitted with MNBcr, MPIGcr and 
# PTcr competing causes, and  empirical rejection rate under H0 of the likelihood ratio test, 
# for data generated with PTcr competing causes under different sample sizes
#===================================================

name00 = paste("Table 3")
write.table(resumo_all ,file=paste(name00,".txt"))
write.csv(resumo_all ,file=paste(name00,".csv"))

save(list = ls(),file = "MC-Study-Misspecification-MIG-MGA-POI.RData")










 

  
  
 