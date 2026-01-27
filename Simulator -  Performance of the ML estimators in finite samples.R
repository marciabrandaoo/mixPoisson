
#==========================================================================================================================
# Simulation Study
#-----------------------------------------------
# Table 1 and 2: Empirical mean, median, standard deviation (SD), bias, standard error (SE), root mean square error (RMSE), and
#         coverage probability (CP) of the maximum likelihood estimators for the WEI distribution parameters in the concurrent
#           causes regression model under the MPIGcr and MNBcr scenarios.
#==========================================================================================================================
 
#------------------------------------------------------------------------
# setting seed used in the manuscript fixed configuration 
set.seed(104729)

#------------------------------------------------------------------------
# sample size
sample_size = c(200,300,500,1000)
num_n = length(sample_size)
 
media    = matrix(NA,num_n,11)
mediana  = matrix(NA,num_n,11)
SD       = matrix(NA,num_n,11)
bias     = matrix(NA,num_n,11)
SE       = matrix(NA,num_n,11)
rmse     = matrix(NA,num_n,11)
cover    = matrix(NA,num_n,11)

for (w in 1:num_n){
  
  n = sample_size[w]
  n=1000
#-------------------------------------------
## Global variables
#-------------------------------------------
NREP    <- 1000  
ALPHANS <- 0.05 
#-------------------------------------------
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
# mixture model 
# model=1: MIG; model=2: MGA
model = 1
#-------------------------------------------
## Sample loop
#-------------------------------------------
## MLE
beta0_mle <- matrix(NA,NREP,length(n))
beta1_mle <- matrix(NA,NREP,length(n))
beta2_mle <- matrix(NA,NREP,length(n))
beta3_mle <- matrix(NA,NREP,length(n))
beta4_mle <- matrix(NA,NREP,length(n))
beta5_mle <- matrix(NA,NREP,length(n))
beta6_mle <- matrix(NA,NREP,length(n))
beta7_mle <- matrix(NA,NREP,length(n))
alpha_mle <- matrix(NA,NREP,length(n))
nu_mle    <- matrix(NA,NREP,length(n))
m_mle     <- matrix(NA,NREP,length(n))
nu_mle2    <- matrix(NA,NREP,length(n))
m_mle2     <- matrix(NA,NREP,length(n))
#-------------------------------------------
# SE 
#-------------------------------------------
SE_beta0_mle  <- matrix(NA,NREP,length(n))
SE_beta1_mle  <- matrix(NA,NREP,length(n))
SE_beta2_mle  <- matrix(NA,NREP,length(n))
SE_beta3_mle  <- matrix(NA,NREP,length(n))
SE_beta4_mle  <- matrix(NA,NREP,length(n))
SE_beta5_mle  <- matrix(NA,NREP,length(n))
SE_beta6_mle  <- matrix(NA,NREP,length(n))
SE_beta7_mle  <- matrix(NA,NREP,length(n))
SE_alpha_mle  <- matrix(NA,NREP,length(n))
SE_nu_mle     <- matrix(NA,NREP,length(n))
SE_m_mle      <- matrix(NA,NREP,length(n))
SE_nu_mle2     <- matrix(NA,NREP,length(n))
SE_m_mle2      <- matrix(NA,NREP,length(n))
#-------------------------------------------
#IC
#-------------------------------------------
IC_beta0_mle  <- matrix(NA,NREP,2)  
IC_beta1_mle  <- matrix(NA,NREP,2)
IC_beta2_mle  <- matrix(NA,NREP,2)  
IC_beta3_mle  <- matrix(NA,NREP,2)  
IC_beta4_mle  <- matrix(NA,NREP,2)  
IC_beta5_mle  <- matrix(NA,NREP,2)  
IC_beta6_mle  <- matrix(NA,NREP,2)  
IC_beta7_mle  <- matrix(NA,NREP,2)  
IC_alpha_mle  <- matrix(NA,NREP,2)  
IC_nu_mle     <- matrix(NA,NREP,2)
IC_m_mle      <- matrix(NA,NREP,2)
#-------------------------------------------
# CP
#-------------------------------------------
coverage_beta0_mle  <- matrix(NA,NREP,length(n))
coverage_beta1_mle  <- matrix(NA,NREP,length(n))
coverage_beta2_mle  <- matrix(NA,NREP,length(n))
coverage_beta3_mle  <- matrix(NA,NREP,length(n))
coverage_beta4_mle  <- matrix(NA,NREP,length(n))
coverage_beta5_mle  <- matrix(NA,NREP,length(n))
coverage_beta6_mle  <- matrix(NA,NREP,length(n))
coverage_beta7_mle  <- matrix(NA,NREP,length(n))
coverage_alpha_mle  <- matrix(NA,NREP,length(n))
coverage_nu_mle     <- matrix(NA,NREP,length(n))
coverage_m_mle      <- matrix(NA,NREP,length(n))
#-------------------------------------------
## Bias
#-------------------------------------------
bias_beta0_mle  <- matrix(NA,NREP,length(n))
bias_beta1_mle  <- matrix(NA,NREP,length(n))
bias_beta2_mle  <- matrix(NA,NREP,length(n))
bias_beta3_mle  <- matrix(NA,NREP,length(n))
bias_beta4_mle  <- matrix(NA,NREP,length(n))
bias_beta5_mle  <- matrix(NA,NREP,length(n))
bias_beta6_mle  <- matrix(NA,NREP,length(n))
bias_beta7_mle  <- matrix(NA,NREP,length(n))
bias_alpha_mle  <- matrix(NA,NREP,length(n))
bias_nu_mle     <- matrix(NA,NREP,length(n))
bias_m_mle      <- matrix(NA,NREP,length(n))
#-------------------------------------------
# MC Study
jj=1
while(jj<=NREP){
#------------------------------------------- 
#generating data
  Data   <- data.generator(n,beta.real,alpha.real,nu.real,m.real,model,max.censura)
  t      <- Data$t
  delta  <- Data$delta
  x      <- Data$x
  
#-------------------------------------------  
  i     = 1
  dif   = 10
  lower = c(-Inf,0)
  if(dist==3 || dist==4){lower=c(0,0)}
#-------------------------------------------
  psi         = psi.real
#-------------------------------------------
#  set.seed(3072022)
#------------------------------------------- 
# EM algorithm
#-------------------------------------------
  while(i<=10000 && dif>1e-4)
  {
    latentes  = passoE(psi,Data,dist,model)
    Y         = latentes$Y
    Z         = latentes$Z
    N         = latentes$N
    A         = latentes$A
    B         = latentes$B
    beta      = psi[1:r]
    lambda    = psi[(r+1):(r+2)]
    m         = psi[r+3]
    maximo1   = optim(lambda,Q1,method = "L-BFGS-B",lower=lower,N=N,Data=Data,dist=dist,hessian=F)
    lambda    = maximo1$par  
    maximo2   = optim(beta,Q2,method = "BFGS",N=N,Z=Z,Data=Data,hessian=F)
    beta      = maximo2$par 
    maximo3   = optim(m,Q3,method = "Brent",lower=0.0000001,upper=0.9999999,Y=Y,A=A,B=B,model=model,hessian=F)
    m         = maximo3$par 
    psi.aux   = c(beta,lambda, m) 
    dif       = max(abs(psi-psi.aux))  
    psi       = psi.aux 
    i=i+1
    cat("iteration=", i, "\n")
    cat("replicate=", jj, "\n")
  } 

#-------------------------------------------
psi      = matrix(psi,nrow=1);colnames(psi) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");round(psi,7)
psi.real = matrix(psi.real,nrow=1);colnames(psi.real) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");psi.real

psi.log = c(psi[1:9], log(psi[10]),plogis(psi[11]))

  #----------------------------------------------------
  #  hessian PRACMA
  hessian.psi    =  solve(hessian(llikeobserved,x0=psi,Data=Data,dist=1,model=model)) 
  se  = sqrt(diag(matrix(hessian.psi,ncol=r+3,nrow=r+3)));
  
   if(all(!is.nan(se)) &&   all(!is.na(se)) )
 {
  colnames( hessian.psi) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m")
  rownames( hessian.psi) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m")
  
  estimates              = matrix(c(psi,sqrt(diag( hessian.psi))),ncol=2)
  colnames(estimates)    = c("Estimate","SE")
  rownames(estimates)    = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m")
  
  z_phi1_95      <- qnorm(ALPHANS/2)
  z_phi2_95      <- qnorm(1-ALPHANS/2)
  
  beta0_mle[jj,]   <- estimates[1,1] 
  beta1_mle[jj,]   <- estimates[2,1]
  beta2_mle[jj,]   <- estimates[3,1]
  beta3_mle[jj,]   <- estimates[4,1]
  beta4_mle[jj,]   <- estimates[5,1]
  beta5_mle[jj,]   <- estimates[6,1]
  beta6_mle[jj,]   <- estimates[7,1]
  beta7_mle[jj,]   <- estimates[8,1]
  alpha_mle[jj,]   <- estimates[9,1]
  nu_mle[jj,]      <- estimates[10,1]
  m_mle[jj,]       <- estimates[11,1]
  
  SE_beta0_mle[jj,] <- estimates[1,2]
  SE_beta1_mle[jj,] <- estimates[2,2]
  SE_beta2_mle[jj,] <- estimates[3,2]
  SE_beta3_mle[jj,] <- estimates[4,2]
  SE_beta4_mle[jj,] <- estimates[5,2]
  SE_beta5_mle[jj,] <- estimates[6,2]
  SE_beta6_mle[jj,] <- estimates[7,2]
  SE_beta7_mle[jj,] <- estimates[8,2]
  SE_alpha_mle[jj,] <- estimates[9,2]
  SE_nu_mle[jj,]    <- estimates[10,2]
  SE_m_mle[jj,]     <- estimates[11,2]

  bias_beta0_mle[jj,]  <- beta0_mle[jj,] - beta.real[1]
  bias_beta1_mle[jj,]  <- beta1_mle[jj,] - beta.real[2]
  bias_beta2_mle[jj,]  <- beta2_mle[jj,] - beta.real[3] 
  bias_beta3_mle[jj,]  <- beta3_mle[jj,] - beta.real[4] 
  bias_beta4_mle[jj,]  <- beta4_mle[jj,] - beta.real[5] 
  bias_beta5_mle[jj,]  <- beta5_mle[jj,] - beta.real[6] 
  bias_beta6_mle[jj,]  <- beta6_mle[jj,] - beta.real[7] 
  bias_beta7_mle[jj,]  <- beta7_mle[jj,] - beta.real[8] 
  bias_alpha_mle[jj,]  <- alpha_mle[jj,] - alpha.real   
  bias_nu_mle[jj,]     <- nu_mle[jj,]    - nu.real      
  bias_m_mle [jj,]     <- m_mle[jj,]     -  m.real    

  IC_beta0_mle[jj,] <- round(c(beta0_mle[jj,] + z_phi1_95 * SE_beta0_mle[jj,], beta0_mle[jj,] + z_phi2_95 * SE_beta0_mle[jj,]), 3)
  IC_beta1_mle[jj,] <- round(c(beta1_mle[jj,] + z_phi1_95 * SE_beta1_mle[jj,], beta1_mle[jj,] + z_phi2_95 * SE_beta1_mle[jj,]), 3)
  IC_beta2_mle[jj,] <- round(c(beta2_mle[jj,] + z_phi1_95 * SE_beta2_mle[jj,], beta2_mle[jj,] + z_phi2_95 * SE_beta2_mle[jj,]), 3)
  IC_beta3_mle[jj,] <- round(c(beta3_mle[jj,] + z_phi1_95 * SE_beta3_mle[jj,], beta3_mle[jj,] + z_phi2_95 * SE_beta3_mle[jj,]), 3)
  IC_beta4_mle[jj,] <- round(c(beta4_mle[jj,] + z_phi1_95 * SE_beta4_mle[jj,], beta4_mle[jj,] + z_phi2_95 * SE_beta4_mle[jj,]), 3)
  IC_beta5_mle[jj,] <- round(c(beta5_mle[jj,] + z_phi1_95 * SE_beta5_mle[jj,], beta5_mle[jj,] + z_phi2_95 * SE_beta5_mle[jj,]), 3)
  IC_beta6_mle[jj,] <- round(c(beta6_mle[jj,] + z_phi1_95 * SE_beta6_mle[jj,], beta6_mle[jj,] + z_phi2_95 * SE_beta6_mle[jj,]), 3)
  IC_beta7_mle[jj,] <- round(c(beta7_mle[jj,] + z_phi1_95 * SE_beta7_mle[jj,], beta7_mle[jj,] + z_phi2_95 * SE_beta7_mle[jj,]), 3)
  IC_alpha_mle[jj,] <- round(c(alpha_mle[jj,] + z_phi1_95 * SE_alpha_mle[jj,], alpha_mle[jj,] + z_phi2_95 * SE_alpha_mle[jj,]), 3)
  IC_nu_mle[jj,]    <- round(c(nu_mle[jj,]    + z_phi1_95 * SE_nu_mle[jj,]   , nu_mle[jj,]    + z_phi2_95 * SE_nu_mle[jj,]), 3)
  IC_m_mle[jj,]     <- round(c(m_mle[jj,]     + z_phi1_95 * SE_m_mle[jj,]    ,  m_mle[jj,]    + z_phi2_95 * SE_m_mle[jj,]), 3)

  coverage_beta0_mle[jj,] <- ifelse(IC_beta0_mle[jj,1] < beta.real[1] &  beta.real[1] <  IC_beta0_mle[jj,2], 1, 0)
  coverage_beta1_mle[jj,] <- ifelse(IC_beta1_mle[jj,1] < beta.real[2] &  beta.real[2] <  IC_beta1_mle[jj,2], 1, 0)
  coverage_beta2_mle[jj,] <- ifelse(IC_beta2_mle[jj,1] < beta.real[3] &  beta.real[3] <  IC_beta2_mle[jj,2], 1, 0)
  coverage_beta3_mle[jj,] <- ifelse(IC_beta3_mle[jj,1] < beta.real[4] &  beta.real[4] <  IC_beta3_mle[jj,2], 1, 0)
  coverage_beta4_mle[jj,] <- ifelse(IC_beta4_mle[jj,1] < beta.real[5] &  beta.real[5] <  IC_beta4_mle[jj,2], 1, 0)
  coverage_beta5_mle[jj,] <- ifelse(IC_beta5_mle[jj,1] < beta.real[6] &  beta.real[6] <  IC_beta5_mle[jj,2], 1, 0)
  coverage_beta6_mle[jj,] <- ifelse(IC_beta6_mle[jj,1] < beta.real[7] &  beta.real[7] <  IC_beta6_mle[jj,2], 1, 0)
  coverage_beta7_mle[jj,] <- ifelse(IC_beta7_mle[jj,1] < beta.real[8] &  beta.real[8] <  IC_beta7_mle[jj,2], 1, 0)
  coverage_alpha_mle[jj,] <- ifelse(IC_alpha_mle[jj,1] < alpha.real   &  alpha.real   <  IC_alpha_mle[jj,2],1, 0)
  coverage_nu_mle[jj,]    <- ifelse(IC_nu_mle[jj,1]    < nu.real      &  nu.real      <  IC_nu_mle[jj,2],1, 0)
  coverage_m_mle[jj,]     <- ifelse(IC_m_mle[jj,1]     < m.real       &  m.real       <  IC_m_mle[jj,2],1, 0)

   jj<-jj+1
  
  }
  
  
}


##-------------------------------------------------------------------------------
## Numerical Results
#------------------------------------------------------------------------------- 

#-------------------------------------------------------------------------------
#mean MLE
#-------------------------------------------------------------------------------
mean_beta0_mle  <- mean(beta0_mle) ; mean_beta0_mle
mean_beta1_mle  <- mean(beta1_mle) ; mean_beta1_mle
mean_beta2_mle  <- mean(beta2_mle) ; mean_beta2_mle
mean_beta3_mle  <- mean(beta3_mle) ; mean_beta3_mle
mean_beta4_mle  <- mean(beta4_mle) ; mean_beta4_mle
mean_beta5_mle  <- mean(beta5_mle) ; mean_beta5_mle
mean_beta6_mle  <- mean(beta6_mle) ; mean_beta6_mle
mean_beta7_mle  <- mean(beta7_mle) ; mean_beta7_mle
mean_alpha_mle  <- mean(alpha_mle) ; mean_alpha_mle
mean_nu_mle     <- mean(nu_mle)    ; mean_nu_mle
mean_m_mle      <- mean(m_mle)     ; mean_m_mle
#-------------------------------------------------------------------------------
#median MLE
median_beta0_mle  <- median(beta0_mle) ; median_beta0_mle
median_beta1_mle  <- median(beta1_mle) ; median_beta1_mle
median_beta2_mle  <- median(beta2_mle) ; median_beta2_mle
median_beta3_mle  <- median(beta3_mle) ; median_beta3_mle
median_beta4_mle  <- median(beta4_mle) ; median_beta4_mle
median_beta5_mle  <- median(beta5_mle) ; median_beta5_mle
median_beta6_mle  <- median(beta6_mle) ; median_beta6_mle
median_beta7_mle  <- median(beta7_mle) ; median_beta7_mle
median_alpha_mle  <- median(alpha_mle) ; median_alpha_mle
median_nu_mle     <- median(nu_mle)    ; median_nu_mle
median_m_mle      <- median(m_mle)     ; median_m_mle
#-------------------------------------------------------------------------------
## Variance
sd_beta0_mle  <- sd(beta0_mle); sd_beta0_mle
sd_beta1_mle  <- sd(beta1_mle); sd_beta1_mle
sd_beta2_mle  <- sd(beta2_mle); sd_beta2_mle
sd_beta3_mle  <- sd(beta3_mle); sd_beta3_mle
sd_beta4_mle  <- sd(beta4_mle); sd_beta4_mle
sd_beta5_mle  <- sd(beta5_mle); sd_beta5_mle
sd_beta6_mle  <- sd(beta6_mle); sd_beta6_mle
sd_beta7_mle  <- sd(beta7_mle); sd_beta7_mle
sd_alpha_mle  <- sd(alpha_mle); sd_alpha_mle
sd_nu_mle     <- sd(nu_mle)   ; sd_nu_mle
sd_m_mle      <- sd(m_mle)    ; sd_m_mle
#-------------------------------------------------------------------------------
## Bias
bias_beta0_mle  <- mean(beta0_mle) - beta.real[1]; bias_beta0_mle
bias_beta1_mle  <- mean(beta1_mle) - beta.real[2]; bias_beta1_mle
bias_beta2_mle  <- mean(beta2_mle) - beta.real[3]; bias_beta2_mle
bias_beta3_mle  <- mean(beta3_mle) - beta.real[4]; bias_beta3_mle
bias_beta4_mle  <- mean(beta4_mle) - beta.real[5]; bias_beta4_mle
bias_beta5_mle  <- mean(beta5_mle) - beta.real[6]; bias_beta5_mle
bias_beta6_mle  <- mean(beta6_mle) - beta.real[7]; bias_beta6_mle
bias_beta7_mle  <- mean(beta7_mle) - beta.real[8]; bias_beta7_mle
bias_alpha_mle  <- mean(alpha_mle) - alpha.real  ; bias_alpha_mle
bias_nu_mle     <- mean(nu_mle)    - nu.real     ; bias_nu_mle
bias_m_mle      <- mean(m_mle)    -  m.real      ; bias_m_mle
#-------------------------------------------------------------------------------
# mean SE
mean_SE_beta0_mle <-mean(SE_beta0_mle) ;mean_SE_beta0_mle 
mean_SE_beta1_mle <-mean(SE_beta1_mle) ;mean_SE_beta1_mle 
mean_SE_beta2_mle <-mean(SE_beta2_mle) ;mean_SE_beta2_mle
mean_SE_beta3_mle <-mean(SE_beta3_mle) ;mean_SE_beta3_mle
mean_SE_beta4_mle <-mean(SE_beta4_mle) ;mean_SE_beta4_mle
mean_SE_beta5_mle <-mean(SE_beta5_mle) ;mean_SE_beta5_mle
mean_SE_beta6_mle <-mean(SE_beta6_mle) ;mean_SE_beta6_mle
mean_SE_beta7_mle <-mean(SE_beta7_mle) ;mean_SE_beta7_mle
mean_SE_alpha_mle <-mean(SE_alpha_mle) ;mean_SE_alpha_mle
mean_SE_nu_mle    <-mean(SE_nu_mle)    ;mean_SE_nu_mle
mean_SE_m_mle     <-mean(SE_m_mle)     ;mean_SE_m_mle
#-------------------------------------------------------------------------------
## root of MSE
rmse_beta0_mle  <- sqrt(var(beta0_mle) + bias_beta0_mle^2); rmse_beta0_mle
rmse_beta1_mle  <- sqrt(var(beta1_mle) + bias_beta1_mle^2); rmse_beta1_mle
rmse_beta2_mle  <- sqrt(var(beta2_mle) + bias_beta2_mle^2); rmse_beta2_mle
rmse_beta3_mle  <- sqrt(var(beta3_mle) + bias_beta3_mle^2); rmse_beta3_mle
rmse_beta4_mle  <- sqrt(var(beta4_mle) + bias_beta4_mle^2); rmse_beta4_mle
rmse_beta5_mle  <- sqrt(var(beta5_mle) + bias_beta5_mle^2); rmse_beta5_mle
rmse_beta6_mle  <- sqrt(var(beta6_mle) + bias_beta6_mle^2); rmse_beta6_mle
rmse_beta7_mle  <- sqrt(var(beta7_mle) + bias_beta7_mle^2); rmse_beta7_mle
rmse_alpha_mle  <- sqrt(var(alpha_mle) + bias_alpha_mle^2); rmse_alpha_mle
rmse_nu_mle     <- sqrt(var(nu_mle)    + bias_nu_mle^2)   ; rmse_nu_mle
rmse_m_mle      <- sqrt(var(m_mle)     + bias_m_mle^2)    ; rmse_m_mle
#-------------------------------------------------------------------------------
#coverage of IC
mean_coverage_beta0_mle <- mean(coverage_beta0_mle); mean_coverage_beta0_mle
mean_coverage_beta1_mle <- mean(coverage_beta1_mle); mean_coverage_beta1_mle 
mean_coverage_beta2_mle <- mean(coverage_beta2_mle); mean_coverage_beta2_mle 
mean_coverage_beta3_mle <- mean(coverage_beta3_mle); mean_coverage_beta3_mle 
mean_coverage_beta4_mle <- mean(coverage_beta4_mle); mean_coverage_beta4_mle 
mean_coverage_beta5_mle <- mean(coverage_beta5_mle); mean_coverage_beta5_mle 
mean_coverage_beta6_mle <- mean(coverage_beta6_mle); mean_coverage_beta6_mle 
mean_coverage_beta7_mle <- mean(coverage_beta7_mle); mean_coverage_beta7_mle 
mean_coverage_alpha_mle <- mean(coverage_alpha_mle); mean_coverage_alpha_mle 
mean_coverage_nu_mle    <- mean(coverage_nu_mle)   ; mean_coverage_nu_mle
mean_coverage_m_mle     <- mean(coverage_m_mle)    ; mean_coverage_m_mle

#---------------------------------------------------


media[w,]   = c(mean_beta0_mle, mean_beta1_mle, mean_beta2_mle,mean_beta3_mle,mean_beta4_mle,mean_beta5_mle,mean_beta6_mle,mean_beta7_mle, mean_alpha_mle, mean_nu_mle ,mean_m_mle)
mediana[w,] = c(median_beta0_mle, median_beta1_mle, median_beta2_mle,median_beta3_mle,median_beta4_mle,median_beta5_mle,median_beta6_mle,median_beta7_mle, median_alpha_mle, median_nu_mle ,median_m_mle)
SD[w,]      = c(sd_beta0_mle, sd_beta1_mle, sd_beta2_mle, sd_beta3_mle,sd_beta4_mle,sd_beta5_mle,sd_beta6_mle,sd_beta7_mle,sd_alpha_mle, sd_nu_mle  ,sd_m_mle )
bias[w,]   = c(bias_beta0_mle, bias_beta1_mle, bias_beta2_mle,bias_beta3_mle,bias_beta4_mle,bias_beta5_mle,bias_beta6_mle,bias_beta7_mle,bias_alpha_mle, bias_nu_mle ,bias_m_mle)
SE[w,]     = c(mean_SE_beta0_mle, mean_SE_beta1_mle, mean_SE_beta2_mle, mean_SE_beta3_mle,mean_SE_beta4_mle,mean_SE_beta5_mle,mean_SE_beta6_mle,mean_SE_beta7_mle,mean_SE_alpha_mle, mean_SE_nu_mle  ,mean_SE_m_mle )
rmse[w,]    = c(rmse_beta0_mle, rmse_beta1_mle, rmse_beta2_mle,  rmse_beta3_mle,  rmse_beta4_mle,  rmse_beta5_mle,  rmse_beta6_mle,  rmse_beta7_mle, rmse_alpha_mle,rmse_nu_mle, rmse_m_mle )
cover[w,]   = c(mean_coverage_beta0_mle,mean_coverage_beta1_mle,mean_coverage_beta2_mle,mean_coverage_beta3_mle,mean_coverage_beta4_mle,mean_coverage_beta5_mle,mean_coverage_beta6_mle,mean_coverage_beta7_mle,mean_coverage_alpha_mle,mean_coverage_nu_mle,mean_coverage_m_mle) 
#---------

}



n200  = matrix( rbind(media[1,],mediana[1,],SD[1,],bias[1,] ,SE[1,], rmse[1,], cover[1,]),nrow=7)
n300  = matrix( rbind(media[2,],mediana[2,],SD[2,],bias[2,] ,SE[2,], rmse[2,], cover[2,]),nrow=7)
n500  = matrix( rbind(media[3,],mediana[3,],SD[3,],bias[3,] ,SE[3,], rmse[3,], cover[3,]),nrow=7)
n1000 = matrix( rbind(media[4,],mediana[4,],SD[4,],bias[4,] ,SE[4,], rmse[4,], cover[4,]),nrow=7)
ns = c(rep(200,7), rep(300,7),rep(500,7),rep(1000,7))





#==========================================================================================================================
# Simulation Study
#-----------------------------------------------
# Table 1: Empirical mean, median, standard deviation (SD), bias, standard error (SE), root mean square error (RMSE), and
#         coverage probability (CP) of the maximum likelihood estimators for the WEI distribution parameters in the concurrent
#           causes regression model under the MPIGcr scenario.
#==========================================================================================================================

sim_MC_results1 = round(cbind(ns,matrix(rbind(n200,n300,n500, n1000),nrow=28)),4)
dimnames(sim_MC_results1) = list( rep(c("Mean", "Median", "SD", "Bias","SE", "RMSE","CP"),4),c("Sample_Size","beta_0","beta_1", "beta_2", "beta_3","beta_4","beta_5","beta_6","beta_7","alpha","nu","m" ))


sim_MC_results1=rbind(c("n",c(psi.real)),sim_MC_results1)
sim_MC_results1 = data.frame(sim_MC_results1)

name1 = paste("Table_MC_results-MODEL_MPIGcr")
write.table(sim_MC_results1 ,file=paste(name1,".txt"))
write.csv(sim_MC_results1 ,file=paste(name1,".csv"))





################################################################################################################################

sample_size = c(200,300,500,1000)
num_n = length(sample_size)

media_2    = matrix(NA,num_n,11)
mediana_2  = matrix(NA,num_n,11)
SD_2       = matrix(NA,num_n,11)
bias_2     = matrix(NA,num_n,11)
SE_2       = matrix(NA,num_n,11)
rmse_2     = matrix(NA,num_n,11)
cover_2    = matrix(NA,num_n,11)

for (w in 1:num_n){
  
  n = sample_size[w]
  n=1000
  #-------------------------------------------
  ## Global variables
  #-------------------------------------------
  NREP    <- 1000  
  ALPHANS <- 0.05 
  #-------------------------------------------
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
  # mixture model 
  # model=1: MIG; model=2: MGA
  model = 2
  #-------------------------------------------
  ## Sample loop
  #-------------------------------------------
  ## MLE
  beta0_mle <- matrix(NA,NREP,length(n))
  beta1_mle <- matrix(NA,NREP,length(n))
  beta2_mle <- matrix(NA,NREP,length(n))
  beta3_mle <- matrix(NA,NREP,length(n))
  beta4_mle <- matrix(NA,NREP,length(n))
  beta5_mle <- matrix(NA,NREP,length(n))
  beta6_mle <- matrix(NA,NREP,length(n))
  beta7_mle <- matrix(NA,NREP,length(n))
  alpha_mle <- matrix(NA,NREP,length(n))
  nu_mle    <- matrix(NA,NREP,length(n))
  m_mle     <- matrix(NA,NREP,length(n))
  nu_mle2    <- matrix(NA,NREP,length(n))
  m_mle2     <- matrix(NA,NREP,length(n))
  #-------------------------------------------
  # SE 
  #-------------------------------------------
  SE_beta0_mle  <- matrix(NA,NREP,length(n))
  SE_beta1_mle  <- matrix(NA,NREP,length(n))
  SE_beta2_mle  <- matrix(NA,NREP,length(n))
  SE_beta3_mle  <- matrix(NA,NREP,length(n))
  SE_beta4_mle  <- matrix(NA,NREP,length(n))
  SE_beta5_mle  <- matrix(NA,NREP,length(n))
  SE_beta6_mle  <- matrix(NA,NREP,length(n))
  SE_beta7_mle  <- matrix(NA,NREP,length(n))
  SE_alpha_mle  <- matrix(NA,NREP,length(n))
  SE_nu_mle     <- matrix(NA,NREP,length(n))
  SE_m_mle      <- matrix(NA,NREP,length(n))
  SE_nu_mle2     <- matrix(NA,NREP,length(n))
  SE_m_mle2      <- matrix(NA,NREP,length(n))
  #-------------------------------------------
  #IC
  #-------------------------------------------
  IC_beta0_mle  <- matrix(NA,NREP,2)  
  IC_beta1_mle  <- matrix(NA,NREP,2)
  IC_beta2_mle  <- matrix(NA,NREP,2)  
  IC_beta3_mle  <- matrix(NA,NREP,2)  
  IC_beta4_mle  <- matrix(NA,NREP,2)  
  IC_beta5_mle  <- matrix(NA,NREP,2)  
  IC_beta6_mle  <- matrix(NA,NREP,2)  
  IC_beta7_mle  <- matrix(NA,NREP,2)  
  IC_alpha_mle  <- matrix(NA,NREP,2)  
  IC_nu_mle     <- matrix(NA,NREP,2)
  IC_m_mle      <- matrix(NA,NREP,2)
  #-------------------------------------------
  # CP
  #-------------------------------------------
  coverage_beta0_mle  <- matrix(NA,NREP,length(n))
  coverage_beta1_mle  <- matrix(NA,NREP,length(n))
  coverage_beta2_mle  <- matrix(NA,NREP,length(n))
  coverage_beta3_mle  <- matrix(NA,NREP,length(n))
  coverage_beta4_mle  <- matrix(NA,NREP,length(n))
  coverage_beta5_mle  <- matrix(NA,NREP,length(n))
  coverage_beta6_mle  <- matrix(NA,NREP,length(n))
  coverage_beta7_mle  <- matrix(NA,NREP,length(n))
  coverage_alpha_mle  <- matrix(NA,NREP,length(n))
  coverage_nu_mle     <- matrix(NA,NREP,length(n))
  coverage_m_mle      <- matrix(NA,NREP,length(n))
  #-------------------------------------------
  ## Bias
  #-------------------------------------------
  bias_beta0_mle  <- matrix(NA,NREP,length(n))
  bias_beta1_mle  <- matrix(NA,NREP,length(n))
  bias_beta2_mle  <- matrix(NA,NREP,length(n))
  bias_beta3_mle  <- matrix(NA,NREP,length(n))
  bias_beta4_mle  <- matrix(NA,NREP,length(n))
  bias_beta5_mle  <- matrix(NA,NREP,length(n))
  bias_beta6_mle  <- matrix(NA,NREP,length(n))
  bias_beta7_mle  <- matrix(NA,NREP,length(n))
  bias_alpha_mle  <- matrix(NA,NREP,length(n))
  bias_nu_mle     <- matrix(NA,NREP,length(n))
  bias_m_mle      <- matrix(NA,NREP,length(n))
  #-------------------------------------------
  # MC Study
  jj=1
  while(jj<=NREP){
    #------------------------------------------- 
    #gerando os dados
    Data   <- data.generator(n,beta.real,alpha.real,nu.real,m.real,model,max.censura)
    t      <- Data$t
    delta  <- Data$delta
    x      <- Data$x
 
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
      latentes  = passoE(psi,Data,dist,model)
      #-------------------------------------------
      Y         = latentes$Y
      Z         = latentes$Z
      N         = latentes$N
      A         = latentes$A
      B         = latentes$B
      beta      = psi[1:r]
      lambda    = psi[(r+1):(r+2)]
      m         = psi[r+3]
      maximo1   = optim(lambda,Q1,method = "L-BFGS-B",lower=lower,N=N,Data=Data,dist=dist,hessian=F)
      lambda    = maximo1$par 
      maximo2   = optim(beta,Q2,method = "BFGS",N=N,Z=Z,Data=Data,hessian=F)
      beta      = maximo2$par 
      maximo3   = optim(m,Q3,method = "Brent",lower=0.0000001,upper=0.9999999,Y=Y,A=A,B=B,model=model,hessian=F)
      m         = maximo3$par 
      psi.aux   = c(beta,lambda, m) 
      dif       = max(abs(psi-psi.aux))  
      psi       = psi.aux 
      cat("iteration=", i, "\n")
      cat("replicate=", jj, "\n")
    } 
    
    #-------------------------------------------
    psi      = matrix(psi,nrow=1);colnames(psi) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");round(psi,7)
    psi.real = matrix(psi.real,nrow=1);colnames(psi.real) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");psi.real
     
    #----------------------------------------------------
    #  hessian PRACMA
    hessian.psi    =  solve(hessian(llikeobserved,x0=psi,Data=Data,dist=1,model=model)) 
    se  = sqrt(diag(matrix(hessian.psi,ncol=r+3,nrow=r+3)));
      if(all(!is.nan(se)) &&   all(!is.na(se)) )
    {
      colnames( hessian.psi) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m")
      rownames( hessian.psi) = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m")
      
      estimates              = matrix(c(psi,sqrt(diag( hessian.psi))),ncol=2)
      colnames(estimates)    = c("Estimate","SE")
      rownames(estimates)    = c(paste("beta",0:(r-1),sep=""),"alpha","nu","m")
      
       
      #----------------------------------------------------
      z_phi1_95      <- qnorm(ALPHANS/2)
      z_phi2_95      <- qnorm(1-ALPHANS/2)
      
      beta0_mle[jj,]   <- estimates[1,1] 
      beta1_mle[jj,]   <- estimates[2,1]
      beta2_mle[jj,]   <- estimates[3,1]
      beta3_mle[jj,]   <- estimates[4,1]
      beta4_mle[jj,]   <- estimates[5,1]
      beta5_mle[jj,]   <- estimates[6,1]
      beta6_mle[jj,]   <- estimates[7,1]
      beta7_mle[jj,]   <- estimates[8,1]
      alpha_mle[jj,]   <- estimates[9,1]
      nu_mle[jj,]      <- estimates[10,1]
      m_mle[jj,]       <- estimates[11,1]

      SE_beta0_mle[jj,] <- estimates[1,2]
      SE_beta1_mle[jj,] <- estimates[2,2]
      SE_beta2_mle[jj,] <- estimates[3,2]
      SE_beta3_mle[jj,] <- estimates[4,2]
      SE_beta4_mle[jj,] <- estimates[5,2]
      SE_beta5_mle[jj,] <- estimates[6,2]
      SE_beta6_mle[jj,] <- estimates[7,2]
      SE_beta7_mle[jj,] <- estimates[8,2]
      SE_alpha_mle[jj,] <- estimates[9,2]
      SE_nu_mle[jj,]    <- estimates[10,2]
      SE_m_mle[jj,]     <- estimates[11,2]

      bias_beta0_mle[jj,]  <- beta0_mle[jj,] - beta.real[1]
      bias_beta1_mle[jj,]  <- beta1_mle[jj,] - beta.real[2]
      bias_beta2_mle[jj,]  <- beta2_mle[jj,] - beta.real[3] 
      bias_beta3_mle[jj,]  <- beta3_mle[jj,] - beta.real[4] 
      bias_beta4_mle[jj,]  <- beta4_mle[jj,] - beta.real[5] 
      bias_beta5_mle[jj,]  <- beta5_mle[jj,] - beta.real[6] 
      bias_beta6_mle[jj,]  <- beta6_mle[jj,] - beta.real[7] 
      bias_beta7_mle[jj,]  <- beta7_mle[jj,] - beta.real[8] 
      bias_alpha_mle[jj,]  <- alpha_mle[jj,] - alpha.real   
      bias_nu_mle[jj,]     <- nu_mle[jj,]    - nu.real      
      bias_m_mle [jj,]     <- m_mle[jj,]     -  m.real    

      IC_beta0_mle[jj,] <- round(c(beta0_mle[jj,] + z_phi1_95 * SE_beta0_mle[jj,], beta0_mle[jj,] + z_phi2_95 * SE_beta0_mle[jj,]), 3)
      IC_beta1_mle[jj,] <- round(c(beta1_mle[jj,] + z_phi1_95 * SE_beta1_mle[jj,], beta1_mle[jj,] + z_phi2_95 * SE_beta1_mle[jj,]), 3)
      IC_beta2_mle[jj,] <- round(c(beta2_mle[jj,] + z_phi1_95 * SE_beta2_mle[jj,], beta2_mle[jj,] + z_phi2_95 * SE_beta2_mle[jj,]), 3)
      IC_beta3_mle[jj,] <- round(c(beta3_mle[jj,] + z_phi1_95 * SE_beta3_mle[jj,], beta3_mle[jj,] + z_phi2_95 * SE_beta3_mle[jj,]), 3)
      IC_beta4_mle[jj,] <- round(c(beta4_mle[jj,] + z_phi1_95 * SE_beta4_mle[jj,], beta4_mle[jj,] + z_phi2_95 * SE_beta4_mle[jj,]), 3)
      IC_beta5_mle[jj,] <- round(c(beta5_mle[jj,] + z_phi1_95 * SE_beta5_mle[jj,], beta5_mle[jj,] + z_phi2_95 * SE_beta5_mle[jj,]), 3)
      IC_beta6_mle[jj,] <- round(c(beta6_mle[jj,] + z_phi1_95 * SE_beta6_mle[jj,], beta6_mle[jj,] + z_phi2_95 * SE_beta6_mle[jj,]), 3)
      IC_beta7_mle[jj,] <- round(c(beta7_mle[jj,] + z_phi1_95 * SE_beta7_mle[jj,], beta7_mle[jj,] + z_phi2_95 * SE_beta7_mle[jj,]), 3)
      IC_alpha_mle[jj,] <- round(c(alpha_mle[jj,] + z_phi1_95 * SE_alpha_mle[jj,], alpha_mle[jj,] + z_phi2_95 * SE_alpha_mle[jj,]), 3)
      IC_nu_mle[jj,]    <- round(c(nu_mle[jj,]    + z_phi1_95 * SE_nu_mle[jj,]   , nu_mle[jj,]    + z_phi2_95 * SE_nu_mle[jj,]), 3)
      IC_m_mle[jj,]     <- round(c(m_mle[jj,]     + z_phi1_95 * SE_m_mle[jj,]    ,  m_mle[jj,]    + z_phi2_95 * SE_m_mle[jj,]), 3)
      
      coverage_beta0_mle[jj,] <- ifelse(IC_beta0_mle[jj,1] < beta.real[1] &  beta.real[1] <  IC_beta0_mle[jj,2], 1, 0)
      coverage_beta1_mle[jj,] <- ifelse(IC_beta1_mle[jj,1] < beta.real[2] &  beta.real[2] <  IC_beta1_mle[jj,2], 1, 0)
      coverage_beta2_mle[jj,] <- ifelse(IC_beta2_mle[jj,1] < beta.real[3] &  beta.real[3] <  IC_beta2_mle[jj,2], 1, 0)
      coverage_beta3_mle[jj,] <- ifelse(IC_beta3_mle[jj,1] < beta.real[4] &  beta.real[4] <  IC_beta3_mle[jj,2], 1, 0)
      coverage_beta4_mle[jj,] <- ifelse(IC_beta4_mle[jj,1] < beta.real[5] &  beta.real[5] <  IC_beta4_mle[jj,2], 1, 0)
      coverage_beta5_mle[jj,] <- ifelse(IC_beta5_mle[jj,1] < beta.real[6] &  beta.real[6] <  IC_beta5_mle[jj,2], 1, 0)
      coverage_beta6_mle[jj,] <- ifelse(IC_beta6_mle[jj,1] < beta.real[7] &  beta.real[7] <  IC_beta6_mle[jj,2], 1, 0)
      coverage_beta7_mle[jj,] <- ifelse(IC_beta7_mle[jj,1] < beta.real[8] &  beta.real[8] <  IC_beta7_mle[jj,2], 1, 0)
      coverage_alpha_mle[jj,] <- ifelse(IC_alpha_mle[jj,1] < alpha.real   &  alpha.real   <  IC_alpha_mle[jj,2],1, 0)
      coverage_nu_mle[jj,]    <- ifelse(IC_nu_mle[jj,1]    < nu.real      &  nu.real      <  IC_nu_mle[jj,2],1, 0)
      coverage_m_mle[jj,]     <- ifelse(IC_m_mle[jj,1]     < m.real       &  m.real       <  IC_m_mle[jj,2],1, 0)
 
      jj<-jj+1
      
    }
    
    
  }
  
  
  ##-------------------------------------------------------------------------------
  ## Numerical Results
  #------------------------------------------------------------------------------- 
  
  #-------------------------------------------------------------------------------
  #mean MLE
  #-------------------------------------------------------------------------------
  mean_beta0_mle  <- mean(beta0_mle) ; mean_beta0_mle
  mean_beta1_mle  <- mean(beta1_mle) ; mean_beta1_mle
  mean_beta2_mle  <- mean(beta2_mle) ; mean_beta2_mle
  mean_beta3_mle  <- mean(beta3_mle) ; mean_beta3_mle
  mean_beta4_mle  <- mean(beta4_mle) ; mean_beta4_mle
  mean_beta5_mle  <- mean(beta5_mle) ; mean_beta5_mle
  mean_beta6_mle  <- mean(beta6_mle) ; mean_beta6_mle
  mean_beta7_mle  <- mean(beta7_mle) ; mean_beta7_mle
  mean_alpha_mle  <- mean(alpha_mle) ; mean_alpha_mle
  mean_nu_mle     <- mean(nu_mle)    ; mean_nu_mle
  mean_m_mle      <- mean(m_mle)     ; mean_m_mle
  #-------------------------------------------------------------------------------
  #median MLE
  median_beta0_mle  <- median(beta0_mle) ; median_beta0_mle
  median_beta1_mle  <- median(beta1_mle) ; median_beta1_mle
  median_beta2_mle  <- median(beta2_mle) ; median_beta2_mle
  median_beta3_mle  <- median(beta3_mle) ; median_beta3_mle
  median_beta4_mle  <- median(beta4_mle) ; median_beta4_mle
  median_beta5_mle  <- median(beta5_mle) ; median_beta5_mle
  median_beta6_mle  <- median(beta6_mle) ; median_beta6_mle
  median_beta7_mle  <- median(beta7_mle) ; median_beta7_mle
  median_alpha_mle  <- median(alpha_mle) ; median_alpha_mle
  median_nu_mle     <- median(nu_mle)    ; median_nu_mle
  median_m_mle      <- median(m_mle)     ; median_m_mle
  #-------------------------------------------------------------------------------
  ## Variance
  sd_beta0_mle  <- sd(beta0_mle); sd_beta0_mle
  sd_beta1_mle  <- sd(beta1_mle); sd_beta1_mle
  sd_beta2_mle  <- sd(beta2_mle); sd_beta2_mle
  sd_beta3_mle  <- sd(beta3_mle); sd_beta3_mle
  sd_beta4_mle  <- sd(beta4_mle); sd_beta4_mle
  sd_beta5_mle  <- sd(beta5_mle); sd_beta5_mle
  sd_beta6_mle  <- sd(beta6_mle); sd_beta6_mle
  sd_beta7_mle  <- sd(beta7_mle); sd_beta7_mle
  sd_alpha_mle  <- sd(alpha_mle); sd_alpha_mle
  sd_nu_mle     <- sd(nu_mle)   ; sd_nu_mle
  sd_m_mle      <- sd(m_mle)    ; sd_m_mle
  #-------------------------------------------------------------------------------
  ## Bias
  bias_beta0_mle  <- mean(beta0_mle) - beta.real[1]; bias_beta0_mle
  bias_beta1_mle  <- mean(beta1_mle) - beta.real[2]; bias_beta1_mle
  bias_beta2_mle  <- mean(beta2_mle) - beta.real[3]; bias_beta2_mle
  bias_beta3_mle  <- mean(beta3_mle) - beta.real[4]; bias_beta3_mle
  bias_beta4_mle  <- mean(beta4_mle) - beta.real[5]; bias_beta4_mle
  bias_beta5_mle  <- mean(beta5_mle) - beta.real[6]; bias_beta5_mle
  bias_beta6_mle  <- mean(beta6_mle) - beta.real[7]; bias_beta6_mle
  bias_beta7_mle  <- mean(beta7_mle) - beta.real[8]; bias_beta7_mle
  bias_alpha_mle  <- mean(alpha_mle) - alpha.real  ; bias_alpha_mle
  bias_nu_mle     <- mean(nu_mle)    - nu.real     ; bias_nu_mle
  bias_m_mle      <- mean(m_mle)    -  m.real      ; bias_m_mle
  #-------------------------------------------------------------------------------
  # mean SE
  mean_SE_beta0_mle <-mean(SE_beta0_mle) ;mean_SE_beta0_mle 
  mean_SE_beta1_mle <-mean(SE_beta1_mle) ;mean_SE_beta1_mle 
  mean_SE_beta2_mle <-mean(SE_beta2_mle) ;mean_SE_beta2_mle
  mean_SE_beta3_mle <-mean(SE_beta3_mle) ;mean_SE_beta3_mle
  mean_SE_beta4_mle <-mean(SE_beta4_mle) ;mean_SE_beta4_mle
  mean_SE_beta5_mle <-mean(SE_beta5_mle) ;mean_SE_beta5_mle
  mean_SE_beta6_mle <-mean(SE_beta6_mle) ;mean_SE_beta6_mle
  mean_SE_beta7_mle <-mean(SE_beta7_mle) ;mean_SE_beta7_mle
  mean_SE_alpha_mle <-mean(SE_alpha_mle) ;mean_SE_alpha_mle
  mean_SE_nu_mle    <-mean(SE_nu_mle)    ;mean_SE_nu_mle
  mean_SE_m_mle     <-mean(SE_m_mle)     ;mean_SE_m_mle
  #-------------------------------------------------------------------------------
  ## root of MSE
  rmse_beta0_mle  <- sqrt(var(beta0_mle) + bias_beta0_mle^2); rmse_beta0_mle
  rmse_beta1_mle  <- sqrt(var(beta1_mle) + bias_beta1_mle^2); rmse_beta1_mle
  rmse_beta2_mle  <- sqrt(var(beta2_mle) + bias_beta2_mle^2); rmse_beta2_mle
  rmse_beta3_mle  <- sqrt(var(beta3_mle) + bias_beta3_mle^2); rmse_beta3_mle
  rmse_beta4_mle  <- sqrt(var(beta4_mle) + bias_beta4_mle^2); rmse_beta4_mle
  rmse_beta5_mle  <- sqrt(var(beta5_mle) + bias_beta5_mle^2); rmse_beta5_mle
  rmse_beta6_mle  <- sqrt(var(beta6_mle) + bias_beta6_mle^2); rmse_beta6_mle
  rmse_beta7_mle  <- sqrt(var(beta7_mle) + bias_beta7_mle^2); rmse_beta7_mle
  rmse_alpha_mle  <- sqrt(var(alpha_mle) + bias_alpha_mle^2); rmse_alpha_mle
  rmse_nu_mle     <- sqrt(var(nu_mle)    + bias_nu_mle^2)   ; rmse_nu_mle
  rmse_m_mle      <- sqrt(var(m_mle)     + bias_m_mle^2)    ; rmse_m_mle
  #-------------------------------------------------------------------------------
  #coverage of IC
  mean_coverage_beta0_mle <- mean(coverage_beta0_mle); mean_coverage_beta0_mle
  mean_coverage_beta1_mle <- mean(coverage_beta1_mle); mean_coverage_beta1_mle 
  mean_coverage_beta2_mle <- mean(coverage_beta2_mle); mean_coverage_beta2_mle 
  mean_coverage_beta3_mle <- mean(coverage_beta3_mle); mean_coverage_beta3_mle 
  mean_coverage_beta4_mle <- mean(coverage_beta4_mle); mean_coverage_beta4_mle 
  mean_coverage_beta5_mle <- mean(coverage_beta5_mle); mean_coverage_beta5_mle 
  mean_coverage_beta6_mle <- mean(coverage_beta6_mle); mean_coverage_beta6_mle 
  mean_coverage_beta7_mle <- mean(coverage_beta7_mle); mean_coverage_beta7_mle 
  mean_coverage_alpha_mle <- mean(coverage_alpha_mle); mean_coverage_alpha_mle 
  mean_coverage_nu_mle    <- mean(coverage_nu_mle)   ; mean_coverage_nu_mle
  mean_coverage_m_mle     <- mean(coverage_m_mle)    ; mean_coverage_m_mle
  
  #---------------------------------------------------
  
  
  media_2[w,]   = c(mean_beta0_mle, mean_beta1_mle, mean_beta2_mle,mean_beta3_mle,mean_beta4_mle,mean_beta5_mle,mean_beta6_mle,mean_beta7_mle, mean_alpha_mle, mean_nu_mle ,mean_m_mle)
  mediana_2[w,] = c(median_beta0_mle, median_beta1_mle, median_beta2_mle,median_beta3_mle,median_beta4_mle,median_beta5_mle,median_beta6_mle,median_beta7_mle, median_alpha_mle, median_nu_mle ,median_m_mle)
  SD_2[w,]      = c(sd_beta0_mle, sd_beta1_mle, sd_beta2_mle, sd_beta3_mle,sd_beta4_mle,sd_beta5_mle,sd_beta6_mle,sd_beta7_mle,sd_alpha_mle, sd_nu_mle  ,sd_m_mle )
  bias_2[w,]   = c(bias_beta0_mle, bias_beta1_mle, bias_beta2_mle,bias_beta3_mle,bias_beta4_mle,bias_beta5_mle,bias_beta6_mle,bias_beta7_mle,bias_alpha_mle, bias_nu_mle ,bias_m_mle)
  SE_2[w,]     = c(mean_SE_beta0_mle, mean_SE_beta1_mle, mean_SE_beta2_mle, mean_SE_beta3_mle,mean_SE_beta4_mle,mean_SE_beta5_mle,mean_SE_beta6_mle,mean_SE_beta7_mle,mean_SE_alpha_mle, mean_SE_nu_mle  ,mean_SE_m_mle )
  rmse_2[w,]    = c(rmse_beta0_mle, rmse_beta1_mle, rmse_beta2_mle,  rmse_beta3_mle,  rmse_beta4_mle,  rmse_beta5_mle,  rmse_beta6_mle,  rmse_beta7_mle, rmse_alpha_mle,rmse_nu_mle, rmse_m_mle )
  cover_2[w,]   = c(mean_coverage_beta0_mle,mean_coverage_beta1_mle,mean_coverage_beta2_mle,mean_coverage_beta3_mle,mean_coverage_beta4_mle,mean_coverage_beta5_mle,mean_coverage_beta6_mle,mean_coverage_beta7_mle,mean_coverage_alpha_mle,mean_coverage_nu_mle,mean_coverage_m_mle) 
  #---------
  
}



n200_2  = matrix( rbind(media_2[1,],mediana_2[1,],SD_2[1,],bias_2[1,] ,SE_2[1,], rmse_2[1,], cover_2[1,]),nrow=7)
n300_2  = matrix( rbind(media_2[2,],mediana_2[2,],SD_2[2,],bias_2[2,] ,SE_2[2,], rmse_2[2,], cover_2[2,]),nrow=7)
n500_2  = matrix( rbind(media_2[3,],mediana_2[3,],SD_2[3,],bias_2[3,] ,SE_2[3,], rmse_2[3,], cover_2[3,]),nrow=7)
n1000_2 = matrix( rbind(media_2[4,],mediana_2[4,],SD_2[4,],bias_2[4,] ,SE_2[4,], rmse_2[4,], cover_2[4,]),nrow=7)
ns = c(rep(200,7), rep(300,7),rep(500,7),rep(1000,7))



#==========================================================================================================================
# Simulation Study
#-----------------------------------------------
# Table 1: Empirical mean, median, standard deviation (SD), bias, standard error (SE), root mean square error (RMSE), and
#          coverage probability (CP) of the maximum likelihood estimators for the WEI distribution parameters in the concurrent
#           causes regression model under the MNBcr scenario.
#==========================================================================================================================

sim_MC_results2 = round(cbind(ns,matrix(rbind(n200_2,n300_2,n500_2, n1000_2),nrow=28)),4)
dimnames(sim_MC_results2) = list( rep(c("Mean", "Median", "SD", "Bias","SE", "RMSE","CP"),4),c("Sample_Size","beta_0","beta_1", "beta_2", "beta_3","beta_4","beta_5","beta_6","beta_7","alpha","nu","m" ))

sim_MC_results2=rbind(c("n",c(psi.real)),sim_MC_results2)
sim_MC_results2 = data.frame(sim_MC_results2)

name2 = paste("Table1")
write.table(sim_MC_results2 ,file=paste(name2,".txt"))
write.csv(sim_MC_results2 ,file=paste(name2,".csv"))




#------------------------------------------------------------------
save(list = ls(),file = "Table2_MC-Study-MIG MGA.RData")


  
  
 