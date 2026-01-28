#============================#
# Distributions for W (dist) #
# 1: Weibull                 #
# 2: Log-Normal              #
# 3: Gamma                   #      
#============================#
# Mixture Model              #
# 1: MPIGcr case             #
# 2: MNBcr case              #
#============================#
require(Bessel)
require(pracma)
require(VGAM)
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
    S    <-  pbisa(t, shape = alpha, scale = nu,  lower.tail = FALSE)
  } else {
    stop("Distribuição invalida. Escolha 1 (Weibull), 2 (Log-normal), 3 (Gamma) ou 4 (BS).")
  }
  #-----------------------------------------------------------------------------
  A<- B<- G <- H <- NULL
  #-MPIGcr case   (mix=1) -----------------------------------------------------------
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
  #-MNBcr case   (mix=2)------------------------------------------------------------
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
  } else if (dist == 4) { # Birnbaum-Saunders
       S    <-  pbisa(t, shape = alpha, scale = 1/nu,lower.tail = FALSE)
      logf <- dbisa(t, shape = alpha, scale = 1/nu, log = TRUE)
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
#-MPIGcr case---------------------------------------------------------------------- 
 if (model == 1) {
    q3_function <- (3/2)*(log(m)*(1-B[, "omega"]) +log(m+1)*B[, "omega"] ) +
                   (1/2)*log((1-m)/m^2) - ((1-m)/(2*m^2))*(m*B[,"Z.0"]*(1-B[, "omega"]) + (m+1)*B[,"Z.1"]*B[, "omega"] -
                     2*((m^2)*(1-B[, "omega"]) + ((m+1)^2)*B[, "omega"])  +
                      (m^3)*B[,"kappa.0"]*(1-B[, "omega"]) +((m+1)^3)*B[,"kappa.1"]*B[, "omega"]     ) +
                       Y*log(1-m)+(1-Y)*log(m)
    q3_function=-sum(q3_function)
    return(q3_function)
#-MNBcr case   Case----------------------------------------------------------------------
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
  } else if (dist == 4) { # Birnbaum-Saunders
    S    <- pbisa(t, shape = alpha, scale = 1/nu, lower.tail = FALSE)
    logf <- dbisa(t, shape = alpha, scale = 1/nu, log = TRUE)
  } else {
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




#===============================================================================
llikeobserved2<-function(psi,Data,dist,model,real=FALSE)
  #-------------------------------------------------------------------------------
{
  t     <- Data$t
  delta <- Data$delta
  x     <- matrix(Data$x,nrow=r)
  beta  <- matrix(psi[1:r],ncol=1)
  alpha <- psi[r+1]
  nu    <- exp(psi[r+2])
  m     <- plogis(psi[r+3])
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
  } else if (dist == 4) { # Birnbaum-Saunders
    S    <- pbisa(t, shape = alpha, scale = 1/nu, lower.tail = FALSE)
    logf <- dbisa(t, shape = alpha, scale = 1/nu, log = TRUE)
  } else {
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


#===============================================================================
llikeobserved_PO <- function(psi, data,  dist) {
#===============================================================================
  t <- data$t
  delta <- data$delta
  z <- matrix(data$x, nrow = r)
  beta <- matrix(psi[1:r], ncol = 1)
  alpha <- psi[r + 1]
  sigma <- psi[r + 2]
  
  if (dist == 1) {
    F    <- 2 * pnorm(t/sigma) - 1 - exp((1/2) * (alpha * 
       log(2) - log(pi)) + alpha * (log(sigma) - log(t)) + 
      lgamma((alpha + 1)/2) + pgamma(t^2/(2 * sigma^2), 
            shape = (alpha + 1)/2, rate = 1, log.p = TRUE))
    logf <- log(alpha) + (1/2) * (alpha * log(2) - log(pi)) + 
      alpha * log(sigma) + lgamma((alpha + 1)/2) - 
      (alpha + 1) * log(t) + pgamma(t^2/(2 * sigma^2), 
                                    shape = (alpha + 1)/2, rate = 1, log.p = TRUE)
  }
  if (dist == 2) {
    F <- 1 - exp(-(t/sigma)^alpha)
    logf <- log(alpha) - log(sigma) + (alpha - 1) * (log(t) - 
                                                       log(sigma)) - (t/sigma)^alpha
  }
  if (dist == 3) {
    F <- pgamma(t, shape = alpha, scale = sigma)
    logf <- dgamma(t, shape = alpha, scale = sigma, log = TRUE)
  }
  if (dist == 4) {
    F <- pbisa(t, shape = alpha, scale = sigma)
    logf <- dbisa(t, shape = alpha, scale = sigma, log = TRUE)
  }
  if (dist == 5) {
    F <- plnorm(t, meanlog = alpha, sdlog = sigma)
    logf <- dlnorm(t, meanlog = alpha, sdlog = sigma, 
                   log = TRUE)
  }
  S <- 1 - F
  theta = exp(t(z) %*% beta)
  mu = theta * S
  log.A = theta
  log.AS = mu
  log.fpop = log(theta) + logf - theta * F
  
  
  log.Spop = log.AS - log.A
  -sum(delta * log.fpop + (1 - delta) * log.Spop)
}


#===============================================================================
S_f_hpop<-function(psi,Data,dist,model,real=FALSE)
#-------------------------------------------------------------------------------
{
  t     <- Data$t
  delta <- Data$delta
  x     <- matrix(Data$x,nrow=r)
  beta  <- matrix(psi[1:r],ncol=1)
  alpha <- psi[r+1]
  nu    <- psi[r+2]
  m     <- psi[r+3]
  #-----------------------------------------------------------------------------
  # survival function for W (dist)
  #-----------------------------------------------------------------------------
  if (dist == 1) {                 # Weibull (shape=nu, scale=exp(-alpha/nu))
    S <- exp(-exp(alpha) * t^nu)
    ft <- dweibull(t, shape = nu, scale = exp(-alpha/nu))
    
  } else if (dist == 2) {          # Log-normal
    S <- plnorm(t, meanlog = alpha, sdlog = nu, lower.tail = FALSE)
    ft <- dlnorm(t, meanlog = alpha, sdlog = nu)
    
  } else if (dist == 3) {          # Gamma (shape=alpha, rate=nu)
    S <- pgamma(t, shape = alpha, rate = nu, lower.tail = FALSE)
    ft <- dgamma(t, shape = alpha, rate = nu)
    
  } else if (dist == 4) {          # Birnbaum-Saunders
    S <- pbisa(t, shape = alpha, scale = 1/nu, lower.tail = FALSE)
    ft <- dbisa(t, shape = alpha, scale = 1/nu)
    
  }
  #-----------------------------------------------------------------------------
  theta <- exp(t(x)%*%beta)  #parameter from Poisson (on mixture)
  #------------------------
  
  #-MIG case--------------------------------------------------------------------  
  if (model == 1) {   
  #-----------------------------------------------------------------------------
  aux1 =  exp((1 - m) * (1 - sqrt(1 -(2 * m * theta*(S-1)) / (1 - m)))) #modifiquei o + dentro da raiz para -
  aux2 =  exp(((1 + m) * (1 - m^2))/m^2 * (1 - sqrt(1 - (2 * m^2 * theta*(S-1)) / (1 - m^2))))
  Spop <-  m*aux1 + (1 - m) *aux2 
  fpop <- theta*ft * ( (m^2/sqrt(1- (2*m*theta*(S-1))/(1-m)))*aux1  + ((1-m^2)/sqrt(1- (2*m^2*theta*(S-1))/(1-m^2)) )*aux2       )
  hpop <- fpop/Spop
  return(list(Spop = Spop, fpop = fpop, hpop = hpop))
  #-MGA case--------------------------------------------------------------------
  } else if (model == 2) { 
  #-----------------------------------------------------------------------------
    aux1 <- (1- (m*theta*(S-1))/(1-m))
    aux2 <- (1- ((m^2)*theta*(S-1))/(1-m^2))
    Spop <-  m*aux1^(-(1-m))  + (1-m)*aux2^(- ( (((m+1)^2)*(1-m))/m^2 ) )  
    fpop <-  theta*ft*( m^2 *aux1^(-(1-m)-1) + (1-m^2)*aux2^(- ( (((m+1)^2)*(1-m))/m^2 ) - 1 )  )
    hpop <- fpop/Spop
    return(list(Spop = Spop, fpop = fpop, hpop = hpop))
  } else {
    stop("Modelo invalido: use 1 (MIG) ou 2 (MGA).")
  }
  #------------------------
}#end Spop

#===============================================================================





  
  
  
  
  
  
  
  

#===============================================================================
Spop.PO <-function(psi_PO,Data,model,dist,real=FALSE)
#-------------------------------------------------------------------------------
{
  t     <- Data$t
  delta <- Data$delta
  x     <- matrix(Data$x,nrow=r)
  beta  <- matrix(psi_PO[1:r],ncol=1)
  alpha <- psi_PO[r+1]
  sigma <- psi_PO[r+2]
  #-----------------------------------------------------------------------------
  # survival function and log density for W (dist)
  #-----------------------------------------------------------------------------
  if (dist == 1) {
    F    <- 2 * pnorm(t/sigma) - 1 - exp((1/2) * (alpha * 
          log(2) - log(pi)) + alpha * (log(sigma) - log(t)) + 
          lgamma((alpha + 1)/2) + pgamma(t^2/(2 * sigma^2), 
          shape = (alpha + 1)/2, rate = 1, log.p = TRUE))
    logf <- log(alpha) + (1/2) * (alpha * log(2) - log(pi)) + 
         alpha * log(sigma) + lgamma((alpha + 1)/2) - 
         (alpha + 1) * log(t) + pgamma(t^2/(2 * sigma^2), 
         shape = (alpha + 1)/2, rate = 1, log.p = TRUE)
  }
  if (dist == 2) {
    F    <- 1 - exp(-(t/sigma)^alpha)
    logf <- log(alpha) - log(sigma) + (alpha - 1) * (log(t) - log(sigma)) - (t/sigma)^alpha
  }
  if (dist == 3) {
    F    <- pgamma(t, shape = alpha, scale = sigma)
    logf <- dgamma(t, shape = alpha, scale = sigma, log = TRUE)
  }
  if (dist == 4) {
    F    <- pbisa(t, shape = alpha, scale = sigma)
    logf <- dbisa(t, shape = alpha, scale = sigma, log = TRUE)
  }
  if (dist == 5) {
    F    <- plnorm(t, meanlog = alpha, sdlog = sigma)
    logf <- dlnorm(t, meanlog = alpha, sdlog = sigma, log = TRUE)
  }
  S <- 1 - F
  #-----------------------------------------------------------------------------
  if (model == 1 || model == 4) {
    theta = exp(t(x) %*% beta)
  }
  if (model == 2 || model == 3 || model == 5 || model == 
      6) {
    theta = plogis(t(x) %*% beta)
  }
  if (model == 1) {
    mu = theta * S
    log.A = theta
    log.AS = mu
  }
  if (model == 2) {
    mu = theta * S
    log.A = log(-log(1 - theta)) - log(theta)
    log.AS = log(-log(1 - mu)) - log(mu)
  }
  if (model == 3) {
    mu = theta * S
    log.A = -q * log1p(-theta)
    log.AS = -q * log1p(-mu)
  }
  if (model == 4) {
    mu = theta * S
    log.A = q * log(1 + theta)
    log.AS = q * log(1 + mu)
  }
  if (model == 5) {
    mu = theta * S
    log.A = log(polyloga(theta, q)) + log(S)
    log.AS = log(polyloga(mu, q))
  }
  if (model == 6) {
    mu = theta * S
    log.A = -2 * log1p(-theta)
    log.AS = -2 * log1p(-mu)
  }
  log.Spop = log.AS - log.A
  exp(log.Spop)
}#end Spop_PO
#===============================================================================
Spop.BN<-function(psi_BN,Data,model,dist=2,real=FALSE)
  #-------------------------------------------------------------------------------
{
  t     <- Data$t
  delta <- Data$delta
  x     <- matrix(Data$x,nrow=r)
  beta  <- matrix(psi_BN[1:r],ncol=1)
  alpha <- psi_BN[r+1]
  sigma <- psi_BN[r+2]
  q     <- psi_BN[r+3]
  #-----------------------------------------------------------------------------
  # survival function and log density for W (dist)
  #-----------------------------------------------------------------------------
  if (dist == 1) {
    F    <- 2 * pnorm(t/sigma) - 1 - exp((1/2) * (alpha * 
              log(2) - log(pi)) + alpha * (log(sigma) - log(t)) + 
             lgamma((alpha + 1)/2) + pgamma(t^2/(2 * sigma^2), 
             shape = (alpha + 1)/2, rate = 1, log.p = TRUE))
    logf <- log(alpha) + (1/2) * (alpha * log(2) - log(pi)) + 
             alpha * log(sigma) + lgamma((alpha + 1)/2) - 
            (alpha + 1) * log(t) + pgamma(t^2/(2 * sigma^2), 
            shape = (alpha + 1)/2, rate = 1, log.p = TRUE)
  }
  if (dist == 2) {
    F    <- 1 - exp(-(t/sigma)^alpha)
    logf <- log(alpha) - log(sigma) + (alpha - 1) * (log(t) - log(sigma)) - (t/sigma)^alpha
  }
  if (dist == 3) {
    F    <- pgamma(t, shape = alpha, scale = sigma)
    logf <- dgamma(t, shape = alpha, scale = sigma, log = TRUE)
  }
  if (dist == 4) {
    F    <- pbisa(t, shape = alpha, scale = sigma)
    logf <- dbisa(t, shape = alpha, scale = sigma, log = TRUE)
  }
  if (dist == 5) {
    F    <- plnorm(t, meanlog = alpha, sdlog = sigma)
    logf <- dlnorm(t, meanlog = alpha, sdlog = sigma, log = TRUE)
  }
  S <- 1 - F
  #-----------------------------------------------------------------------------
  if (model == 1 || model == 4) {
    theta = exp(t(x) %*% beta)
  }
  if (model == 2 || model == 3 || model == 5 || model == 
      6) {
    theta = plogis(t(x) %*% beta)
  }
  if (model == 1) {
    mu = theta * S
    log.A = theta
    log.AS = mu
  }
  if (model == 2) {
    mu = theta * S
    log.A = log(-log(1 - theta)) - log(theta)
    log.AS = log(-log(1 - mu)) - log(mu)
  }
  if (model == 3) {
    mu = theta * S
    log.A = -q * log1p(-theta)
    log.AS = -q * log1p(-mu)
  }
  if (model == 4) {
    mu = theta * S
    log.A = q * log(1 + theta)
    log.AS = q * log(1 + mu)
  }
  if (model == 5) {
    mu = theta * S
    log.A = log(polyloga(theta, q)) + log(S)
    log.AS = log(polyloga(mu, q))
  }
  if (model == 6) {
    mu = theta * S
    log.A = -2 * log1p(-theta)
    log.AS = -2 * log1p(-mu)
  }
  log.Spop = log.AS - log.A
  exp(log.Spop)
}#end Spop_BN

#===============================================================================
Spop.Bern<-function(psi_Bern,Data,model=4,dist=2,real=FALSE)
  #-------------------------------------------------------------------------------
{
  t     <- Data$t
  delta <- Data$delta
  x     <- matrix(Data$x,nrow=r)
  beta  <- matrix(psi_Bern[1:r],ncol=1)
  alpha <- psi_Bern[r+1]
  sigma <- psi_Bern[r+2]
  q = 1
  #-----------------------------------------------------------------------------
  # survival function and log density for W (dist)
  #-----------------------------------------------------------------------------
  if (dist == 1) {
    F    <- 2 * pnorm(t/sigma) - 1 - exp((1/2) * (alpha * 
                log(2) - log(pi)) + alpha * (log(sigma) - log(t)) + 
            lgamma((alpha + 1)/2) + pgamma(t^2/(2 * sigma^2), 
            shape = (alpha + 1)/2, rate = 1, log.p = TRUE))
    logf <- log(alpha) + (1/2) * (alpha * log(2) - log(pi)) + 
      alpha * log(sigma) + lgamma((alpha + 1)/2) - 
      (alpha + 1) * log(t) + pgamma(t^2/(2 * sigma^2), 
         shape = (alpha + 1)/2, rate = 1, log.p = TRUE)
  }
  if (dist == 2) {
    F    <- 1 - exp(-(t/sigma)^alpha)
    logf <- log(alpha) - log(sigma) + (alpha - 1) * (log(t) - log(sigma)) - (t/sigma)^alpha
  }
  if (dist == 3) {
    F    <- pgamma(t, shape = alpha, scale = sigma)
    logf <- dgamma(t, shape = alpha, scale = sigma, log = TRUE)
  }
  if (dist == 4) {
    F    <- pbisa(t, shape = alpha, scale = sigma)
    logf <- dbisa(t, shape = alpha, scale = sigma, log = TRUE)
  }
  if (dist == 5) {
    F    <- plnorm(t, meanlog = alpha, sdlog = sigma)
    logf <- dlnorm(t, meanlog = alpha, sdlog = sigma, log = TRUE)
  }
  S <- 1 - F
  #-----------------------------------------------------------------------------
  if (model == 1 || model == 4) {
    theta = exp(t(x) %*% beta)
  }
  if (model == 2 || model == 3 || model == 5 || model == 
      6) {
    theta = plogis(t(x) %*% beta)
  }
  if (model == 1) {
    mu = theta * S
    log.A = theta
    log.AS = mu
  }
  if (model == 2) {
    mu = theta * S
    log.A = log(-log(1 - theta)) - log(theta)
    log.AS = log(-log(1 - mu)) - log(mu)
  }
  if (model == 3) {
    mu = theta * S
    log.A = -q * log1p(-theta)
    log.AS = -q * log1p(-mu)
  }
  if (model == 4) {
    mu = theta * S
    log.A = q * log1p( theta)
    log.AS = q * log1p(mu)
  }
  if (model == 5) {
    mu = theta * S
    log.A = log(polyloga(theta, q)) + log(S)
    log.AS = log(polyloga(mu, q))
  }
  if (model == 6) {
    mu = theta * S
    log.A = -2 * log1p(-theta)
    log.AS = -2 * log1p(-mu)
  }
  log.Spop = log.AS - log.A
  exp(log.Spop)
}#end Spop_BERN


#  hessian(llikeobserved,x0=psi,Data=Data,dist=1)










#=========================================================
est.NA<-function(t,delta)
{
  sfitall = survfit(Surv(t, delta)~1)
  c = coxph(formula=Surv(t, delta)~1)
  n=basehaz(c)
  cbind(n[,2],n[,1])
}

RQ.res<-function(Spop,nrep=1000)
{
  nrep = 1000
  mqresid = NULL
  u = delta * (1 - Spop) + (1 - delta) * runif(length(t), 1 - Spop)
  #for (i in 1:nrep) 
  #{
  #	qresid = sort(qnorm(runif(u)))
  #	mqresid = cbind(mqresid, qresid)
  #}
  #qresid = apply(mqresid, 1, median)
  qnorm(u)
}

#########################################################
## QQ plot - Quantile - Residuals
#########################################################
envelopeDS <- function(x){
  U	         <- x
  n	         <- length(x)
  d2s 	     <- sort(U)
  xq2 	     <- qnorm(ppoints(n))
  Xsim 	     <- matrix(0, 100, n)
  for(i in 1:100){
    u2       <- rnorm(n)
    Xsim[i,] <- u2
  }
  Xsim2      <- apply(Xsim, 1, sort)
  d21        <- matrix(0, n, 1)
  d22        <- matrix(0, n, 1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,], 0.025)
    d22[i]  <- quantile(Xsim2[i,], 0.975)
  }
  d2med      <- apply(Xsim2, 1, mean)
  fy         <- range(d2s, d21, d22)
  plot(xq2, d2s, xlab = quote("Theoretical quantile"),
       ylab = quote("Empirical quantile"), 
       pch = 20, ylim = fy)
  par(new = T)
  plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "")
  par(new = T)
  plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "")
  par(new = T)
  plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "")
}


#########################################################
## QQ plot - Cox-Snell - Residuals
######################################################### 

envelopeCS <- function(x){
  U	         <- x
  n	         <- length(x)
  d2s 	     <- sort(U)
  xq2 	     <- qexp(ppoints(n))
  Xsim 	     <- matrix(0, 100, n)
  for(i in 1:100){
    u2       <- rexp(n)
    Xsim[i,] <- u2
  }
  Xsim2      <- apply(Xsim, 1, sort)
  d21        <- matrix(0, n, 1)
  d22        <- matrix(0, n, 1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,], 0.025)
    d22[i]  <- quantile(Xsim2[i,], 0.975)
  }
  d2med      <- apply(Xsim2, 1, mean)
  fy         <- range(d2s, d21, d22)
  plot(xq2, d2s, xlab = quote("qe"),
       ylab = quote("qr"), 
       pch = 20, ylim = fy)
  par(new = T)
  plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "")
  par(new = T)
  plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "")
  par(new = T)
  plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "")
}



#===============================================================================
passoE1<-function(psi,Data,dist,model)
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
    
    logBsab0p =   log(besselK(sab0,p,expon.scaled = TRUE)) - sab0
    logBsab0p1 =   log(besselK(sab0,p+1,expon.scaled = TRUE)) - sab0
    logBsab1p =   log(besselK(sab1,p,expon.scaled = TRUE)) - sab1
    logBsab1p1 =   log(besselK(sab1,p+1,expon.scaled = TRUE)) - sab1
    
    aux1_omega =  (1-m)*exp(((m+1)*(1-m^2))/(m^2))*((m+1)^(3/2))*(
      exp(logBsab1p- log((ay1/by1)^(p/2)) ) )
     
    #---------------------------------------------------------------------------
    omega = aux1_omega/(aux0_omega + aux1_omega)
    #--Conditional expectations-------------------------------------------------
    Y     = omega
    Z = (sqrt(by0)/sqrt(ay0))*exp(  logBsab0p1 -   logBsab0p )*(1-omega) + 
      (sqrt(by1)/sqrt(ay1))*exp( logBsab1p1 -  logBsab1p  )*omega 
    N     = delta + theta*S*Z 

    Z.0 = (sqrt(by0)/sqrt(ay0))*exp(  logBsab0p1 -   logBsab0p )
    Z.1 = (sqrt(by1)/sqrt(ay1))*exp( logBsab1p1 -  logBsab1p  )
    
    kappa.0  = (sqrt(ay0)/sqrt(by0))*exp(  logBsab0p1 -   logBsab0p ) - ((2*p)/by0)
    kappa.1  = (sqrt(ay1)/sqrt(by1))*exp( logBsab1p1 -  logBsab1p  )  - ((2*p)/by1)
    
    B = cbind(omega, Z.0, Z.1 , kappa.0, kappa.1) - ((2*p)/by1)
    colnames(B)<-c("omega","Z.0","Z.1","kappa.0","kappa.1")
  }
  #-MGA CASE (mix=2)------------------------------------------------------------
  if(model==2){
    #--auxiliar values for conditional expectations-----------------------------
    aux0_nu = m*( (((1-m)^delta)*(((1-m)/m)^(1-m))) /( ((1-m)/m + theta*(1-S))^(1-m+delta) ) )
    aa = (((m+1)*(1-m^2))/m^2 + delta)
    bb = ((1-m^2)/(m^2))
    cc = theta*(1-S)
    aux1_nu = (1-m)*( ((m+1)^delta) * ( (  bb/ ( bb+ cc ))^aa    )) 
    
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




################################################################################ 
# data generator for MPIGcr or MNBcr
################################################################################ 

data.generator<-function(n,beta,alpha,nu,m,model,max.censura)
{
  r=length(beta.real)
  beta<-matrix(beta.real,ncol=1)
  u1 = sample(1:4, size=n, replace=TRUE, prob=c(0.25,0.37,0.27, 1-(0.25+0.37+0.27)))
  u2 = sample(1:2, size=n, replace=TRUE, prob=c(0.23,0.77))
  u3 = sample(1:2, size=n, replace=TRUE, prob=c(0.55,0.45))
  u4 = sample(1:2, size=n, replace=TRUE, prob=c(0.33,0.67))
  id = round(rnorm(n, mean = 56.3, sd = 13.62))
  id = pmax(12, pmin(103, id))
  x<-t(model.matrix(~as.factor(u1) + as.factor(u2) + as.factor(u3) + as.factor(u4) + id  ))
  #--------------------------------------- 
  #  regression estructure on theta
  theta=exp(t(x)%*%beta) #parameter of Poisson (on mixture)
  #------------------------------------------
  ## model=1 MGA; model=2 MIG
  RIG<-function(n, m, gamma) rinvgauss(n, mean=m, shape=m^3/gamma)
  RGA<-function(n, m, gamma) rgamma(n, shape=m^2/gamma, rate=m/gamma)
  D<-get(ifelse(model==1,"RIG","RGA"))
  gamma<-m^2/(1-m)
  p=gamma/(gamma+m)
  Y<-rbinom(n, 1, 1-p)
  Z<-D(n, m+Y, gamma)
  N<-rpois(n, theta*Z)
  #---------------------------------------------- 
  W<-rep(Inf,n)  
  for(i in 1:n)
  {
    if(N[i]>0) 
    {
      U<-runif(N[i])
      W[i]<-(exp(-alpha)*min(-log(U)))^(1/nu)
    }
  }
  t<-c();
  delta<-rep(1,n) 
  for(i in 1:n)
  {
    Cens = max.censura  
    t[i]=Cens;
    delta[i]<-0
    if(W[i]<=Cens)
    {
      t[i]=W[i];
      delta[i]=1
    }
  }
  list(t=t,delta=delta,x=x)
}


################################################################################ 
# data generator for Poisson 
################################################################################ 
data.generator_PTcr<-function(n,beta,alpha,nu,max.censura)
{
  r=length(beta.real)
  beta<-matrix(beta.real,ncol=1)
  u1 = sample(1:4, size=n, replace=TRUE, prob=c(0.25,0.37,0.27, 1-(0.25+0.37+0.27)))
  u2 = sample(1:2, size=n, replace=TRUE, prob=c(0.23,0.77))
  u3 = sample(1:2, size=n, replace=TRUE, prob=c(0.55,0.45))
  u4 = sample(1:2, size=n, replace=TRUE, prob=c(0.33,0.67))
  id = round(rnorm(n, mean = 56.3, sd = 13.62))
  id = pmax(12, pmin(103, id))
  x<-t(model.matrix(~as.factor(u1) + as.factor(u2) + as.factor(u3) + as.factor(u4) + id  ))
  #--------------------------------------- 
  #  regression estructure on theta
  theta=exp(t(x)%*%beta) #parameter of Poisson (on mixture)
  #------------------------------------------
 
  N<-rpois(n, theta)
  #---------------------------------------------- 
  W<-rep(Inf,n)  
  for(i in 1:n)
  {
    if(N[i]>0) 
    {
      U<-runif(N[i])
      W[i]<-(exp(-alpha)*min(-log(U)))^(1/nu)
    }
  }
  t<-c();
  delta<-rep(1,n) 
  for(i in 1:n)
  {
    Cens = max.censura  
    t[i]=Cens;
    delta[i]<-0
    if(W[i]<=Cens)
    {
      t[i]=W[i];
      delta[i]=1
    }
  }
  list(t=t,delta=delta,x=x)
}



