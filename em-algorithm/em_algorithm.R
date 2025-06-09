library(dplyr)
library(rmutil)
library(Rmpfr)

##########################################################################################
########################## Define E[h(u_i) | Y_i, T_i] terms ############################# 
##########################################################################################
E_h_step <- function(Data, Indv , time=TRUE, gamma, nu2, tau2, alpha, beta){
  
  # Gauss-Hermite quadrature points and weights with 2 points
  gg = gauss.hermite(2)
  zz = gg[,1] # quadrature nodes (points)
  ww = gg[,2] # quadrature weights
  
  # Select data for individual Indv
  ii = Data$id %in% Indv
  
  # Check if individual exists in the data
  if(sum(ii) == 0){stop("there is no this ID in the data") }
  
  y_i = Data$y[ii] # response vector for individual i
  m_i <- sum(ii) # number of observations for individual i
  
  t_i = Data$t[ii] # time points for individual i
  r_i <- diff(t_i) # time differences between observations
  r_i <- c(t_i[1],r_i) # include first time point
  
  # Choose columns to exclude based on whether 'time' is considered or not
  ifelse(time==FALSE, dro <- c("y","id","t"), dro <- c("y","id"))
  cov_i <- as.matrix(Data[ii, !(names(Data) %in% dro),drop=F],nrow=length(y_i)) # covariate matrix
  
  # Design matrix with intercept plus covariates
  x_i <- matrix(c(rep(1,m_i),cov_i),nrow=m_i,ncol=ncol(cov_i)+1) # Design matrix with intercept + covariates
  
  # Calculate conditional variance component dd_i
  dd_i <-  nu2*t(rep(1,m_i)) %*% solve((nu2*matrix(rep(1),ncol=m_i,nrow=m_i))+(tau2*diag(m_i)))
  
  var_i <- as.numeric( nu2 - dd_i %*% (nu2*rep(1,m_i)) )# variance term
  
  # Initialize numerator and denominator sums for expected values
  E_h1_ui_num = 0
  E_h2_ui_num = 0
  E_h3_ui_num = 0
  E_h_ui_denom = 0
  
  # Loop over quadrature points for numerical integration
  for(k in 1:(length(zz))){
    
    f1_i <- c()
    f2_i <- c()
    
    # Loop over observations for individual i
    for(j in 1:m_i){
      
      f1_ij <- as.numeric( exp(alpha+beta*sqrt(2*var_i)*zz[k]-r_i[j]*exp(alpha+beta*sqrt(2*var_i)*zz[k])) )
      
      f2_ij <- as.numeric( exp( (-1/(2*var_i)) *( y_i[j]^2 - 2*y_i[j]*(x_i[j,]%*%gamma + sqrt(2*var_i)*zz[k]) + (x_i[j,]%*%gamma)^2 + 2*x_i[j,]%*%gamma*sqrt(2*var_i)*zz[k] )))
      
      f1_i <- c(f1_i, f1_ij)
      f2_i <- c(f2_i, f2_ij)
    }
    
    # Accumulate weighted sums for numerator and denominator of expectations
    E_h1_ui_num = E_h1_ui_num + ww[k]*sqrt(2*var_i)*zz[k]*prod(f1_i)*ifelse(prod(f2_i)==0,1e-100,prod(f2_i))
    E_h2_ui_num = E_h2_ui_num + ww[k]*((sqrt(2*var_i)*zz[k])^2)*prod(f1_i)*ifelse(prod(f2_i)==0,1e-100,prod(f2_i))
    E_h3_ui_num = E_h3_ui_num + ww[k]*exp(beta*sqrt(2*var_i)*zz[k])*prod(f1_i)*ifelse(prod(f2_i)==0,1e-100,prod(f2_i))
    
    E_h_ui_denom = E_h_ui_denom + ww[k]*prod(f1_i)*ifelse(prod(f2_i)==0,1e-100,prod(f2_i))
  }
  
  # Calculate expected values dividing numerator by denominator
  E_h1_ui = E_h1_ui_num / E_h_ui_denom
  E_h2_ui = E_h2_ui_num / E_h_ui_denom
  E_h3_ui = E_h3_ui_num / E_h_ui_denom
  
  # Return list of expectations
  return(list(E_h1_ui,E_h2_ui,E_h3_ui))
}

##########################################################################################
######################## Define negative expected log-likelihood ######################### 
##########################################################################################
neg_E_loglik <- function(dd, Time=TRUE, parameters = list(gamma , nu2, tau2, alpha, beta )){
  
  # Extract variables from data
  y <- as.matrix(dd$y)
  t <- as.matrix(dd$t)
  id <- as.matrix(dd$id)
  ii <- unique(id)
  
  # Drop appropriate columns based on 'Time' argument
  ifelse(Time==FALSE, dro <- c("y","id","t"), dro <- c("y","id"))
  cov <- as.matrix(dd[, !(names(dd) %in% dro),drop=F],nrow=length(y))
  
  n <- length(ii)# Number of unique individuals
  N <- length(y) # Total number of observations
  
  # Design matrix with intercept and covariates
  x <- matrix(c(rep(1,N),cov),nrow=N,ncol=ncol(cov)+1)
  
  # Extract parameters from input list
  gg <- parameters[[1]]
  nn2 <- parameters[[2]]
  tt2 <- parameters[[3]]
  aa <- parameters[[4]]
  bb <- parameters[[5]]
  
  ffloglik <- 0
  E_h1_step <- c()
  E_h2_step <- c()
  E_h3_step <- c()
  
  # Loop over individuals
  for(i in 1:n){
    
    ID = ii[i]
    
    jj = id == ID
    
    # Compute expected values E[h(u_i)|data] for individual i
    E_h_step_i = E_h_step(Data=dd, Indv = ID , time=Time, gamma =gg, nu2=nn2, tau2=tt2, alpha=aa, beta=bb)
    
    E_h1_step_i = E_h_step_i[[1]]
    E_h2_step_i = E_h_step_i[[2]]
    E_h3_step_i = E_h_step_i[[3]]
    
    m_i <- length(id[jj])
    x_i <- matrix(x[jj],ncol=ncol(x))
    y_i <- y[jj]
    t_i <- t[jj]
    r_i <- diff(t_i)
    r_i <- c(t_i[1],r_i)
    
    ffloglik_i <- 0
    
    # Loop over observations within individual i
    for(j in 1:m_i){
      
      vv = y_i[j] - x_i[j,]%*%gg
      
      # Calculate individual contribution to expected log-likelihood
      ffloglik_i = ffloglik_i + aa + bb * E_h1_step_i - r_i[j]*exp(aa) * E_h3_step_i - log(sqrt(tt2*2*pi)) - ((vv)^2 + E_h2_step_i - 2*vv*E_h1_step_i)/(2*tt2)
      
    }
    
    # Accumulate total expected log-likelihood
    ffloglik = ffloglik -log(sqrt(nn2*2*pi))-E_h2_step_i/(2*nn2) + ffloglik_i
    
    # Store expected values for diagnostics or further use
    E_h1_step <- c(E_h1_step,E_h1_step_i)
    E_h2_step <- c(E_h2_step,E_h2_step_i)
    E_h3_step <- c(E_h3_step,E_h3_step_i)
    
  }
  
  loglik_optim <- -ffloglik # Negative log-likelihood to minimize
  
  return(loglik_optim)
  
}

##########################################################################################
########################## Define parameter estimators ################################### 
##########################################################################################
param_hat = function(DATA, TIME=TRUE, init=list(GAMMA, NU2, TAU2, ALPHA, BETA)){
  
  # Extract variables
  y <- as.matrix(DATA$y) # response variable vector
  t <- as.matrix(DATA$t) # time points vector
  id <- as.matrix(DATA$id) # individual IDs vector
  ii <- unique(id) # unique individuals
  
  # Drop appropriate columns based on 'TIME' argument
  ifelse(TIME==FALSE, dro <- c("y","id","t"), dro <- c("y","id"))
  cov <- as.matrix(DATA[, !(names(DATA) %in% dro),drop=F],nrow=length(y))
  
  n <- length(ii)# Number of unique individuals
  N <- length(y)# Total number of observations
  
  # Design matrix with intercept + covariates
  x <- matrix(c(rep(1,N),cov),nrow=N,ncol=ncol(cov)+1)
  
  # Initialize parameter estimates
  nu2_hat = 0
  tau2_hat = 0
  alpha_hat = 0
  
  E_U1 = c() # empty vector to store expected latent variable
  
  # Loop over individuals to calculate estimators
  for(i in 1:n){
    
    ID = ii[i]
    
    jj = id == ID
    
    # Calculate expected values via function E_h_step
    E_h_step_i = E_h_step(Data=DATA, Indv = ID , time=TIME, gamma =init[[1]], nu2=init[[2]], tau2=init[[3]], alpha=init[[4]], beta=init[[5]])
    
    E_h1_step_i = E_h_step_i[[1]]
    E_h2_step_i = E_h_step_i[[2]]
    E_h3_step_i = E_h_step_i[[3]]
    
    m_i <- length(id[jj])
    x_i <- matrix(x[jj],ncol=ncol(x))
    y_i <- y[jj]
    t_i <- t[jj]
    r_i <- diff(t_i)
    r_i <- c(t_i[1],r_i)
    
    tau2_hat_i = 0
    alpha_hat_i = 0
    
    # Loop over observations of individual i
    for(j in 1:m_i){
      
      vv = y_i[j] - x_i[j,]%*%init[[1]]
      
      tau2_hat_i = tau2_hat_i + as.numeric((vv)^2 + E_h2_step_i - 2*vv*E_h1_step_i)
      
      alpha_hat_i = alpha_hat_i + r_i[j]*E_h3_step_i
      
    }
    
    nu2_hat = nu2_hat + E_h2_step_i
    
    tau2_hat = tau2_hat + tau2_hat_i
    
    E_U1 = c(E_U1 , rep(E_h1_step_i , m_i))
    
    alpha_hat = alpha_hat + alpha_hat_i
    
  }
  
  gamma_hat = solve(t(x)%*%x)%*%t(x)%*%(y - E_U1)
  nu2_hat = nu2_hat/n 
  tau2_hat = as.numeric(tau2_hat/N)
  alpha_hat = log(N/alpha_hat)
  
  # Negative expected log-likelihood function for optimisation of beta
  negEloglik = function(BETA){
    return(as.numeric(neg_E_loglik(dd = DATA , Time=TIME, parameters = list(gamma=init[[1]], nu2=init[[2]], tau2=init[[3]], alpha=init[[4]], beta=BETA ))[[1]]))
  }
  
  beta_hat = optim(par=init[[5]], fn=negEloglik, method="CG" )$par # Optimise beta 
  
  # Return updated parameter estimates and design matrix
  return(list(gamma_new=gamma_hat, nu2_new=nu2_hat , tau2_new=tau2_hat , alpha_new=alpha_hat, beta_new=beta_hat, E_U=E_U1, X=x))
  
}

############################################################################################
######## Define approximation to the integral within the likelihood per individual #########
############################################################################################
approx <- function(Dat, time=TRUE, Indiv , Gamma, Nu2, Tau2, Alpha, Beta){
  
  GG = gauss.hermite(2) # Gauss-Hermite quadrature points and weights
  ZZ = GG[,1] # Quadrature nodes
  WW = GG[,2] # Quadrature weights
  
  ii = Dat$id %in% Indiv
  if(sum(ii) == 0){stop("there is no this ID in the data") }
  
  y_i = Dat$y[ii]
  m_i <- sum(ii)
  
  t_i = Dat$t[ii]
  r_i <- diff(t_i)
  r_i <- c(t_i[1],r_i)
  
  # Determine columns to exclude depending on 'time' 
  ifelse(time==FALSE, dro <- c("y","id","t"), dro <- c("y","id"))
  cov_i <- as.matrix(Dat[ii, !(names(Dat) %in% dro),drop=F],nrow=length(y_i))
  
  # Design matrix with intercept and covariates for individual i
  x_i <- matrix(c(rep(1,m_i),cov_i),nrow=m_i,ncol=ncol(cov_i)+1)
  
  integrand_k <- 0 # Initialise integral sum
  
  # Loop over Gauss-Hermite nodes
  for(k in 1:(length(ZZ))){
    
    mm_i <- c() # Initialise vector to store product terms
    
    # Loop over observations for individual i
    for(j in 1:m_i){
      
      f1_ij <- as.numeric( exp(Alpha+Beta*sqrt(2*Nu2)*ZZ[k]-r_i[j]*exp(Alpha+Beta*sqrt(2*Nu2)*ZZ[k])) )
      
      f2_ij <- as.numeric( exp(-(y_i[j]-x_i[j,]%*%Gamma-sqrt(2*Nu2)*ZZ[k])^2/(2*Tau2)) )
      
      mm_ij <- 1/sqrt(2*pi*Tau2) * f2_ij * f1_ij * 1/sqrt(2*Nu2)
      
      mm_i <- c(mm_i, mm_ij)
      
    }
    
    # Accumulate weighted product of integrand components for node k
    integrand_k <-  integrand_k + WW[k]/sqrt(2*pi*Nu2)*ifelse(prod(mm_i)==0,1e-100,prod(mm_i))
    
  }
  
  return(integrand_k) # Return approximated integral value
  
}

###############################################################################################
######################### Define the likelihood of the observed data ##########################
###############################################################################################
likelihood_obs <- function(DD, Time=TRUE, param=list(G, N2, T2, A, B)){
  
  id <- as.matrix(DD$id)
  ii=unique(id)
  n <- length(ii)
  
  L <- 1 # Initialise likelihood product
  
  # Loop over all individuals
  for(i in 1:n){
    
    ID = ii[i] # Current individual ID
    
    # Compute approximate integral for individual i using 'approx' function
    L_i = mpfr(approx(Dat=DD, Indiv = ID , time=Time, Gamma=param[[1]], Nu2=param[[2]], Tau2=param[[3]], Alpha=param[[4]], Beta=param[[5]]), precBits=100, scientific=T)
    
    L <- L * L_i # Update likelihood product
    
  }
  
  return(L) # Return overall likelihood value
}

##########################################################################################
################################## Define EM algorithm ###################################
##########################################################################################
EM_alg <- function(Data2, time=TRUE, param_inic = list(gamma.0, nu2.0, tau2.0, alpha.0, beta.0)){
  
  param_ant <- param_inic # Initialise parameters with initial values
  lik_ant <- likelihood_obs(DD=Data2, Time=time, param=param_ant) # Compute initial likelihood
  
  diff_lik <- 100 # Initialise difference in likelihood for convergence check
  
  # Iterate until convergence criterion is met
  while(diff_lik > 0.001){
    
    # Update parameters using parameter estimator function
    param_new <- param_hat(DATA=Data2, TIME=time, init=param_ant)
    
    # Extract updated parameter list
    param_new2 <- list(param_new$gamma_new, param_new$nu2_new, param_new$tau2_new, param_new$alpha_new, param_new$beta_new)
    
    # Calculate new likelihood with updated parameters
    lik_new <- likelihood_obs(DD= Data2, Time=time, param=param_new2)
    
    diff_lik <- abs(lik_new - lik_ant) # Compute absolute difference in likelihood
    
    lik_ant <- lik_new # Update likelihood for next iteration
    
    param_ant <- param_new2 # Update parameters for next iteration
    
  }
  return(list(param_new=param_new, lik_new=lik_new, diff_lik=diff_lik)) # Return final estimates and likelihood
}