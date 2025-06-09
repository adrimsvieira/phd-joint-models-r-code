library(dplyr)
library(rmutil)
library(nlme)

###################### Simulation of the random effect U_i ##########################
Usim<-function(M,by,nusq){
  s<-seq(0,M,by) # Generate a sequence of time points from 0 to M in steps of 'by'
  
  list(U=rnorm(1, mean=0, sd=sqrt(nusq)),S=s) # Simulate a single random effect U ~ N(0, nusq) and return along with time points
}

# Function to simulate longitudinal outcomes given random effect U and time points S
Ysim <- function(U,S,alpha,beta,tausq,gamma1,gamma2){
  t <- NULL # Initialise vector to store observed times
  y <- NULL # Initialise vector to store longitudinal observations
  tau <- sqrt(tausq) # Standard deviation of measurement error
  
  # Define intensity function depending on random effect U
  lambda <- function(alpha,beta,U){
    exp(alpha+beta*U)
  }
  
  for (i in 1:length(S)) {
    lambda.current<-lambda(alpha,beta,U)
    p<-1-exp(-lambda.current)
    x = rbinom(1 , prob= p , size=1)
    if(x==1){
      t <- c(t,S[i])
      # Simulate longitudinal measurement at event time with measurement error and random effect
      y <- c(y,gamma1+gamma2*S[i]+U+rnorm(1,mean=0 , sd = sqrt(tausq)))
    }
  }
  list(T=t,Y=y) # Return event times and corresponding longitudinal measurements
}  

# Function to simulate dataset with n individuals
sim_data_set<- function(n , M , by , nusq , alpha, beta, tausq , gamma1 , gamma2){
  
  dd = data.frame(id=NULL , time = NULL , Y = NULL)
  
  for(l in 1:n){
    
    # Simulate data for individual l
    Y.l = Ysim(Usim(M,by,nusq)$U, Usim(M,by,nusq)$S, alpha, beta , tausq , gamma1, gamma2)
    
    id = rep(l , length(Y.l$Y)) # Replicate individual id for all observations
    
    dd.l = cbind(id ,  y=Y.l$Y, t=Y.l$T ) # Combine id, longitudinal measurements, and times into a dataframe
    
    dd = rbind(dd , dd.l) # Append individual data to overall dataset
  }
  
  list(data =dd ) # Return the full dataset
}

# Function to run simulation study with n.base repetitions
simul_study = function(n.base, n, M, by, nusq, alpha, beta, tausq, gamma1, gamma2){
  
  # Initialise vectors to store estimates from each simulation
  g1hat_vect=NULL
  g2hat_vect=NULL
  nu2hat_vect=NULL
  tau2hat_vect=NULL
  alphahat_vect=NULL
  betahat_vect=NULL
  
  g1hat_vect_lme=NULL
  g2hat_vect_lme=NULL
  nu2hat_vect_lme=NULL
  tau2hat_vect_lme=NULL
  
  ii = 1
  while(ii <=n.base){
    
    print(paste0("n.base ", ii)) # Print progress
 
    # Simulate dataset
    DADOSsimul <- sim_data_set(n, M, by, nusq, alpha, beta, tausq, gamma1, gamma2)$data
    
    # Fit joint model using EM algorithm (defined em_algorithm function)
    mm_fit <- EM_alg(Data2=DADOSsimul, time=TRUE, param_inic = list(gamma.0=c(gamma1, gamma2), nu2.0=nusq, tau2.0=tausq, alpha.0=alpha, beta.0=beta))
    
    # Fit linear mixed-effects model using nlme
    mm_fit_lme <- lme(y~t,random=~1|id, data=DADOSsimul)
    
    if(is.numeric(mm_fit)){ # If EM algorithm returns numeric (error or no convergence), skip this iteration
      }else{
      ii = ii + 1 # Increment simulation counter
      
      # Store parameter estimates from EM algorithm
      g1hat_vect=c(g1hat_vect, mm_fit$param_new$gamma_new[1,])
      g2hat_vect=c(g2hat_vect, mm_fit$param_new$gamma_new[2,])
      nu2hat_vect=c(nu2hat_vect, mm_fit$param_new$nu2_new)
      tau2hat_vect=c(tau2hat_vect, mm_fit$param_new$tau2_new)
      alphahat_vect=c(alphahat_vect, mm_fit$param_new$alpha_new)
      betahat_vect=c(betahat_vect, mm_fit$param_new$beta_new)
      
      # Store parameter estimates from linear mixed-effects model
      g1hat_vect_lme=c(g1hat_vect_lme, mm_fit_lme$coefficients$fixed[1])
      g2hat_vect_lme=c(g2hat_vect_lme, mm_fit_lme$coefficients$fixed[2])
      nu2hat_vect_lme=c(nu2hat_vect_lme, as.numeric(VarCorr(mm_fit_lme)[1,1]))
      tau2hat_vect_lme=c(tau2hat_vect_lme, as.numeric(VarCorr(mm_fit_lme)[2,1]))
      
    }		
  }
  
  # Return all collected estimates
  list(gamma1_est = g1hat_vect, gamma2_est=g2hat_vect, nusq_est=nu2hat_vect, tausq_est=tau2hat_vect, alpha_est=alphahat_vect,beta_est=betahat_vect, gamma1_est_lme=g1hat_vect_lme, gamma2_est_lme=g2hat_vect_lme, nusq_est_lme=nu2hat_vect_lme, tausq_est_lme=tau2hat_vect_lme)
  
}

# Function to calculate Mean Squared Error (MSE), Bias and Variance for parameter estimates
mse_bias <- function(n.base, n, M, by, nusq, alpha, beta, tausq, gamma1, gamma2){
  
  ss = simul_study(n.base, n, M, by, nusq, alpha, beta, tausq, gamma1, gamma2)
  
  # Calculate MSE, Bias and Variance for each parameter from EM algorithm
  g1_MSE = mean((ss$gamma1_est-gamma1)^2)
  g1_BIAS = mean(ss$gamma1_est)-gamma1
  g1_VAR = var(ss$gamma1_est)
  
  g2_MSE = mean((ss$gamma2_est-gamma2)^2)
  g2_BIAS = mean(ss$gamma2_est)-gamma2
  g2_VAR = var(ss$gamma2_est)
  
  nusq_MSE = mean((ss$nusq_est-nusq)^2)
  nusq_BIAS = mean(ss$nusq_est)-nusq
  nusq_VAR = var(ss$nusq_est)
  
  tausq_MSE = mean((ss$tausq_est-tausq)^2)
  tausq_BIAS = mean(ss$tausq_est)-tausq
  tausq_VAR = var(ss$tausq_est)
  
  alpha_MSE = mean((ss$alpha_est-alpha)^2)
  alpha_BIAS = mean(ss$alpha_est)-alpha
  alpha_VAR = var(ss$alpha_est)
  
  beta_MSE = mean((ss$beta_est-beta)^2)
  beta_BIAS = mean(ss$beta_est)-beta
  beta_VAR = var(ss$beta_est)
  
  # Calculate MSE, Bias and Variance for linear mixed-effects model estimates
  g1_lme_MSE = mean((ss$gamma1_est_lme-gamma1)^2)
  g1_lme_BIAS = mean(ss$gamma1_est_lme)-gamma1
  g1_lme_VAR = var(ss$gamma1_est_lme)
  
  g2_lme_MSE = mean((ss$gamma2_est_lme-gamma2)^2)
  g2_lme_BIAS = mean(ss$gamma2_est_lme)-gamma2
  g2_lme_VAR = var(ss$gamma2_est_lme)
  
  nusq_lme_MSE = mean((ss$nusq_est_lme-nusq)^2)
  nusq_lme_BIAS = mean(ss$nusq_est_lme)-nusq
  nusq_lme_VAR = var(ss$nusq_est_lme)
  
  tausq_lme_MSE = mean((ss$tausq_est_lme-tausq)^2)
  tausq_lme_BIAS = mean(ss$tausq_est_lme)-tausq
  tausq_lme_VAR = var(ss$tausq_est_lme)
  
  # Return all calculated metrics
  list(g1_MSE=g1_MSE, g1_BIAS=g1_BIAS,  g1_VAR = g1_VAR ,g2_MSE=g2_MSE, g2_BIAS=g2_BIAS, g2_VAR = g2_VAR , nusq_MSE=nusq_MSE, nusq_BIAS=nusq_BIAS, nusq_VAR = nusq_VAR ,  tausq_MSE=tausq_MSE, tausq_BIAS=tausq_BIAS,tausq_VAR = tausq_VAR , alpha_MSE=alpha_MSE, alpha_BIAS=alpha_BIAS,  alpha_VAR = alpha_VAR, beta_MSE=beta_MSE, beta_BIAS=beta_BIAS,  beta_VAR = beta_VAR, 
       g1_MSE_lme=g1_lme_MSE, g1_BIAS_lme=g1_lme_BIAS,  g1_VAR_lme = g1_lme_VAR ,g2_MSE_lme=g2_lme_MSE, g2_BIAS_lme=g2_lme_BIAS, g2_VAR_lme = g2_lme_VAR , nusq_MSE_lme=nusq_lme_MSE, nusq_BIAS_lme=nusq_lme_BIAS, nusq_VAR_lme = nusq_lme_VAR ,  tausq_MSE_lme=tausq_lme_MSE, tausq_BIAS_lme=tausq_lme_BIAS,tausq_VAR_lme = tausq_lme_VAR, SS=ss) # Also return all estimates for further analysis if needed
  
}

# Examples of running the mse_bias function with specified parameters
mm100 <- mse_bias(n.base=100, n=100, M=20, by=1, nusq=5, alpha=1, beta=1, tausq=10, gamma1=7, gamma2=0.5)

mm200 <- mse_bias(n.base=100, n=200, M=20, by=1, nusq=5, alpha=1, beta=1, tausq=10, gamma1=7, gamma2=0.5)
