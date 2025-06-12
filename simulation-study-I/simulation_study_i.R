library(geoR)
library(nlme)

# Function to simulate W(s) using exponential covariance structure
W.proc<-function(M,by,sigmasq,phi){
  s<-seq(0,M,by)
  g<-cbind(s,rep(0,length(s)))
  
  # simulate Gaussian random field with exponential covariance
  list(W=grf(M,grid=g,cov.model="exponential",cov.pars=c(sigmasq,phi),messages=FALSE)$data,S=s)
}

# Wrapper to extract W and S directly
Wsim<-function(M,by,sigmasq,phi){
  W.ss<-W.proc(M,by,sigmasq,phi)
  W<-W.ss$W
  S<-W.ss$S
  
  list(W=W,S=S)
}

#####################################################################################################
########################################### Intensity 1 #############################################
#####################################################################################################

# Simulation with intensity scenario 1 and a random effect U 
Ysim1u<-function(W,S,alpha,beta,tausq,gamma1,gamma2,nusq){
  t<-NULL
  y<-NULL
  tau<-sqrt(tausq)
  
  lambda<-function(alpha,beta,W){
    exp(alpha+beta*W)
  }
  
  for (i in 1:length(S)) {
    lambda.current<-lambda(alpha,beta,W[i])
    p<-1-exp(-lambda.current)
    x = rbinom(1 , prob= p , size=1)
    if(x==1){
      t<-c(t,S[i])
      u = rnorm(1,mean = 0, sd=sqrt(nusq)) # subject-specific random effect
      y=c(y,gamma1+gamma2*S[i]+u+rnorm(1,mean=0 , sd = sqrt(tausq)))
    }
  }
  list(T=t,Y=y)
}

# Simulation with intensity scenario 1 (without U)
Ysim1<-function(W,S,alpha,beta,tausq,gamma1,gamma2){
  t<-NULL
  y<-NULL
  tau<-sqrt(tausq)
  
  lambda<-function(alpha,beta,W){
    exp(alpha+beta*W)
  }
  
  for (i in 1:length(S)) {
    lambda.current<-lambda(alpha,beta,W[i])
    p<-1-exp(-lambda.current)
    x = rbinom(1 , prob= p , size=1)
    if(x==1){
      t<-c(t,S[i])
      y=c(y,gamma1+gamma2*S[i]+W[i]+rnorm(1,mean=0 , sd = sqrt(tausq)))
    }
  }
  list(T=t,Y=y)
}

# Simulation with intensity scenario 1 (without U)
sim_data_set1<- function(n , M , by , sigmasq , phi , alpha, beta, nusq , tausq , gamma1 , gamma2, method=c("standard","random")){
  
  dd = data.frame(id=NULL , time = NULL , Y = NULL)
  
  for(l in 1:n){
    if(method=="standard"){
      
      W.true = Wsim(M,by,sigmasq,phi)
      W = W.true$W
      S = W.true$S
      
      Y.l = Ysim1(W, S, alpha, beta , tausq , gamma1, gamma2)
    }else{
      
      W.true = Wsim(M,by,sigmasq,phi)
      W = W.true$W
      S = W.true$S
      
      Y.l = Ysim1u(W, S, alpha, beta , tausq , nusq, gamma1, gamma2)
    }
    
    id = rep(l , length(Y.l$Y))
    
    dd.l = cbind(id , time=Y.l$T , Y=Y.l$Y)
    
    dd = rbind(dd , dd.l)
  }
  
  list(data =dd )
}

# Model fitting function depending on simulation method
model_fit1 <- function(dados, method=c("standard","random")){
  
  if(method=="standard"){
    
    m <- gls(Y~time,correlation=corExp(form=~ time|id ,nugget=TRUE),data=dados$data)
    
    g1_hat = as.numeric(coefficients(m))[1]
    g2_hat = as.numeric(coefficients(m))[2]			
    phi_hat = as.numeric(coef(m$modelStruct$corStruct,unconstrained=FALSE)[1])
    rstderror = as.numeric(sigma(m))
    nugget = as.numeric(coef(m$modelStruct$corStruct,unconstrained=FALSE)[2])
    tau2_hat = (rstderror^2)*nugget
    sigma2_hat= (rstderror^2)*(1-nugget)
    nu2_hat=NULL
    
  }else{
    
    
    m <- lme(Y~time,random=~1|id, correlation=corExp(form=~ time|id ,nugget=TRUE), data=dados$data)
    
    g1_hat = as.numeric(m$coefficients$fixed[1])
    g2_hat = as.numeric(m$coefficients$fixed[2])
    phi_hat = as.numeric(coef(m$modelStruct$corStruct,unconstrained=FALSE)[1])
    nu2_hat = as.numeric(VarCorr(m)[1,1])
    rstderror = as.numeric(VarCorr(m)[2,2])
    nugget = as.numeric(coef(m$modelStruct$corStruct,unconstrained=FALSE)[2])
    tau2_hat = (rstderror^2)*nugget
    sigma2_hat= (rstderror^2)*(1-nugget)
    
  }
  
  list(gamma1_hat=g1_hat,gamma2_hat=g2_hat,phi_hat=phi_hat,nusq_hat=nu2_hat,sigmasq_hat=sigma2_hat,tausq_hat=tau2_hat)
}

# Simulation study function
simul_study1 = function(n.base, n, M, by, sigmasq, phi, alpha, beta, nusq, tausq, gamma1, gamma2, method.sim.data.set=c("standard","random"), method.model.fit=c("standard","random")){
  
  g1hat_vect=NULL
  g2hat_vect=NULL
  phihat_vect=NULL
  nu2hat_vect=NULL
  sigma2hat_vect=NULL
  tau2hat_vect=NULL
  
  ii = 1
  while(ii <=n.base){
    
    print(paste0("n.base ", ii)) # progress message
    
    dd = sim_data_set1(n, M, by, sigmasq, phi, alpha, beta, nusq, tausq, gamma1, gamma2, method=method.sim.data.set)	
    
    mm_fit = tryCatch({
      
      model_fit1(dd, method=method.model.fit)
      
    } , error=function(e){print(paste0("ERROR na base", ii ,conditionMessage(e),"\n")) ; 0 }
    )
    
    if(is.numeric(mm_fit)){  }else{
      ii = ii + 1
      g1hat_vect=c(g1hat_vect, mm_fit$gamma1_hat)
      g2hat_vect=c(g2hat_vect, mm_fit$gamma2_hat)
      phihat_vect=c(phihat_vect, mm_fit$phi_hat)
      nu2hat_vect=c(nu2hat_vect, mm_fit$nusq_hat)
      sigma2hat_vect=c(sigma2hat_vect, mm_fit$sigmasq_hat)
      tau2hat_vect=c(tau2hat_vect, mm_fit$tausq_hat)
    }		
  }
  
  list(gamma1_est = g1hat_vect, gamma2_est=g2hat_vect, phi_est=phihat_vect, nusq_est=nu2hat_vect, sigmasq_est=sigma2hat_vect, tausq_est=tau2hat_vect)
  
}

# Function to compute MSE, Bias and Variance across parameter combinations
mse_bias1 <- function(n.base, n, M, by, sigmasq, phi, alpha, beta, nusq, tausq, gamma1, gamma2, method.sim.data.set, method.model.fit){
  
  g1_MSE <- matrix(ncol=length(alpha) , nrow =length(beta))
  g1_BIAS <- matrix(ncol=length(alpha) , nrow =length(beta))
  g1_VAR <- matrix(ncol=length(alpha) , nrow =length(beta))
  g2_MSE <- matrix(ncol=length(alpha) , nrow =length(beta))
  g2_BIAS <- matrix(ncol=length(alpha) , nrow =length(beta))
  g2_VAR <- matrix(ncol=length(alpha) , nrow =length(beta))
  phi_MSE <- matrix(ncol=length(alpha) , nrow =length(beta))
  phi_BIAS <- matrix(ncol=length(alpha) , nrow =length(beta))
  phi_VAR <- matrix(ncol=length(alpha) , nrow =length(beta))
  nusq_MSE <- matrix(ncol=length(alpha) , nrow =length(beta))
  nusq_BIAS <- matrix(ncol=length(alpha) , nrow =length(beta))
  nusq_VAR <- matrix(ncol=length(alpha) , nrow =length(beta))
  sigmasq_MSE <- matrix(ncol=length(alpha) , nrow =length(beta))
  sigmasq_BIAS <- matrix(ncol=length(alpha) , nrow =length(beta))
  sigmasq_VAR <- matrix(ncol=length(alpha) , nrow =length(beta))
  tausq_MSE <- matrix(ncol=length(alpha) , nrow =length(beta))
  tausq_BIAS <- matrix(ncol=length(alpha) , nrow =length(beta))
  tausq_VAR <- matrix(ncol=length(alpha) , nrow =length(beta))
  
  ss_all = vector("list" , length(alpha)*length(beta))
  
  # simulation loops over all alpha and beta combinations
  for(i in 1:length(alpha)){
    print(paste0("alpha ", i))
    aa = alpha[i]
    
    for(j in 1:length(beta)){
      print(paste0("beta ", j))
      bb = beta[j]
      
      ss = simul_study1(n.base, n, M, by, sigmasq, phi, alpha=aa, beta=bb, nusq, tausq, gamma1, gamma2, method.sim.data.set, method.model.fit)
      
      ss_all[[length(beta)*(i-1) + j]] = ss
      
      g1_MSE[j,i] = mean((ss$gamma1_est-gamma1)^2)
      g1_BIAS[j,i] = mean(ss$gamma1_est)-gamma1
      g1_VAR[j,i] = var(ss$gamma1_est)
      
      g2_MSE[j,i] = mean((ss$gamma2_est-gamma2)^2)
      g2_BIAS[j,i] = mean(ss$gamma2_est)-gamma2
      g2_VAR[j,i] = var(ss$gamma2_est)
      
      phi_MSE[j,i] = mean((ss$phi_est-phi)^2)
      phi_BIAS[j,i] = mean(ss$phi_est)-phi
      phi_VAR[j,i] = var(ss$phi_est)
      
      if(method.model.fit=="standard"){
        nusq_MSE[j,i] = NA
        nusq_BIAS[j,i] = NA
        nusq_VAR[j,i] = NA
      }else{
        nusq_MSE[j,i] = mean((ss$nusq_est-nusq)^2)
        nusq_BIAS[j,i] = mean(ss$nusq_est)-nusq
        nusq_VAR[j,i] = var(ss$nusq_est)
      }
      
      sigmasq_MSE[j,i] = mean((ss$sigmasq_est-sigmasq)^2)
      sigmasq_BIAS[j,i] = mean(ss$sigmasq_est)-sigmasq
      sigmasq_VAR[j,i] = var(ss$sigmasq_est)
      
      tausq_MSE[j,i] = mean((ss$tausq_est-tausq)^2)
      tausq_BIAS[j,i] = mean(ss$tausq_est)-tausq
      tausq_VAR[j,i] = var(ss$tausq_est)
      
    }
  }
  
  list(g1_MSE=g1_MSE, g1_BIAS=g1_BIAS,  g1_VAR = g1_VAR ,g2_MSE=g2_MSE, g2_BIAS=g2_BIAS, g2_VAR = g2_VAR , phi_MSE=phi_MSE, phi_BIAS=phi_BIAS, phi_VAR = phi_VAR , nusq_MSE=nusq_MSE, nusq_BIAS=nusq_BIAS, nusq_VAR = nusq_VAR , sigmasq_MSE=sigmasq_MSE, sigmasq_BIAS=sigmasq_BIAS, sigmasq_VAR =sigmasq_VAR , tausq_MSE=tausq_MSE, tausq_BIAS=tausq_BIAS,tausq_VAR = tausq_VAR , SS=ss_all)
  
}

# 
# The mse_bias function may take some time to run.
# Therefore, we recommend running mse_bias1 once and saving the workspace image with:
#    # save.image("BiasStudy.Rdata")


stst <- mse_bias1(n.base=100, n=50, M=120, by=1, sigmasq=20, phi=7, alpha=-2:5, beta=-10:10, nusq=5, tausq=10, gamma1=7, gamma2=0.5, method.sim.data.set="standard", method.model.fit="standard")
#save.image("BiasStudy.Rdata")

library(fields)

par(mfrow=c(3,5),mar=c(7,4,4,1),oma=c(2,0,0,0))
image.plot(-2:5,-10:10,t(stst$g1_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(stst$g2_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(stst$phi_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ phi), line=1)
image.plot(-2:5,-10:10,t(stst$sigmasq_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(stst$tausq_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ tau^2), line=1)


image.plot(-2:5,-10:10,t(stst$g1_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(stst$g2_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(stst$phi_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ phi), line=1)
image.plot(-2:5,-10:10,t(stst$sigmasq_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(stst$tausq_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ tau^2), line=1)


image.plot(-2:5,-10:10,t(stst$g1_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(stst$g2_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(stst$phi_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=1.6)
title(main=expression("VAR" ~ phi), line=1)
image.plot(-2:5,-10:10,t(stst$sigmasq_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(stst$tausq_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ tau^2), line=1)

strd <- mse_bias1(n.base=100, n=50, M=120, by=1, sigmasq=20, phi=7, alpha=-2:5, beta=-10:10, nusq=5, tausq=10, gamma1=7, gamma2=0.5, method.sim.data.set="standard", method.model.fit="random")
#save.image("BiasStudy.Rdata")

par(mfrow=c(3,5),mar=c(7,4,4,1),oma=c(2,0,0,0))
image.plot(-2:5,-10:10,t(strd$g1_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(strd$g2_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(strd$phi_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ phi), line=1)
image.plot(-2:5,-10:10,t(strd$sigmasq_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(strd$tausq_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ tau^2), line=1)


image.plot(-2:5,-10:10,t(strd$g1_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(strd$g2_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(strd$phi_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ phi), line=1)
image.plot(-2:5,-10:10,t(strd$sigmasq_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(strd$tausq_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ tau^2), line=1)


image.plot(-2:5,-10:10,t(strd$g1_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(strd$g2_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(strd$phi_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=1.6)
title(main=expression("VAR" ~ phi), line=1)
image.plot(-2:5,-10:10,t(strd$sigmasq_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(strd$tausq_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ tau^2), line=1)

rdst <- mse_bias1(n.base=100, n=50, M=120, by=1, sigmasq=20, phi=7, alpha=-2:5, beta=-10:10, nusq=5, tausq=10, gamma1=7, gamma2=0.5, method.sim.data.set="random", method.model.fit="standard")
#save.image("BiasStudy.Rdata")

par(mfrow=c(3,5),mar=c(7,4,4,1),oma=c(2,0,0,0))
image.plot(-2:5,-10:10,t(rdst$g1_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(rdst$g2_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(rdst$phi_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ phi), line=1)
image.plot(-2:5,-10:10,t(rdst$sigmasq_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(rdst$tausq_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ tau^2), line=1)


image.plot(-2:5,-10:10,t(rdst$g1_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(rdst$g2_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(rdst$phi_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ phi), line=1)
image.plot(-2:5,-10:10,t(rdst$sigmasq_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(rdst$tausq_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ tau^2), line=1)


image.plot(-2:5,-10:10,t(rdst$g1_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(rdst$g2_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(rdst$phi_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ phi), line=1)
image.plot(-2:5,-10:10,t(rdst$sigmasq_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(rdst$tausq_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ tau^2), line=1)

rdrd <- mse_bias1(n.base=100, n=50, M=120, by=1, sigmasq=20, phi=7, alpha=-2:5, beta=-10:10, nusq=5, tausq=10, gamma1=7, gamma2=0.5, method.sim.data.set="random", method.model.fit="random")
#save.image("BiasStudy.Rdata")

par(mfrow=c(3,6),mar=c(7,4,4,1),oma=c(2,0,0,0))
image.plot(-2:5,-10:10,t(rdrd$g1_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(rdrd$g2_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(rdrd$nusq_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ nu^2), line=1)
image.plot(-2:5,-10:10,t(rdrd$phi_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ phi), line=1)
image.plot(-2:5,-10:10,t(rdrd$sigmasq_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(rdrd$tausq_MSE),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("MSE" ~ tau^2), line=1)


image.plot(-2:5,-10:10,t(rdrd$g1_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(rdrd$g2_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(rdrd$nusq_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ nu^2), line=1)
image.plot(-2:5,-10:10,t(rdrd$phi_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ phi), line=1)
image.plot(-2:5,-10:10,t(rdrd$sigmasq_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(rdrd$tausq_BIAS),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("BIAS" ~ tau^2), line=1)


image.plot(-2:5,-10:10,t(rdrd$g1_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ gamma[1]), line=1)
image.plot(-2:5,-10:10,t(rdrd$g2_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ gamma[2]), line=1)
image.plot(-2:5,-10:10,t(rdrd$phi_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=1.6)
title(main=expression("VAR" ~ phi), line=1)
image.plot(-2:5,-10:10,t(rdrd$nusq_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ nu^2), line=1)
image.plot(-2:5,-10:10,t(rdrd$sigmasq_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ sigma^2), line=1)
image.plot(-2:5,-10:10,t(rdrd$tausq_VAR),xlab="",ylab="",main="",horizontal=T)
title(xlab=expression(alpha),line=2)
title(ylab=expression(beta),line=2)
title(main=expression("VAR" ~ tau^2), line=1)

#####################################################################################################
########################################### Intensity 2 #############################################
#####################################################################################################

#
# Repeat for Intensity 2 #
#


#Simulation with intensity scenario 2, and W=U random effect
Ysim2u<-function(W,S,alpha,beta,tausq,gamma1,gamma2,nusq){
  
  t<-NULL
  y<-NULL
  tau<-sqrt(tausq)
  
  lambda<-function(alpha,beta,ww){
    exp(alpha+beta*sum(ww))
  }
  
  for (i in 1:length(S)) {
    
    if(i>=5){
      w<-NULL
      for(j in (i-4):i){
        w<-c(w,dnorm(S[j],mean=S[i],sd=sd(S[(i-4):i])))
      }
      
      w<-w/sum(w)
      
      lambda.current<-lambda(alpha,beta,(W[(i-4):i]*(w)))
      
    }else{
      w<-NULL
      for(j in 1:i){
        if(length(S[1:i])==1){SD = 1}else{SD=sd(S[1:i])}
        w<-c(w,dnorm(S[j],mean=S[i],sd=SD))
      }
      w<-w/sum(w)
      lambda.current<-lambda(alpha,beta,(W[1:i]*(w)))
      
    }
    
    
    p<-1-exp(-lambda.current)
    x = rbinom(1 , prob= p , size=1)
    if(x==1){
      t<-c(t,S[i])
      u = rnorm(1,mean = 0, sd=sqrt(nusq))
      y=c(y,gamma1+gamma2*S[i]+u+rnorm(1,mean=0 , sd = sqrt(tausq)))
    }
  }
  
  list(T=t,Y=y)
}


# Simulation with intensity scenario 2 (without U)
Ysim2<-function(W,S,alpha,beta,tausq,gamma1,gamma2){
  
  t<-NULL
  y<-NULL
  tau<-sqrt(tausq)
  
  lambda<-function(alpha,beta,ww){
    exp(alpha+beta*sum(ww))
  }
  
  for (i in 1:length(S)) {
    
    if(i>=5){
      w<-NULL
      for(j in (i-4):i){
        w<-c(w,dnorm(S[j],mean=S[i],sd=sd(S[(i-4):i])))
      }
      
      w<-w/sum(w)
      
      lambda.current<-lambda(alpha,beta,(W[(i-4):i]*(w)))
      
    }else{
      w<-NULL
      for(j in 1:i){
        if(length(S[1:i])==1){SD = 1}else{SD=sd(S[1:i])}
        w<-c(w,dnorm(S[j],mean=S[i],sd=SD))
      }
      w<-w/sum(w)
      lambda.current<-lambda(alpha,beta,(W[1:i]*(w)))
      
    }
    
    
    p<-1-exp(-lambda.current)
    x = rbinom(1 , prob= p , size=1)
    if(x==1){
      t<-c(t,S[i])
      #u = rnorm(1,mean = 0, sd=sqrt(nusq))
      y=c(y,gamma1+gamma2*S[i]+W[i]+rnorm(1,mean=0 , sd = sqrt(tausq)))
    }
  }
  
  list(T=t,Y=y)
}