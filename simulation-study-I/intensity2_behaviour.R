library(geoR)


####################################################################################
################################## Intensity 2 #####################################
####################################################################################

# function to simulate stationary Gaussian process W(s)
Wsim<-function(M,by,sigmasq,phi,k=0.5){
  s<-seq(0,M,by) # generate time grid from 0 to M
  g<-cbind(s,rep(0,length(s)))
  list(W=grf(M,grid=g,cov.pars=c(sigmasq,phi),kappa=k)$data,S=s)
}

# main simulation function
Ysim<-function(W,S,alpha,beta,tausq,gamma1,gamma2){
  
  t<-NULL # initialize observation times
  y<-NULL # store time points
  tau<-sqrt(tausq) # standard deviation of measurement error
  
  # intensity function for Poisson process
  lambda<-function(alpha,beta,ww){
    exp(alpha+beta*sum(ww))
  }
  
  for (i in 1:length(S)) {
    
    if(i>=5){
      w<-NULL # generate event occurrence
      for(j in (i-4):i){
        w<-c(w,dnorm(S[j],mean=S[i],sd=sd(S[(i-4):i])))
      }
      
      w<-w/sum(w)
      
      # compute intensity using weighted sum of past W
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
    x = rbinom(1 , prob= p , size=1) # generate event occurrence
    if(x==1){
      t<-c(t,S[i]) # store observation time
      # generate corresponding Y value with measurement error and fixed effects
      y=c(y,gamma1+gamma2*S[i]+W[i]+rnorm(1,mean=0 , sd = sqrt(tausq)))
    }
  }
  
  list(T=t,Y=y) # compute intensity using weighted sum of past W
}

# to obtain the plots, simulate data with M=120 (121 equidistant time points), sigmasq=1, phi=0.25, tausq=0, gamma0=0, gamma1=0
M<-120
by=1
sigmasq<-1
phi<-0.25
tausq<-0
gamma1 <- 0
gamma2 <- 0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result1<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result2<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result3<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result4<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result5<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result6<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result7<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result8<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result9<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result10<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result11<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result12<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result13<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result14<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result15<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result16<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result17<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result18<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result19<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result20<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result21<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result22<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result23<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result24<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result25<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result26<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result27<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result28<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result29<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result30<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result31<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result32<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result33<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result34<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result35<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result36<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=120 (121 equidistant time points), sigmasq=20, phi=7, tausq=0, gamma0=0, gamma1=0
M<-120
by=1
sigmasq<-20
phi<-7
tausq<-0
gamma1 <- 0
gamma2 <- 0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result1<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result2<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result3<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result4<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result5<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result6<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result7<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result8<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result9<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result10<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result11<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result12<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result13<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result14<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result15<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result16<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result17<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result18<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result19<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result20<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result21<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result22<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result23<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result24<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result25<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result26<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result27<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result28<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result29<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result30<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result31<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result32<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result33<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result34<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result35<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result36<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=320 (321 equidistant time points), sigmasq=1, phi=0.25, tausq=0, gamma0=0, gamma1=0
M<-320
by=1
sigmasq<-1
phi<-0.25
tausq<-0
gamma1 <- 0
gamma2 <- 0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result1<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result2<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result3<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result4<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result5<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result6<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result7<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result8<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result9<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result10<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result11<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result12<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result13<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result14<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result15<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result16<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result17<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result18<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result19<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result20<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result21<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result22<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result23<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result24<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result25<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result26<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result27<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result28<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result29<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result30<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result31<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result32<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result33<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result34<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result35<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result36<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=320 (321 equidistant time points), sigmasq=20, phi=7, tausq=0, gamma0=0, gamma1=0
M<-320
by=1
sigmasq<-20
phi<-7
tausq<-0
gamma1 <- 0
gamma2 <- 0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result1<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result2<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result3<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result4<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result5<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result6<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result7<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result8<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result9<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result10<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result11<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result12<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result13<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result14<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result15<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result16<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result17<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result18<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result19<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result20<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result21<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result22<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result23<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result24<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result25<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result26<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result27<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result28<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result29<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result30<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result31<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result32<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result33<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result34<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result35<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result36<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=520 (521 equidistant time points), sigmasq=1, phi=0.25, tausq=0, gamma0=0, gamma1=0
M<-520
by=1
sigmasq<-1
phi<-0.25
tausq<-0
gamma1 <- 0
gamma2 <- 0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result1<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result2<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result3<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result4<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result5<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result6<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result7<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result8<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result9<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result10<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result11<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result12<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result13<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result14<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result15<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result16<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result17<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result18<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result19<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result20<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result21<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result22<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result23<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result24<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result25<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result26<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result27<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result28<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result29<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result30<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result31<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result32<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result33<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result34<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result35<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result36<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=520 (521 equidistant time points), sigmasq=20, phi=7, tausq=0, gamma0=0, gamma1=0
M<-520
by=1
sigmasq<-20
phi<-7
tausq<-0
gamma1 <- 0
gamma2 <- 0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result1<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result2<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result3<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result4<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result5<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result6<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result7<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result8<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result9<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result10<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result11<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result12<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result13<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result14<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result15<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result16<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result17<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result18<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result19<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result20<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result21<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result22<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result23<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result24<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result25<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result26<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result27<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result28<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result29<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result30<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result31<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result32<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result33<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result34<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result35<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
WW <- Wsim(M,by,sigmasq,phi,k=0.5)
result36<-Ysim(WW$W,WW$S,alpha,beta,tausq,gamma1,gamma2)
plot(WW$S,WW$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)