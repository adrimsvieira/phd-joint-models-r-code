library(geoR)

####################################################################################
################################## Intensity 1 #####################################
####################################################################################

# function to simulate stationary Gaussian process W(s)
Wsim <- function(M,by,sigmasq,phi,k=0.5){
  s <- seq(0,M,by) # generate time grid from 0 to M
  g <- cbind(s,rep(0,length(s)))
  list(W=grf(M,grid=g,cov.pars=c(sigmasq,phi),kappa=k)$data,S=s)
}

# main simulation function
Ysim <- function(M,by,sigmasq,phi,k,alpha,beta,tausq){
  t <- NULL # initialize observation times
  y <- NULL # store time points
  tau <- sqrt(tausq) # standard deviation of measurement error
  
  W.ss < -Wsim(M,by,sigmasq,phi,k) # simulate latent Gaussian process
  W <- W.ss$W # extract simulated W(s)
  S <- W.ss$S # extract time points 
  
  # intensity function lambda(s)
  lambda.a <- function(alpha,beta,W){
    exp(alpha+beta*W)
  }
  
  # loop over time grid
  for (i in 1:length(S)) { 
    lambda.current <- lambda.a(alpha,beta,W[i]) # compute lambda at time S[i]
    p <- 1-exp(-lambda.current) 
    x = rbinom(1 , prob= p , size=1) # simulate event occurrence
    
    if(x==1){ # simulate event occurrence
      t <- c(t,S[i]) # store observation time
      y <- c(y,W[i]+tau*rnorm(1)) # store longitudinal process
    }
  }
  list(W=W,S=S,T=t,Y=y) # return simulated latent process, time grid, event times, and longitudinal data
}

######################################################
# result# gives the simulated data, uncomment to run #
######################################################

# to obtain the plots, simulate data with M=120 (121 equidistant time points), sigmasq=1, phi=0.25, and tausq=0
M<-120
by=1
sigmasq<-1
phi<-0.25
tausq<-0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
#result1<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result1$S,result1$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
#result2<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result2$S,result2$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
#result3<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result3$S,result3$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
#result4<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result4$S,result4$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
#result5<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result5$S,result5$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
#result6<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result6$S,result6$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
#result7<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result7$S,result7$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
#result8<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result8$S,result8$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
#result9<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result9$S,result9$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
#result10<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result10$S,result10$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
#result11<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result11$S,result11$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
#result12<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result12$S,result12$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
#result13<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result13$S,result13$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
#result14<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result14$S,result14$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
#result15<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result15$S,result15$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
#result16<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result16$S,result16$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
#result17<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result17$S,result17$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
#result18<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result18$S,result18$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
#result19 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result19$S,result19$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
#result20 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result20$S,result20$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
#result21 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result21$S,result21$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
#result22 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result22$S,result22$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
#result23 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result23$S,result23$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
#result24 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result24$S,result24$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
#result25 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result25$S,result25$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
#result26 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result26$S,result26$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
#result27 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result27$S,result27$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
#result28 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result28$S,result28$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
#result29 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result29$S,result29$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
#result30 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result30$S,result30$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
#result31 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result31$S,result31$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
#result32 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result32$S,result32$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
#result33 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result33$S,result33$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
#result34 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result34$S,result34$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
#result35 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result35$S,result35$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
#result36 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result36$S,result36$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=120 (121 equidistant time points), sigmasq=20, phi=7, and tausq=0
M<-120
by=1
sigmasq<-20
phi<-7
tausq<-0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
#result1<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result1$S,result1$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
#result2<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result2$S,result2$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
#result3<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result3$S,result3$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
#result4<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result4$S,result4$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
#result5<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result5$S,result5$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
#result6<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result6$S,result6$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
#result7<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result7$S,result7$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
#result8<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result8$S,result8$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
#result9<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result9$S,result9$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
#result10<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result10$S,result10$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
#result11<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result11$S,result11$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
#result12<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result12$S,result12$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
#result13<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result13$S,result13$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
#result14<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result14$S,result14$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
#result15<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result15$S,result15$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
#result16<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result16$S,result16$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
#result17<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result17$S,result17$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
#result18<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result18$S,result18$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
#result19 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result19$S,result19$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
#result20 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result20$S,result20$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
#result21 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result21$S,result21$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
#result22 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result22$S,result22$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
#result23 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result23$S,result23$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
#result24 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result24$S,result24$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
#result25 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result25$S,result25$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
#result26 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result26$S,result26$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
#result27 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result27$S,result27$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
#result28 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result28$S,result28$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
#result29 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result29$S,result29$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
#result30 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result30$S,result30$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
#result31 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result31$S,result31$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
#result32 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result32$S,result32$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
#result33 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result33$S,result33$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
#result34 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result34$S,result34$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
#result35 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result35$S,result35$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
#result36 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result36$S,result36$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=320 (321 equidistant time points), sigmasq=1, phi=0.25, and tausq=0
M<-320
by=1
sigmasq<-1
phi<-0.25
tausq<-0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
#result1<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result1$S,result1$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
#result2<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result2$S,result2$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
#result3<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result3$S,result3$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
#result4<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result4$S,result4$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
#result5<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result5$S,result5$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
#result6<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result6$S,result6$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
#result7<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result7$S,result7$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
#result8<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result8$S,result8$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
#result9<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result9$S,result9$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
#result10<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result10$S,result10$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
#result11<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result11$S,result11$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
#result12<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result12$S,result12$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
#result13<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result13$S,result13$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
#result14<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result14$S,result14$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
#result15<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result15$S,result15$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
#result16<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result16$S,result16$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
#result17<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result17$S,result17$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
#result18<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result18$S,result18$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
#result19 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result19$S,result19$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
#result20 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result20$S,result20$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
#result21 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result21$S,result21$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
#result22 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result22$S,result22$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
#result23 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result23$S,result23$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
#result24 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result24$S,result24$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
#result25 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result25$S,result25$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
#result26 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result26$S,result26$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
#result27 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result27$S,result27$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
#result28 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result28$S,result28$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
#result29 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result29$S,result29$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
#result30 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result30$S,result30$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
#result31 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result31$S,result31$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
#result32 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result32$S,result32$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
#result33 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result33$S,result33$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
#result34 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result34$S,result34$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
#result35 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result35$S,result35$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
#result36 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result36$S,result36$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=320 (321 equidistant time points), sigmasq=20, phi=7, and tausq=0
M<-320
by=1
sigmasq<-20
phi<-7
tausq<-0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
#result1<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result1$S,result1$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
#result2<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result2$S,result2$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
#result3<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result3$S,result3$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
#result4<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result4$S,result4$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
#result5<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result5$S,result5$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
#result6<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result6$S,result6$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
#result7<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result7$S,result7$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
#result8<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result8$S,result8$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
#result9<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result9$S,result9$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
#result10<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result10$S,result10$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
#result11<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result11$S,result11$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
#result12<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result12$S,result12$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
#result13<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result13$S,result13$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
#result14<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result14$S,result14$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
#result15<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result15$S,result15$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
#result16<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result16$S,result16$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
#result17<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result17$S,result17$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
#result18<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result18$S,result18$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
#result19 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result19$S,result19$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
#result20 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result20$S,result20$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
#result21 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result21$S,result21$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
#result22 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result22$S,result22$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
#result23 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result23$S,result23$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
#result24 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result24$S,result24$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
#result25 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result25$S,result25$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
#result26 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result26$S,result26$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
#result27 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result27$S,result27$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
#result28 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result28$S,result28$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
#result29 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result29$S,result29$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
#result30 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result30$S,result30$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
#result31 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result31$S,result31$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
#result32 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result32$S,result32$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
#result33 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result33$S,result33$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
#result34 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result34$S,result34$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
#result35 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result35$S,result35$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
#result36 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result36$S,result36$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=520 (521 equidistant time points), sigmasq=1, phi=0.25, and tausq=0
M<-520
by=1
sigmasq<-1
phi<-0.25
tausq<-0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
#result1<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result1$S,result1$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
#result2<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result2$S,result2$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
#result3<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result3$S,result3$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
#result4<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result4$S,result4$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
#result5<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result5$S,result5$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
#result6<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result6$S,result6$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
#result7<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result7$S,result7$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
#result8<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result8$S,result8$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
#result9<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result9$S,result9$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
#result10<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result10$S,result10$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
#result11<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result11$S,result11$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
#result12<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result12$S,result12$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
#result13<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result13$S,result13$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
#result14<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result14$S,result14$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
#result15<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result15$S,result15$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
#result16<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result16$S,result16$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
#result17<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result17$S,result17$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
#result18<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result18$S,result18$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
#result19 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result19$S,result19$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
#result20 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result20$S,result20$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
#result21 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result21$S,result21$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
#result22 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result22$S,result22$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
#result23 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result23$S,result23$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
#result24 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result24$S,result24$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
#result25 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result25$S,result25$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
#result26 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result26$S,result26$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
#result27 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result27$S,result27$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
#result28 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result28$S,result28$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
#result29 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result29$S,result29$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
#result30 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result30$S,result30$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
#result31 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result31$S,result31$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
#result32 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result32$S,result32$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
#result33 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result33$S,result33$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
#result34 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result34$S,result34$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
#result35 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result35$S,result35$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
#result36 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result36$S,result36$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)

# to obtain the plots, simulate data with M=520 (521 equidistant time points), sigmasq=20, phi=7, and tausq=0
M<-520
by=1
sigmasq<-20
phi<-7
tausq<-0
k=0.5

par(mfrow=c(6,3))
alpha<--10
beta<--10
#result1<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result1$S,result1$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result1$Y))) )
points(result1$T,result1$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result1$Y)

alpha<--10
beta<--5
#result2<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result2$S,result2$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result2$Y))) )
points(result2$T,result2$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result2$Y)

alpha<--10
beta<--2
#result3<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result3$S,result3$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result3$Y))))
points(result3$T,result3$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result3$Y)

alpha<--10
beta<-2
#result4<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result4$S,result4$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result4$Y))))
points(result4$T,result4$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result4$Y)

alpha<--10
beta<-5
#result5<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result5$S,result5$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result5$Y))))
points(result5$T,result5$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result5$Y)

alpha<--10
beta<-10
#result6<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result6$S,result6$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result6$Y))))
points(result6$T,result6$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result6$Y)

alpha<- -5
beta<- -10
#result7<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result7$S,result7$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result7$Y))))
points(result7$T,result7$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result7$Y)

alpha<- -5
beta<- -5
#result8<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result8$S,result8$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result8$Y))))
points(result8$T,result8$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result8$Y)

alpha<- -5
beta<- -2
#result9<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result9$S,result9$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result9$Y))))
points(result9$T,result9$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result9$Y)

alpha<- -5
beta<- 2
#result10<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result10$S,result10$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result10$Y))))
points(result10$T,result10$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result10$Y)

alpha<- -5
beta<- 5
#result11<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result11$S,result11$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result11$Y))))
points(result11$T,result11$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result11$Y)

alpha<- -5
beta<- 10
#result12<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result12$S,result12$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result12$Y))))
points(result12$T,result12$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result12$Y)

alpha<- -2
beta<- -10
#result13<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result13$S,result13$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result13$Y))))
points(result13$T,result13$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result13$Y)

alpha<- -2
beta<- -5
#result14<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result14$S,result14$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result14$Y))))
points(result14$T,result14$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result14$Y)

alpha<- -2
beta<- -2
#result15<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result15$S,result15$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result15$Y))))
points(result15$T,result15$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result15$Y)

alpha<- -2
beta<- 2
#result16<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result16$S,result16$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result16$Y))))
points(result16$T,result16$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result16$Y)

alpha<- -2
beta<- 5
#result17<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result17$S,result17$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result17$Y))))
points(result17$T,result17$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result17$Y)

alpha<- -2
beta<- 10
#result18<-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result18$S,result18$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result18$Y))))
points(result18$T,result18$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result18$Y)

par(mfrow=c(6,3))
alpha<- 2
beta<- -10
#result19 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result19$S,result19$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result19$Y))))
points(result19$T,result19$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result19$Y)

alpha<- 2
beta<- -5
#result20 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result20$S,result20$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result20$Y))))
points(result20$T,result20$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result20$Y)

alpha<- 2
beta<- -2
#result21 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result21$S,result21$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result21$Y))))
points(result21$T,result21$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result21$Y)

alpha<- 2
beta<- 2
#result22 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result22$S,result22$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result22$Y))))
points(result22$T,result22$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result22$Y)

alpha<- 2
beta<- 5
#result23 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result23$S,result23$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result23$Y))))
points(result23$T,result23$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result23$Y)

alpha<- 2
beta<- 10
#result24 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result24$S,result24$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result24$Y))))
points(result24$T,result24$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result24$Y)

alpha<- 5
beta<- -10
#result25 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result25$S,result25$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result25$Y))))
points(result25$T,result25$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result25$Y)

alpha<- 5
beta<- -5
#result26 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result26$S,result26$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result26$Y))))
points(result26$T,result26$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result26$Y)

alpha<- 5
beta<- -2
#result27 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result27$S,result27$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result27$Y))))
points(result27$T,result27$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result27$Y)

alpha<- 5
beta<- 2
#result28 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result28$S,result28$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result28$Y))))
points(result28$T,result28$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result28$Y)

alpha<- 5
beta<- 5
#result29 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result29$S,result29$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result29$Y))))
points(result29$T,result29$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result29$Y)

alpha<- 5
beta<- 10
#result30 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result30$S,result30$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result30$Y))))
points(result30$T,result30$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result30$Y)

alpha<- 10
beta<- -10
#result31 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result31$S,result31$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result31$Y))))
points(result31$T,result31$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result31$Y)

alpha<- 10
beta<- -5
#result32 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result32$S,result32$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result32$Y))))
points(result32$T,result32$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result32$Y)

alpha<- 10
beta<- -2
#result33 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result33$S,result33$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result33$Y))))
points(result33$T,result33$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result33$Y)

alpha<- 10
beta<- 2
#result34 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result34$S,result34$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result34$Y))))
points(result34$T,result34$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result34$Y)

alpha<- 10
beta<- 5
#result35 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result35$S,result35$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result35$Y))))
points(result35$T,result35$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result35$Y)

alpha<- 10
beta<- 10
#result36 <-Ysim(M,by,sigmasq,phi,k,alpha,beta,tausq)
plot(result36$S,result36$W,type="l",xlab="Time",ylab="W(s)", main=bquote(~ alpha == .(alpha) ~ "," ~ beta == .(beta) ~ "," ~ k == .(length(result36$Y))))
points(result36$T,result36$Y,pch=20,col="deepskyblue3",lwd=1)
abline(h=0, col="red", lty=2)
length(result36$Y)