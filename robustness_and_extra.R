## autocovar

Covar.curl.sm<-function(s,t,A){
  G<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  
  for (i in 1:401){
    A.beta.curl.sm.s<-t(A[[i]][s,])
    A.beta.curl.sm.t<-t(A[[i]][t,])
    G.i<-t(A.beta.curl.sm.s)%*%A.beta.curl.sm.t
    G<-G+G.i
  }
  return(G)
}


autocov<-function(delta, fnc,A){
  autovar<-matrix(0,nrow=168,ncol=length(beta[[1]]))
  for (t in 1:(168.9-delta)){
    autovar[t,]<-diag(fnc(t,t+delta,A))
  }
  return(colMeans(autovar))
}

delta<-seq(0.5,15, by=0.5)
sigma.delta<-matrix(0,nrow=length(delta),ncol=length(beta[[1]]))
for (i in 1:length(delta)){
  sigma.delta[i,]<-autocov(delta[i],Covar.curl.sm,A.beta.curl)
}
plot(delta,sigma.delta[,2],type='l',main='covariance for covariates',xlab='delta',ylab='mean_t Cov(t,t+delta)',ylim=c(0,0.08))
lines(delta,sigma.delta[,3])
lines(delta,sigma.delta[,4])
lines(delta,sigma.delta[,6])



## contours

t1<-seq(1,168)
t2<-seq(1,168)
cov.ts<-array(0,dim=c(length(t1),length(t2),length(beta[[1]]),length(beta[[1]])))
for (i in 166:length(t1)){
  for (j in 1:length(t2)){
    cov.ts[i,j,,]<-Covar.curl.sm(t1[i],t2[j],A.beta.curl.sm)
  }
}

# plot

require(grDevices) 
filled.contour(t1,t2,cov.ts[,,1,1])

# Plot
persp(t1, t2, cov.ts[,,7,7], main = "Contours of covariance for ISHAK", 
      col="red", theta = 55, phi = 30, r = 40, d = 0.1, zlab='Cov')




## robustness test


############################ Analysis as mean of effects ###############################



xmt<-demo[,c(1:7,9:16)]
beta.rob<-vector("list",24)
vhcid_list<-peg_rib_demo$vhcid[peg_rib_demo$week==3 & peg_rib_demo$n==1 & peg_rib_demo$n2<1 & peg_rib_demo$svr==1]
vhcid_list1<-peg_rib_demo$vhcid[peg_rib_demo$week==3 & peg_rib_demo$n==0 & peg_rib_demo$n2==1 & peg_rib_demo$svr==1]


for (i in 1:24){
  xvmt<-rib_peg[rib_peg$week==i,]
  xmt1<-merge(xvmt,xmt,by="vhcid")
  if (i==3) {
    #xmt1$n2[xmt1$vhcid==vhcid_list1[1]]<-5/7 
    xmt1$n[xmt1$vhcid==vhcid_list[1]]<-0 
  }## switching adherence of a single patient
  beta.rob[[i]]<-glm(svr~n*n2+SEX+RACEW+ISHAK, data=xmt1,family = binomial(link = "logit"),na.action="na.omit")$coefficients
}

time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)
yy<-rep(0,24)

for (i in 1:24){
  y1[i]<-beta.rob[[i]][3]
  y2[i]<-beta.rob[[i]][2]
  y3[i]<-beta.rob[[i]][7]
}

plot(time*7, y1, type='l',lwd=2,lty=1,ylim=c(-8,8),axes='F',ylab="change in log-odds for svr",xlab="Days",main="Adherence to Ribavirin/Peginterferon")
par(new=T)
plot(time*7, y2, type='l',lwd=2,lty=3,ylim=c(-8,8),axes='F',ylab="change in log-odds for svr",xlab='Days')
par(new=T)
plot(time*7, y3, type='l',lty=5,lwd=2,ylim=c(-8,8),axes='F',ylab="change in log-odds for svr",xlab='Days')
legend(15,40 ,lty=c(1,3,5),lwd=c(2,2,2), c('Ribavirin','Peginterferon','Interaction'), text.col = c(1,2,4) )
Axis(side=2)
Axis(side=1,at=c(7,42,84,126,168))
lines(time*7,rep(0,24))





############## final model #################


xmt<-demo[,c(1:7,9:16)]
beta.f<-vector("list",168)
vhcid_list<-peg_rib_demo$vhcid[peg_rib_demo$week==3 & peg_rib_demo$n==1 & peg_rib_demo$n2<1 & peg_rib_demo$svr==1]
vhcid_list1<-peg_rib_demo$vhcid[peg_rib_demo$week==3 & peg_rib_demo$n==0 & peg_rib_demo$n2==1 & peg_rib_demo$svr==1]



for (i in 1:168){
  xvmt<-rib[rib$Day==i,]
  xmt1<-merge(xvmt,xmt,by="vhcid")
  beta.f[[i]]<-glm(svr~n2+SEX+RACEW+ISHAK, data=xmt1,family = binomial(link = "logit"),na.action="na.omit")$coefficients
}

beta.ff<-matrix(unlist(beta.f),ncol=length(beta.f[[1]]))

