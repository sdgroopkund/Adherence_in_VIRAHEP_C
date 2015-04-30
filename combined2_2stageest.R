homdir<-"C:/material/coursework old/coursework/coursework/tracs/data"
#homdir <- "H:/codes/survival/additive_model"
setwd(homdir)
source(file="fun_general.R")

################################################################
######## Data management #######################################
################################################################

demo<-read.table('demo.txt',header=T,sep='\t')
peg<-read.table('pegifn.txt',header=T,sep='\t')
rib<-read.table('rib.txt',header=T,sep='\t')
vload<-read.table('vload.txt',header=T,sep='\t')

demo[demo==""]<-NA
demo$ISHAK[demo$ISHAK=="B"]<-NA
demo$ISHAK <- as.numeric(as.character(demo$ISHAK))
write.table(demo, file = "demo1.txt",col.names = TRUE,sep='\t')
demo<-read.table('demo1.txt',header=T,sep='\t')

rib<-rib[-c(26881,72914), ]

rib<-cbind(rib,n1=(rib$n==1)+(rib$n==2),n2=as.numeric(rib$n==2))


## removing dropouts 

rib<-rib[rib$dropout==0,]
rib<-rib[rib$dropout2==0,]

#######################

L<-data.frame(cbind(vhcid=demo$vhcid,id=seq(1,length(demo$vhcid))))
rib_demo<-merge(rib,merge(L,demo,by='vhcid'),by='vhcid')


rib_demo1<-rib_demo[rib_demo$Day<=168,]
rib_demo21<-rib_demo[rib_demo$Day>168,]
rib_demo2<-rib_demo21[rib_demo21$response24=='Responder',]
rib_demo3<-rib_demo[rib_demo$response24=='Responder',]


peg<-peg[-8449,]


## removing dropouts #
peg2<-peg[is.na(peg$do_not_include)==1,]
peg3<-peg2[is.na(peg2$do_not_include2)==1,]

#######################




peg_demo<-merge(peg3,merge(L,demo,by='vhcid'),by='vhcid')


peg_demo1<-peg_demo[peg_demo$week<=24,]

peg_demo21<-peg_demo[peg_demo$week>24,]
peg_demo2<-peg_demo21[peg_demo21$response24=='Responder',]

peg_demo3<-peg_demo[peg_demo$response24=='Responder',]




############################ transposing dataset ################
rib_f24<-rib[rib$Day<=168,]
rib_cn<-cbind("vhcid"=rib_f24$vhcid,"Day"=rib_f24$Day,"n2"=rib_f24$n2)
rib_cn<-data.frame(rib_cn)
week_day<-rib_cn$Day%%7
week_day[which(week_day==0)]=7
rib_cn<-cbind(rib_cn,"wday"=week_day)
week<-(rib_cn$Day-1)%/%7+1
rib_cn<-cbind(rib_cn,"week"=week)

rib_r<-reshape(rib_cn,v.names="n2",timevar="wday",idvar = c("week","vhcid"),direction="wide",drop="Day")
peg_cn<-data.frame(cbind("vhcid"=peg3$vhcid,"week"=peg3$week,"n"=peg3$n))
peg_cn<-peg_cn[peg_cn$week<=24,]

rib_peg<-merge(peg_cn,rib_r,by=c("vhcid","week"))

rib_peg <- rib_peg[order(rib_peg$vhcid, rib_peg$week),]
rib_peg<-data.frame(cbind(rib_peg,"n2"=rowSums(rib_peg[,4:10])/7))

##########################################################################

L<-data.frame(cbind(vhcid=demo$vhcid,id=seq(1,length(demo$vhcid))))
peg_rib_demo<-merge(rib_peg,merge(L,demo,by='vhcid'),by='vhcid')
data.to.fit<-subset(rib_demo1,select = c(id,Day,svr,n2,SEX,RACEW,ISHAK))
data.to.fit<-cbind(data.to.fit[,1:3],offset=1,data.to.fit[,4:dim(data.to.fit)[2]])

data.to.fit2<-subset(peg_demo1,select = c(id,week,svr,n,SEX,RACEW,ISHAK))
data.to.fit2<-cbind(data.to.fit2[,1:3],offset=1,data.to.fit2[,4:dim(data.to.fit2)[2]])




###########################################
### Combined analysis 2 ###################
###########################################

peg.demo.new <- subset(peg_demo1, select =c(vhcid,id, week, n, svr, SEX, RACEW, ISHAK))
peg.demo.new2<-data.frame(cbind(peg.demo.new[rep(seq_len(nrow(peg.demo.new)), each=7),],weekday=rep(1:7,times=9113)))
Day<-(peg.demo.new2$week-1)*7+peg.demo.new2$weekday
peg.demo.new3<-data.frame(cbind(peg.demo.new2,Day=Day))

peg.rib.demo.new<-merge(rib_cn,peg.demo.new3,by=c("vhcid","Day"),all.x=T)

xmt<-demo[,c(1:7,9:16)]
#xmt<-model.matrix(~SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=demo)


for (i in 1:168){
  xvmt<-peg.rib.demo.new[peg.rib.demo.new$Day==i,]
  beta[[i]]<-glm(svr~n*n2+SEX+RACEW+ISHAK, data=xvmt,family = binomial(link = "logit"),na.action="na.omit")$coefficients
}


beta.m<-matrix(0,168,6)

for (i in 1:168){
  xvmt<-peg.rib.demo.new[peg.rib.demo.new$Day==i,]
  beta.m[i,]<-glm(svr~n+n2+SEX+RACEW+ISHAK, data=xvmt,family = binomial(link = "logit"),na.action="na.omit")$coefficients
}

time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  y1[i]<-beta[[i]][2]
  y2[i]<-beta[[i]][3]
  y3[i]<-beta[[i]][4]
  y4[i]<-beta[[i]][5]
  y5[i]<-beta[[i]][7]
}

plot(time, y1, type='l',ylim=c(-15,15),ylab="Compliance")
par(new=T)
plot(time, y2, type='l',ylim=c(-15,15),col="blue",ylab="Compliance")
par(new=T)
plot(time, y7, type='l',ylim=c(-15,15),col=2,ylab="Compliance")


plot(time, y1, type='l',lwd=2,ylim=c(-25,35),axes='F',ylab="change in log-odds for svr",xlab="Days",main="Adherence to Ribavirin/Peginterferon")
par(new=T)
plot(time, y2, type='l',lty=2,lwd=2,ylim=c(-25,35),axes='F',ylab="change in log-odds for svr",xlab='Days')
par(new=T)
plot(time, y5, type='l',lwd=2,lty=6,ylim=c(-25,35),axes='F',ylab="change in log-odds for svr",xlab='Days')
legend(100,35 ,lty=c(1,2,6),lwd=c(2,2,2), c('Ribavirin','Peginterferon','Interaction') )
Axis(side=2)
Axis(side=1,at=c(7,42,84,126,168))

X<-subset(peg.rib.demo.new,select = c(svr,id,Day,n,n2,SEX,RACEW,ISHAK))
H.t.arr<-array(0,c(ncol(beta.m),ncol(beta.m),168))
A.beta.curl.new<-array(0,c(401,6,168))
At<-array(0,c(401,6,168))
for (t in 1:168){
  ZZ.t<-model.matrix(~svr+id+n+n2+SEX+RACEW+ISHAK ,data=X[X$Day==t,])
  Z.t<-model.matrix(~n+n2+SEX+RACEW+ISHAK ,data=X[X$Day==t,])
  D.beta.t<-apply(Z.t, 2,function(x){x*exp(Z.t%*%beta.m[t,])/((1+exp(Z.t%*%beta.m[t,]))^2)})
  Y.bar.t<- ZZ.t[,2]-exp(Z.t%*%beta.m[t,])/((1+exp(Z.t%*%beta.m[t,])))  
  V.beta.t<-(1+exp(Z.t%*%beta.m[t,]))^2/exp(Z.t%*%beta.m[t,])
  mat<-apply(D.beta.t,1,function(z){z%*%t(z)})
  Ht<-apply(mat,1,function(x){x*V.beta.t})
  H.t.arr[,,t]<-matrix(colSums(Ht),ncol(beta.m),ncol(beta.m))
  vec1<-V.beta.t*Y.bar.t
  At[ZZ.t[,3],,t]<-apply(D.beta.t,2,function(x){x*vec1})
  A.beta.curl.new[,,t]<-At[,,t]%*%solve(H.t.arr[,,t])
}





Var.curl<-function(A,t){
  G<-matrix(0,ncol(beta.m),ncol(beta.m))
  for (i in 1:401){
    A.i<-A[i,,t]
    G.i<-A.i%*%t(A.i)
    G<-G+G.i
  }
  return(G)
}




EV<-array(0,c(6,6,168))
for (t in 1:168){
  EV[,,t]<-sqrt(diag(diag(Var.curl(A.beta.curl.new,t))))
}


B.quant=matrix(0,ncol(beta.m),168)
Quant<-matrix(0,ncol(beta.m),5000)
for (j in 1:5000){
  z<-rnorm(401,0,1)
  ss<-apply(A.beta.curl.new,3,function(x){Sum<-colSums(401*z*x,2)})
  for (t in 1:168) {
    B.quant[,t]=(solve(EV[,,t]))%*%ss[,t]/sqrt(401)
  }
  Quant[,j]<-apply(B.quant,1,function(x){max(abs(x))})
}

Quant.alpha.unsm<-apply(Quant,1,function(x){quantile(x, probs=0.95)})

ucl<-matrix(0,168,6)
lcl<-matrix(0,168,6)

for (t in 1:168){
  ucl[t,]<-beta.m[t,]+EV[,,t]%*%Quant.alpha.unsm/sqrt(401)
  lcl[t,]<-beta.m[t,]-EV[,,t]%*%Quant.alpha.unsm/sqrt(401)
}


plot(time, beta.m[,4], type='l',ylim=c(-3,8),ylab="change in log-odds for svr",xlab="Days",main="Compliance for Ribavirin (Confidence Bands)")
par(new=T)
plot(time, ucl[,4], type='l',ylim=c(-3,8),col="blue",ylab="change in log-odds for svr",xlab="Days")
par(new=T)
plot(time, lcl[,4], type='l',ylim=c(-3,8),col="red",ylab="change in log-odds for svr",xlab="Days")
lines(rep(0,168))



#T2


T2<-function(C,c,W,beta){
  int<-matrix(0,length(c[[1]]),168)
  for (t in 1:168){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    W.t<-W[[t]]
    int[,t]<-(C.t%*%beta.t - c.t)*W.t
  }
  int1<- rowSums(int) - (int[,1]+int[,168])/2
  Z.sm=0
  for (t in 1:168){
    C.t<-C[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    A.t<-A.beta.curl.new[,,t]
    Int.t.sm<-W.t*A.t%*%t(C.t)
    Z.sm<-Z.sm+if(t==1 || t==168){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,1,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  t2<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
  return(pval)
}


contrast<-matrix(0,2,6)
contrast[1,2]<-1
contrast[2,3]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-contrast})
c<-vector(mode="list",168)
null.h<-matrix(0,2,1)
c<-lapply(c,function(x){x<-null.h})

W<-vector(mode="list",168)
W<-lapply(W,function(x){x<-1})



Var.Ct<-function(t,A){
  G.Ct<-matrix(0,length(C[[1]]%*%beta.m[1,]),length(C[[1]]%*%beta.m[1,]))
  
  for (i in 1:401){
    A.Ct.beta.t<-A[i,,t]
    G.i<-A.Ct.beta.t%*%t(A.Ct.beta.t)
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}



A.Ct.beta<-array(0,c(401,2,168))
for (t in 1:168){
  A.Ct.beta[,,t]<-A.beta.curl.new[,,t]%*%t(C[[t]])
}


EV.Ct<-vector("list",168)

for (t in 1:168){
  EV.Ct[[t]]<-sqrt(diag(diag(Var.Ct(t,A.Ct.beta))))
}



T2(C,c,W,beta.m)


T3<-function(C,c,beta){
  Stat<-rep(0,24)
  for (t in 1:24){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t<-(C.t%*%beta.t - c.t)		
    Cov.t<-C.t%*%Var.curl(A.beta.curl.new,t)%*%t(C.t)
    Stat[t]<-t(T.t)%*%solve(Cov.t)%*%T.t
  }
  T3<-max(Stat)
  return(T3)
}



Var.curl.z<-function(t,z,A){
  G<-matrix(0,ncol(beta.m),ncol(beta.m))
  for (i in 1:401){
    A.beta.t<-z[i]*A[i,,t]
    G.i<-A.beta.t%*%t(A.beta.t)
    G<-G+G.i
  }
  return(G)
}




B.Ct=rep(0,168)
Quant.Ct<-rep(0,5000)
Cov.t<-vector("list",168)
z<-matrix(rnorm(401*5000,0,1),5000,401)
for (j in 1:5000){
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(401*z[j,]*x,2)})
  T.Ct=ss/401
  for (t in 1:168){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.new)%*%t(C[[t]])
    B.Ct[t]=t(T.Ct[,t])%*%solve(Cov.t)%*%T.Ct[,t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3(C,c,beta.m)
p.val<-mean(Quant.Ct[1:3952]>=t.stat)






betas<-matrix(unlist(beta),168,length(beta[[1]]),byrow=T)
betas.sm<-apply(betas,2,function(x){tri_smoothed(1:168,x,3)})



beta.sm<-lapply(1:nrow(betas.sm), function(i) betas.sm[i,])








############ cumulative adherence ##############
n2.cum=1/168
n2.cum.id<-rep(0,nrow(data.to.fit))
n2.cum.id[1]<-data.to.fit$n2[1]/168
for (i in 2:nrow(data.to.fit)){
  if (data.to.fit$id[i]==data.to.fit$id[i-1]){
    n2.cum<-n2.cum+data.to.fit$n2[i]
  } else  n2.cum=data.to.fit$n2[i]
  n2.cum.id[i]<-n2.cum/168
}
n.cum=1/24
n.cum.id<-rep(0,nrow(data.to.fit2))
n.cum.id[1]<-data.to.fit2$n[1]/24
for (i in 2:nrow(data.to.fit2)){
  if (data.to.fit2$id[i]==data.to.fit2$id[i-1]){
    n.cum<-n.cum+data.to.fit2$n[i]
  } else  n.cum=data.to.fit2$n[i]
  n.cum.id[i]<-n.cum/24
}

data.to.fit$n2<-n2.cum.id
data.to.fit2$n<-n.cum.id



##########################################
## Two stage estimation ##################
##########################################

#################### Ribavirin #######################

library(rootSolve)


data1<-data.to.fit
data1<-data1[!data1$id%in%c(82,83),]
data1<-data.matrix(data1)

# data<-data[complete.cases(data),]


xmt<-demo[,c(1:7,9:16)]
beta<-vector("list",168)

for (i in 1:168){
  xvmt<-data.to.fit[data.to.fit$Day==i,]
  beta[[i]]<-glm(svr~n2+SEX+RACEW+ISHAK, data=xvmt,family = binomial(link = "logit"),na.action="na.omit")$coefficients
}



time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  y1[i]<-beta[[i]][1]
  y2[i]<-beta[[i]][2]
  y3[i]<-beta[[i]][3]
  y4[i]<-beta[[i]][4]
  y5[i]<-beta[[i]][5]
  
}


beta.init<-matrix(y1,length(y1),1)
theta.init<-c(mean(y2),mean(y3),mean(y4),mean(y5))

estimates<-est2stage(data1, beta.init, est.theta,theta.init)


###################### Peginterferon ##################################

data2<-data.to.fit2
data2<-data2[!data2$id%in%c(82,83),]
data2<-data.matrix(data2)
# data<-data[complete.cases(data),]


xmt<-demo[,c(1:7,9:16)]
beta.p<-vector("list",24)

for (i in 1:24){
  xvmt<-data.to.fit2[data.to.fit2$week==i,]
  beta.p[[i]]<-glm(svr~n+SEX+RACEW+ISHAK, data=xvmt,family = binomial(link = "logit"),na.action="na.omit")$coefficients
}



time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)
y4<-rep(0,24)
y5<-rep(0,24)


for (i in 1:24){
  y1[i]<-beta.p[[i]][1]
  y2[i]<-beta.p[[i]][2]
  y3[i]<-beta.p[[i]][3]
  y4[i]<-beta.p[[i]][4]
  y5[i]<-beta.p[[i]][5]
  
}


beta.init2<-matrix(y1,length(y1),1)
theta.init2<-c(mean(y2),mean(y3),mean(y4),mean(y5))

estimates2<-est2stage(data2, beta.init2,est.theta2 ,theta.init2)



beta_peg.sm<-tri_smoothed(1:24,estimates_peg$beta[,2],3)
beta_riba.sm<-tri_smoothed(1:168,estimates_riba$beta[,2],3)


####################################################################
### Asymptotic Distribution of Theta ###############################
####################################################################


cov.rib<-cov.theta(data,beta,theta)
cov.peg<-cov.theta(data2,beta2,theta2)

var.theta.rib<-solve(cov.rib$A)%*%cov.rib$B%*%solve(cov.rib$A)
var.theta.peg<-solve(cov.peg$A)%*%cov.peg$B%*%solve(cov.peg$A)

#F test

f1<-(395/4)*t(theta)%*%solve(var.theta.rib)%*%theta
f2<-(395/4)*t(theta2)%*%solve(var.theta.peg)%*%theta2

pf(f1,4,395,lower.tail=F)
pf(f2,4,395,lower.tail=F)

# t test

stat<-rep(0,4)
stat2<-rep(0,4)
for (i in 1:length(theta)){
  stat[i]<-sqrt(399)*theta[i]/sqrt(var.theta.rib[i,i])
  stat2[i]<-sqrt(399)*theta2[i]/sqrt(var.theta.peg[i,i])
  pt(stat,398,lower.tail=F)
  pt(stat2,398,lower.tail=F)
}

pt(stat,398,lower.tail=F)
pt(stat2,398,lower.tail=F)

pt(stat,398)
pt(stat2,398)