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

data.to.fit$n2.cum<-n2.cum.id
data.to.fit2$n.cum<-n.cum.id
data.to.fit$week<-ceiling(data.to.fit$Day/7)

dd<-merge(data.to.fit,data.to.fit2[,c(1,2,5,9)],by=c('id','week'),sort=T)

data1<-data.to.fit
data1<-data1[!data1$id%in%c(82,83),]
data1<-data.matrix(data1)
# data<-data[complete.cases(data),]


beta.n2cum.n2<-matrix(0,168,6)
beta.ncum.n<-matrix(0,24,6)
beta.n2cum<-matrix(0,168,5)
beta.cum<-matrix(0,168,6)
beta.cum.int<-matrix(0,168,7)
beta.ncum<-matrix(0,24,5)
for (i in 1:168){
    beta.n2cum.n2[i,]<-glm(svr~n2+n2.cum+SEX+RACEW+ISHAK, data=data.to.fit[data.to.fit$Day==i,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
}

for (i in 1:24){
  beta.ncum.n[i,]<-glm(svr~n+n.cum+SEX+RACEW+ISHAK, data=data.to.fit2[data.to.fit2$week==i,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
}

for (i in 1:168){
  beta.n2cum[i,]<-glm(svr~n2.cum+SEX+RACEW+ISHAK, data=data.to.fit[data.to.fit$Day==i,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
}

for (i in 1:24){
  beta.ncum[i,]<-glm(svr~n.cum+SEX+RACEW+ISHAK, data=data.to.fit2[data.to.fit2$week==i,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
}


for (i in 1:168){
  beta.cum[i,]<-glm(svr~n.cum+n2.cum+SEX+RACEW+ISHAK, data=dd[dd$Day==i,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
}


for (i in 1:168){
  beta.cum.int[i,]<-glm(svr~n.cum*n2.cum+SEX+RACEW+ISHAK, data=dd[dd$Day==i,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
}

beta.n2cum.sm<-apply(beta.n2cum,2,function(x){tri_smoothed(time,x,3)})
beta.ncum.sm<-apply(beta.ncum,2,function(x){tri_smoothed(time,x,3)})
beta.cum.sm<-apply(beta.cum,2,function(x){tri_smoothed(1:168,x,3)})
beta.cum.int.sm<-apply(beta.cum.int,2,function(x){tri_smoothed(1:168,x,3)})



time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  y1[i]<-beta.cum.int.sm[i,2]
  y2[i]<-beta.cum.int.sm[i,3]
  y3[i]<-beta.cum[i,4]
  y4[i]<-beta.cum[i,5]
  y5[i]<-beta.cum.int.sm[i,7]
}

plot(time, y1, type='l',lwd=2,ylim=c(-95,95),axes='F',ylab="change in log-odds for svr",xlab="Days",main="Cumulative adherence to Ribavirin/Peginterferon")
par(new=T)
plot(time, y2, type='l',lty=2,lwd=2,ylim=c(-95,95),axes='F',ylab="change in log-odds for svr",xlab='Days')
par(new=T)
plot(time, y5, type='l',lwd=2,lty=6,ylim=c(-95,95),axes='F',ylab="change in log-odds for svr",xlab='Days')
legend(100,85 ,lty=c(1,2,6),lwd=c(2,2,2), c('Ribavirin','Peginterferon','Interaction') )
Axis(side=2)
Axis(side=1,at=c(7,42,84,126,168))
lines(1:168,rep(0,168))



############## riba ########################

X<-subset(data.to.fit,select = c(id,svr,Day,n2,n2.cum,SEX,RACEW,ISHAK))
X<-X[!X$id%in%c(82,83),]
H.t.arr<-array(0,c(ncol(beta.n2cum.sm),ncol(beta.n2cum.sm),168))
A.beta.curl.new<-array(0,c(399,5,168))
At<-array(0,c(399,5,168))
for (t in 1:168){
  R<-as.numeric(rowSums(is.na(X[X$Day==t,]))==0)
  Z.t<-model.matrix(~n2.cum+SEX+RACEW+ISHAK ,data=X[X$Day==t,])
  D.beta.t<-apply(Z.t, 2,function(x){x*exp(Z.t%*%beta.n2cum.sm[t,])/((1+exp(Z.t%*%beta.n2cum.sm[t,]))^2)})
  Y.bar.t<- X$svr[X$Day==t]-exp(Z.t%*%beta.n2cum.sm[t,])/((1+exp(Z.t%*%beta.n2cum.sm[t,])))  
  V.beta.t<-(1+exp(Z.t%*%beta.n2cum.sm[t,]))^2/exp(Z.t%*%beta.n2cum.sm[t,])
  mat<-apply(D.beta.t,1,function(z){z%*%t(z)})
  vec<-V.beta.t*R
  Ht<-apply(mat,1,function(x){x*vec})
  H.t.arr[,,t]<-matrix(colSums(Ht),ncol(beta.n2cum.sm),ncol(beta.n2cum.sm))
  vec1<-V.beta.t*R*Y.bar.t
  At[,,t]<-apply(D.beta.t,2,function(x){x*vec1})
  A.beta.curl.new[,,t]<-At[,,t]%*%solve(H.t.arr[,,t])
}

################## Confidence Bands ######################################

Var.curl.boot<-function(t,A){return(t(A[,,t])%*%A[,,t])}
Var.boot<-function(t,A,H){return(solve(H[t,,])%*%t(A[,,t])%*%A[,,t]%*%solve(H[t,,]))}

Var.curl.sm<-function(t,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:dim(A)[1]){
    A.beta.curl.sm.t<-t(A[i,,t])
    G.i<-t(A.beta.curl.sm.t)%*%A.beta.curl.sm.t
    G<-G+G.i
  }
  return(G)
}


A.beta.curl.sm<-apply(A.beta.curl.new,c(1,2),function(x){tri_smoothed(time,x,3)})
A.beta.curl.sm1<-aperm(A.beta.curl.sm,c(2,3,1))

EV.curl.sm<-vector("list",168)
for (t in 1:168) EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm1))))

B.quant=matrix(0,ncol(beta.n2cum.sm),168)
Quant<-matrix(0,ncol(beta.n2cum.sm),5000)
for (j in 1:5000){
  z<-rnorm(399,0,1)
  ss<-apply(A.beta.curl.sm,1,function(x){Sum<-colSums(399*z*x,2)})
  for (t in 1:168){
    B.quant[,t]=(solve(EV.curl.sm[[t]]))%*%ss[,t]/sqrt(399)
  }
  Quant[,j]<-apply(B.quant,1,function(x){max(abs(x))})
}


Quant.alpha.sm<-apply(Quant,1,function(x){quantile(x, probs=0.95)})

ucl<-vector("list",168)
lcl<-vector("list",168)

for (t in 1:168){
  
  ucl[[t]]<-beta.n2cum.sm[t,]+EV[[t]]%*%Quant.alpha.sm/sqrt(401)
  lcl[[t]]<-beta.n2cum.sm[t,]-EV[[t]]%*%Quant.alpha.sm/sqrt(401)
}

time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  y1[i]<-beta.n2cum.sm[i,2]
  y2[i]<-lcl[[i]][2]
  y3[i]<-ucl[[i]][2]
  
}

plot(time, y1, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=1,xlab="Days",ylab="Cumulative Adherence",main="Cumulative adherence (Ribavirin analysis)")
par(new=T)
plot(time, y2, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=2,xlab="Days",ylab="Cumulative Adherence")
par(new=T)
plot(time, y3, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=6,xlab="Days",ylab="Cumulative Adherence")
Axis(side=2)
Axis(side=1,at=c(7,42,84,126,168))
lines(time,rep(0,length(time)))


### correction done #########

T2<-function(C,c,W,beta,A){
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
    A.t<-A[,,t]
    Int.t.sm<-W.t*A.t%*%t(C.t)
    Z.sm<-Z.sm+if(t==1 || t==168){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,1,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  t2<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
  return(pval)
}


contrast<-rep(0,5)
contrast[2]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",168)
W<-lapply(W,function(x){x<-1})




Var.Ct<-function(t,A){
  G.Ct<-matrix(0,ncol(A),ncol(A))
  
  for (i in 1:nrow(A)){
    A.Ct.beta.t<-A[i,,t]
    G.i<-A.Ct.beta.t%*%t(A.Ct.beta.t)
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}



A.Ct.beta<-array(0,c(nrow(A.beta.curl.sm1),1,168))
for (t in 1:168){
  A.Ct.beta[,,t]<-A.beta.curl.sm1[,,t]%*%t(C[[t]])
}


EV.Ct<-vector("list",168)

for (t in 1:168){
  EV.Ct[[t]]<-sqrt(diag(diag(Var.Ct(t,A.Ct.beta))))
}



T2(C,c,W,beta.n2cum.sm,A.beta.curl.sm1)


T3.t1<-function(C,c,beta,A){
  Stat<-rep(0,length(t1))
  T.t<-rep(0,length(t1))
  for (t in t1){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t[t]<-(C.t%*%beta.t - c.t)  	
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t]<-t(T.t[t])%*%solve(Cov.t)%*%T.t[t]
  }
  T3<-max(Stat)
  return(T3)
}


T3.t2<-function(C,c,beta,A){
  Stat<-rep(0,length(t2))
  T.t<-rep(0,length(t2))
  for (t in t2){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t[t-t2[1]+1]<-(C.t%*%beta.t - c.t)    
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t+1-t2[1]]<-t(T.t[t-t2[1]+1])%*%solve(Cov.t)%*%T.t[t-t2[1]+1]
  }
  T3<-max(Stat)
  return(T3)
}


Var.curl.z<-function(t,z,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:nrow(A)){
    A.beta.t<-z[i]*A[i,,t]
    G.i<-A.beta.t%*%t(A.beta.t)
    G<-G+G.i
  }
  return(G)
}

#### for t1 


B.Ct=rep(0,length(t1))
Quant.Ct<-rep(0,5000)
Cov.t<-vector("list",length(t1))
z<-matrix(rnorm(399*5000,0,1),5000,399)
for (j in 1:5000){
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(z[j,]*x,2)})
  T.Ct=ss
  for (t in t1){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t]=t(T.Ct[t])%*%solve(Cov.t)%*%T.Ct[t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3.t1(C,c,beta.n2cum.sm,A.beta.curl.sm1)
p.val<-mean(Quant.Ct>=t.stat)


### for t2

B.Ct=rep(0,length(t2))
Quant.Ct<-rep(0,5000)
Cov.t<-vector("list",length(t2))
z<-matrix(rnorm(399*5000,0,1),5000,399)
for (j in 1:5000){
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(z[j,]*x,2)})
  T.Ct=ss
  for (t in t2){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t-t2[1]+1]=t(T.Ct[t])%*%solve(Cov.t)%*%T.Ct[t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3(C,c,beta.n2cum.sm,A.beta.curl.sm1)
p.val<-mean(Quant.Ct>=t.stat)


##### significance in the first/second half of the regime #########

t1<-1:84
t2<-85:168

T2.t1<-function(C,c,W,beta,A){
  int<-matrix(0,length(c[[1]]),length(t1))
  for (t in t1){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    W.t<-W[[t]]
    int[,t]<-(C.t%*%beta.t - c.t)*W.t
  }
  int1<- rowSums(int) - (int[,t1[1]]+int[,t1[length(t1)]])/2
  Z.sm=0
  for (t in t1){
    C.t<-C[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    A.t<-A[,,t]
    Int.t.sm<-W.t*A.t%*%t(C.t)
    Z.sm<-Z.sm+if(t==t1[1] || t==t1[length(t1)]){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,1,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  tt<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(tt,nrow(as.matrix(C[[1]])))
  return(pval)
}


T2.t2<-function(C,c,W,beta,A){
  int<-matrix(0,length(c[[1]]),length(t2))
  for (t in t2){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    W.t<-W[[t]]
    int[,(t-t2[1]+1)]<-(C.t%*%beta.t - c.t)*W.t
  }
  int1<- rowSums(int) - (int[,1]+int[,length(int)])/2
  Z.sm=0
  for (t in t2){
    C.t<-C[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    A.t<-A[,,t]
    Int.t.sm<-W.t*A.t%*%t(C.t)
    Z.sm<-Z.sm+if(t==t1[1] || t==t1[length(t1)]){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,1,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  tt<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(tt,nrow(as.matrix(C[[1]])))
  return(pval)
}




contrast<-rep(0,5)
contrast[2]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",length(t1))
W<-lapply(W,function(x){x<-1})




Var.Ct<-function(t,A){
  G.Ct<-matrix(0,ncol(A),ncol(A))
  
  for (i in 1:nrow(A)){
    A.Ct.beta.t<-A[i,,t]
    G.i<-A.Ct.beta.t%*%t(A.Ct.beta.t)
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}



A.Ct.beta<-array(0,c(nrow(A.beta.curl.sm1),1,168))
for (t in 1:168){
  A.Ct.beta[,,t]<-A.beta.curl.sm1[,,t]%*%t(C[[t]])
}


EV.Ct<-vector("list",168)

for (t in 1:168){
  EV.Ct[[t]]<-sqrt(diag(diag(Var.Ct(t,A.Ct.beta))))
}



T2(C,c,W,beta.n2cum.sm,A.beta.curl.sm1)

T3.t1<-function(C,c,beta,A){
  Stat<-rep(0,length(t1))
  T.t<-rep(0,length(t1))
  for (t in t1){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t[t]<-(C.t%*%beta.t - c.t)
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t]<-t(T.t[t])%*%solve(Cov.t)%*%T.t[t]
  }
  T3<-max(Stat)
  return(T3)
}

T3.t2<-function(C,c,beta,A){
  Stat<-rep(0,length(t2))
  T.t<-rep(0,length(t2))
  for (t in t2){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t[t-t2[1]+1]<-(C.t%*%beta.t - c.t)
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t+1-t2[1]]<-t(T.t[t-t2[1]+1])%*%solve(Cov.t)%*%T.t[t-t2[1]+1]
  }
  T3<-max(Stat)
  return(T3)
}



Var.curl.z<-function(t,z,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:nrow(A)){
    A.beta.t<-z[i]*A[i,,t]
    G.i<-A.beta.t%*%t(A.beta.t)
    G<-G+G.i
  }
  return(G)
}


B.Ct=rep(0,length(t1))
Quant.Ct<-rep(0,5000)
Cov.t<-vector("list",length(t1))
z<-matrix(rnorm(399*5000,0,1),5000,399)
for (j in 1:5000){
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(z[j,]*x,2)})
  T.Ct=ss
  for (t in t1){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t]=t(T.Ct[t])%*%solve(Cov.t)%*%T.Ct[t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3.t1(C,c,beta.n2cum.sm,A.beta.curl.sm1)
p.val<-mean(Quant.Ct>=t.stat)


B.Ct=rep(0,length(t2))
Quant.Ct<-rep(0,5000)
Cov.t<-vector("list",length(t2))
z<-matrix(rnorm(399*5000,0,1),5000,399)
for (j in 1:5000){
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(z[j,]*x,2)})
  T.Ct=ss
  for (t in t2){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t-t2[1]+1]=t(T.Ct[t])%*%solve(Cov.t)%*%T.Ct[t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3.t2(C,c,beta.n2cum.sm,A.beta.curl.sm1)
p.val<-mean(Quant.Ct>=t.stat)






############### peg ###################################

data.to.fit3<-data.to.fit2[!data.to.fit2$id%in%c(82,83),]
LL<-cbind(id=unique(data.to.fit3$id),idn=1:length(unique(data.to.fit3$id)))
new.data.peg=merge(data.to.fit3,LL,by='id')


X<-subset(new.data.peg,select = c(id,idn,svr,week,n,n.cum,SEX,RACEW,ISHAK))
H.t.arr<-array(0,c(ncol(beta.ncum),ncol(beta.ncum),24))
A.beta.curl.new<-array(0,c(397,5,24))
At<-array(0,c(397,5,24))
for (t in 1:24){
  R<-as.numeric(rowSums(is.na(X[X$week==t,]))==0)
  Z.t<-model.matrix(~n.cum+SEX+RACEW+ISHAK ,data=X[X$week==t,])
  ZZ.t<-model.matrix(~id+idn+n.cum+SEX+RACEW+ISHAK ,data=X[X$week==t,])
  D.beta.t<-apply(Z.t, 2,function(x){x*exp(Z.t%*%beta.ncum[t,])/((1+exp(Z.t%*%beta.ncum[t,]))^2)})
  Y.bar.t<- X$svr[X$week==t]-exp(Z.t%*%beta.ncum[t,])/((1+exp(Z.t%*%beta.ncum[t,])))  
  V.beta.t<-(1+exp(Z.t%*%beta.ncum[t,]))^2/exp(Z.t%*%beta.ncum[t,])
  mat<-apply(D.beta.t,1,function(z){z%*%t(z)})
  vec<-V.beta.t*R
  Ht<-apply(mat,1,function(x){x*vec})
  H.t.arr[,,t]<-matrix(colSums(Ht),ncol(beta.ncum),ncol(beta.ncum))
  vec1<-V.beta.t*R*Y.bar.t
  At[ZZ.t[,3],,t]<-apply(D.beta.t,2,function(x){x*vec1})
  A.beta.curl.new[,,t]<-At[,,t]%*%solve(H.t.arr[,,t])
}


Var.curl<-function(t,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:397){
    A.beta.t<-A[i,,t]
    G.i<-A.beta.t%*%t(A.beta.t)
    G<-G+G.i
  }
  return(G)
}

A.beta.curl.sm<-apply(A.beta.curl.new,c(1,2),function(x){tri_smoothed(time,x,3)})
A.beta.curl.sm1<-aperm(A.beta.curl.sm,c(2,3,1))


Var.curl.boot<-function(t,A){return(t(A[,,t])%*%A[,,t])}
Var.boot<-function(t,A,H){return(solve(H[t,,])%*%t(A[,,t])%*%A[,,t]%*%solve(H[t,,]))}
#A.beta.curl.sm<-apply(A.beta.curl.new,c(1,2),function(x){tri_smoothed(1:24,x,3)})
#A.beta.curl.sm1<-aperm(A.beta.curl.sm,c(2,3,1))
EV.curl.sm<-vector("list",24)
for (t in 1:24) EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm))))
B.quant=matrix(0,ncol(beta.ncum),24)
Quant<-matrix(0,ncol(beta.ncum),1000)
for (j in 1:1000){
  z<-rnorm(397,0,1)
  ss<-apply(A.beta.curl.sm,1,function(x){Sum<-colSums(397*z*x,2)})
  for (t in 1:24){
    B.quant[,t]=(solve(EV.curl.sm[[t]]))%*%ss[,t]/sqrt(397)
  }
  Quant[,j]<-apply(B.quant,1,function(x){max(abs(x))})
}
Quant.alpha<-apply(Quant,1,function(x){quantile(x, probs=0.95)})

ucl.ci<-vector("list",24)
lcl.ci<-vector("list",24)
for (t in 1:24){
  lcl.ci[[t]]<-beta.ncum.sm[t,]-EV.curl.sm[[t]]%*%Quant.alpha/sqrt(397)
  ucl.ci[[t]]<-beta.ncum.sm[t,]+EV.curl.sm[[t]]%*%Quant.alpha/sqrt(397)
}


time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)
y4<-rep(0,24)
y5<-rep(0,24)



for (i in 1:24){
  y1[i]<-beta.ncum.sm[i,2]
  y2[i]<-lcl.ci[[i]][2]
  y3[i]<-ucl.ci[[i]][2]
  
}

plot(time*7, y1, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=1,xlab="Days",ylab="Cumulative Adherence",main="Cumulative adherence (Peginterferon analysis)")
par(new=T)
plot(time*7, y2, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=2,xlab="Days",ylab="Cumulative Adherence")
par(new=T)
plot(time*7, y3, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=6,xlab="Days",ylab="Cumulative Adherence")
Axis(side=2)
Axis(side=1,at=c(7,42,84,126,168))
lines(time*7,rep(0,length(time)))



T2<-function(C,c,W,beta,A){
  int<-matrix(0,length(c[[1]]),24)
  for (t in 1:24){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    W.t<-W[[t]]
    int[,t]<-(C.t%*%beta.t - c.t)*W.t
  }
  int1<- rowSums(int) - (int[,1]+int[,24])/2
  Z.sm=0
  for (t in 1:24){
    C.t<-C[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    A.t<-A[,,t]
    Int.t.sm<-W.t*A.t%*%t(C.t)
    Z.sm<-Z.sm+if(t==1 || t==24){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,1,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  t2<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
  return(pval)
}


contrast<-rep(0,5)
contrast[2]<-1
C<-vector(mode="list",24)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",24)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",24)
W<-lapply(W,function(x){x<-1})




Var.Ct<-function(t,A){
  G.Ct<-matrix(0,ncol(A),ncol(A))
  
  for (i in 1:nrow(A)){
    A.Ct.beta.t<-A[i,,t]
    G.i<-A.Ct.beta.t%*%t(A.Ct.beta.t)
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}



A.Ct.beta<-array(0,c(nrow(A.beta.curl.sm1),1,24))
for (t in 1:24){
  A.Ct.beta[,,t]<-A.beta.curl.sm1[,,t]%*%t(C[[t]])
}


EV.Ct<-vector("list",24)

for (t in 1:24){
  EV.Ct[[t]]<-sqrt(diag(diag(Var.Ct(t,A.Ct.beta))))
}



T2(C,c,W,beta.ncum.sm,A.beta.curl.sm1)


T3<-function(C,c,beta,A){
  Stat<-rep(0,24)
  for (t in 1:24){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t<-(C.t%*%beta.t - c.t)    
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t]<-t(T.t)%*%solve(Cov.t)%*%T.t
  }
  T3<-max(Stat)
  return(T3)
}



Var.curl.z<-function(t,z,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:nrow(A)){
    A.beta.t<-z[i]*A[i,,t]
    G.i<-A.beta.t%*%t(A.beta.t)
    G<-G+G.i
  }
  return(G)
}




B.Ct=rep(0,24)
Quant.Ct<-rep(0,5000)
Cov.t<-vector("list",24)
z<-matrix(rnorm(397*5000,0,1),5000,397)
for (j in 1:5000){
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(397*z[j,]*x,2)})
  T.Ct=ss/397
  for (t in 1:24){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t]=t(T.Ct[t])%*%solve(Cov.t)%*%T.Ct[t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3(C,c,beta.ncum.sm,A.beta.curl.sm)
p.val<-mean(Quant.Ct>=t.stat)



##### significance in the first/second half of the regime #########

t1<-1:12
t2<-13:24

T2.t1<-function(C,c,W,beta,A){
  int<-matrix(0,length(c[[1]]),length(t1))
  for (t in t1){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    W.t<-W[[t]]
    int[,t]<-(C.t%*%beta.t - c.t)*W.t
  }
  int1<- rowSums(int) - (int[,t1[1]]+int[,t1[length(t1)]])/2
  Z.sm=0
  for (t in t1){
    C.t<-C[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    A.t<-A[,,t]
    Int.t.sm<-W.t*A.t%*%t(C.t)
    Z.sm<-Z.sm+if(t==t1[1] || t==t1[length(t1)]){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,1,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  tt<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(tt,nrow(as.matrix(C[[1]])))
  return(pval)
}


T2.t2<-function(C,c,W,beta,A){
  int<-matrix(0,length(c[[1]]),length(t2))
  for (t in t2){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    W.t<-W[[t]]
    int[,(t-t2[1]+1)]<-(C.t%*%beta.t - c.t)*W.t
  }
  int1<- rowSums(int) - (int[,1]+int[,length(int)])/2
  Z.sm=0
  for (t in t2){
    C.t<-C[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    A.t<-A[,,t]
    Int.t.sm<-W.t*A.t%*%t(C.t)
    Z.sm<-Z.sm+if(t==t1[1] || t==t1[length(t1)]){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,1,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  tt<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(tt,nrow(as.matrix(C[[1]])))
  return(pval)
}




contrast<-rep(0,5)
contrast[2]<-1
C<-vector(mode="list",24)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",24)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",24)
W<-lapply(W,function(x){x<-1})




Var.Ct<-function(t,A){
  G.Ct<-matrix(0,ncol(A),ncol(A))
  
  for (i in 1:nrow(A)){
    A.Ct.beta.t<-A[i,,t]
    G.i<-A.Ct.beta.t%*%t(A.Ct.beta.t)
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}



A.Ct.beta<-array(0,c(nrow(A.beta.curl.sm1),1,24))
for (t in 1:24){
  A.Ct.beta[,,t]<-A.beta.curl.sm1[,,t]%*%t(C[[t]])
}


EV.Ct<-vector("list",24)

for (t in 1:24){
  EV.Ct[[t]]<-sqrt(diag(diag(Var.Ct(t,A.Ct.beta))))
}



T2(C,c,W,beta.ncum.sm,A.beta.curl.sm1)


Var.curl.z<-function(t,z,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:nrow(A)){
    A.beta.t<-z[i]*A[i,,t]
    G.i<-A.beta.t%*%t(A.beta.t)
    G<-G+G.i
  }
  return(G)
}


T3.t1<-function(C,c,beta,A){
  Stat<-rep(0,length(t1))
  T.t<-rep(0,length(t1))
  for (t in t1){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t[t]<-(C.t%*%beta.t - c.t)
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t]<-t(T.t[t])%*%solve(Cov.t)%*%T.t[t]
  }
  T3<-max(Stat)
  return(T3)
}

T3.t2<-function(C,c,beta,A){
  Stat<-rep(0,length(t2))
  T.t<-rep(0,length(t2))
  for (t in t2){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t[t-t2[1]+1]<-(C.t%*%beta.t - c.t)
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t+1-t2[1]]<-t(T.t[t-t2[1]+1])%*%solve(Cov.t)%*%T.t[t-t2[1]+1]
  }
  T3<-max(Stat)
  return(T3)
}



B.Ct=rep(0,length(t1))
Quant.Ct<-rep(0,5000)
Cov.t<-vector("list",length(t1))
z<-matrix(rnorm(397*5000,0,1),5000,397)
for (j in 1:5000){
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(z[j,]*x,2)})
  T.Ct=ss
  for (t in t1){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t]=t(T.Ct[t])%*%solve(Cov.t)%*%T.Ct[t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3.t1(C,c,beta.ncum.sm,A.beta.curl.sm)
p.val<-mean(Quant.Ct>=t.stat)


B.Ct=rep(0,length(t2))
Quant.Ct<-rep(0,5000)
Cov.t<-vector("list",length(t2))
z<-matrix(rnorm(397*5000,0,1),5000,397)
for (j in 1:5000){
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(z[j,]*x,2)})
  T.Ct=ss
  for (t in t2){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t-t2[1]+1]=t(T.Ct[t])%*%solve(Cov.t)%*%T.Ct[t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3.t2(C,c,beta.cum.sm,A.beta.curl.sm)
p.val<-mean(Quant.Ct>=t.stat)









################## peg/rib together ################################

X<-subset(dd,select = c(id,svr,Day,n,n.cum,n2,n2.cum,SEX,RACEW,ISHAK))
X<-X[!X$id%in%c(82,83),]
H.t.arr<-array(0,c(ncol(beta.cum.sm),ncol(beta.cum.sm),168))
A.beta.curl.new<-array(0,c(401,ncol(beta.cum.sm),168))
At<-array(0,c(401,ncol(beta.cum.sm),168))
for (t in 1:168){
  ZZ.t<-model.matrix(~svr+id+n.cum+n2.cum+SEX+RACEW+ISHAK ,data=X[X$Day==t,])
  R<-as.numeric(rowSums(is.na(X[X$Day==t,]))==0)
  Z.t<-model.matrix(~n.cum+n2.cum+SEX+RACEW+ISHAK ,data=X[X$Day==t,])
  D.beta.t<-apply(Z.t, 2,function(x){x*exp(Z.t%*%beta.cum[t,])/((1+exp(Z.t%*%beta.cum[t,]))^2)})
  Y.bar.t<- X$svr[X$Day==t]-exp(Z.t%*%beta.cum[t,])/((1+exp(Z.t%*%beta.cum[t,])))  
  V.beta.t<-(1+exp(Z.t%*%beta.cum[t,]))^2/exp(Z.t%*%beta.cum[t,])
  mat<-apply(D.beta.t,1,function(z){z%*%t(z)})
  vec<-V.beta.t*R
  Ht<-apply(mat,1,function(x){x*vec})
  H.t.arr[,,t]<-matrix(colSums(Ht),ncol(beta.cum.sm),ncol(beta.cum.sm))
  vec1<-V.beta.t*R*Y.bar.t
  At[ZZ.t[,3],,t]<-apply(D.beta.t,2,function(x){x*vec1})
  A.beta.curl.new[,,t]<-At[,,t]%*%solve(H.t.arr[,,t])
}

################## Confidence Bands ######################################

Var.curl.boot<-function(t,A){return(t(A[,,t])%*%A[,,t])}
Var.boot<-function(t,A,H){return(solve(H[t,,])%*%t(A[,,t])%*%A[,,t]%*%solve(H[t,,]))}

Var.curl.sm<-function(t,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:dim(A)[1]){
    A.beta.curl.sm.t<-t(A[i,,t])
    G.i<-t(A.beta.curl.sm.t)%*%A.beta.curl.sm.t
    G<-G+G.i
  }
  return(G)
}


A.beta.curl.sm<-apply(A.beta.curl.new,c(1,2),function(x){tri_smoothed(1:168,x,3)})
A.beta.curl.sm1<-aperm(A.beta.curl.sm,c(2,3,1))

EV.curl.sm<-vector("list",168)
for (t in 1:168) EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm1))))

B.quant=matrix(0,ncol(beta.cum.sm),168)
Quant<-matrix(0,ncol(beta.cum.sm),5000)
for (j in 1:5000){
  z<-rnorm(401,0,1)
  ss<-apply(A.beta.curl.sm,1,function(x){Sum<-colSums(401*z*x,2)})
  for (t in 1:168){
    B.quant[,t]=(solve(EV.curl.sm[[t]]))%*%ss[,t]/sqrt(401)
  }
  Quant[,j]<-apply(B.quant,1,function(x){max(abs(x))})
}


Quant.alpha.sm<-apply(Quant,1,function(x){quantile(x, probs=0.95)})

ucl<-vector("list",168)
lcl<-vector("list",168)

for (t in 1:168){
  
  ucl[[t]]<-beta.cum.sm[t,]+EV.curl.sm[[t]]%*%Quant.alpha.sm/sqrt(401)
  lcl[[t]]<-beta.cum.sm[t,]-EV.curl.sm[[t]]%*%Quant.alpha.sm/sqrt(401)
}

time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  y1[i]<-beta.cum.sm[i,3]
  y2[i]<-lcl[[i]][3]
  y3[i]<-ucl[[i]][3]
  
}

plot(time, y1, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=1,xlab="Days",ylab="Cumulative Adherence",main="Cumulative adherence (Ribavirin analysis)")
par(new=T)
plot(time, y2, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=2,xlab="Days",ylab="Cumulative Adherence")
par(new=T)
plot(time, y3, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=6,xlab="Days",ylab="Cumulative Adherence")
Axis(side=2)
Axis(side=1,at=c(7,42,84,126,168))
lines(time,rep(0,length(time)))


### correction done #########

T2<-function(C,c,W,beta,A){
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
    A.t<-A[,,t]
    Int.t.sm<-W.t*A.t%*%t(C.t)
    Z.sm<-Z.sm+if(t==1 || t==168){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,1,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  t2<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
  return(pval)
}


contrast<-rep(0,6)
contrast[2]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",168)
W<-lapply(W,function(x){x<-1})



Var.Ct<-function(t,A){
  G.Ct<-matrix(0,ncol(A),ncol(A))
  
  for (i in 1:nrow(A)){
    A.Ct.beta.t<-A[i,,t]
    G.i<-A.Ct.beta.t%*%t(A.Ct.beta.t)
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}



A.Ct.beta<-array(0,c(nrow(A.beta.curl.sm1),1,168))
for (t in 1:168){
  A.Ct.beta[,,t]<-A.beta.curl.sm1[,,t]%*%t(C[[t]])
}


EV.Ct<-vector("list",168)

for (t in 1:168){
  EV.Ct[[t]]<-sqrt(diag(diag(Var.Ct(t,A.Ct.beta))))
}



T2(C,c,W,beta.cum.sm,A.beta.curl.sm1)


T3<-function(C,c,beta,A){
  Stat<-rep(0,168)
  T.t<-rep(0,168)
  for (t in 1:168){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t[t]<-(C.t%*%beta.t - c.t)    
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t]<-t(T.t[t])%*%solve(Cov.t)%*%T.t[t]
  }
  T3<-max(Stat)
  return(T3)
}



Var.curl.z<-function(t,z,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:nrow(A)){
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
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(z[j,]*x,2)})
  T.Ct=ss
  for (t in 1:168){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t]=t(T.Ct[t])%*%solve(Cov.t)%*%T.Ct[t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3(C,c,beta.cum.sm,A.beta.curl.sm1)
p.val<-mean(Quant.Ct>=t.stat)


## testing for both effects


contrast<-matrix(0,6,2)
contrast[2,1]<-1
contrast[3,2]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-matrix(0,2,1)})

W<-vector(mode="list",168)
W<-lapply(W,function(x){x<-1})



Var.Ct<-function(t,A){
  G.Ct<-matrix(0,ncol(A),ncol(A))
  
  for (i in 1:nrow(A)){
    A.Ct.beta.t<-A[i,,t]
    G.i<-A.Ct.beta.t%*%t(A.Ct.beta.t)
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}



A.Ct.beta<-array(0,c(nrow(A.beta.curl.sm1),2,168))
for (t in 1:168){
  A.Ct.beta[,,t]<-A.beta.curl.sm1[,,t]%*%t(C[[t]])
}


EV.Ct<-vector("list",168)

for (t in 1:168){
  EV.Ct[[t]]<-sqrt(diag(diag(Var.Ct(t,A.Ct.beta))))
}



T2(C,c,W,beta.cum.sm,A.beta.curl.sm1)




T3<-function(C,c,beta,A){
  Stat<-rep(0,168)
  T.t<-matrix(0,2,168)
  for (t in 1:168){
    C.t<-C[[t]]
    beta.t<-t(beta[t,,drop=F])
    c.t<-c[[t]]
    T.t[,t]<-(C.t%*%beta.t - c.t)    
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t]<-t(T.t[,t])%*%solve(Cov.t)%*%T.t[,t]
  }
  T3<-max(Stat)
  return(T3)
}



Var.curl.z<-function(t,z,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:nrow(A)){
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
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(z[j,]*x,2)})
  T.Ct=ss
  for (t in 1:168){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t]=t(T.Ct[,t])%*%solve(Cov.t)%*%T.Ct[,t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3(C,c,beta.cum.sm,A.beta.curl.sm1)
p.val<-mean(Quant.Ct>=t.stat)







##################### Testing for interaction ########


X<-subset(dd,select = c(id,svr,Day,n,n.cum,n2,n2.cum,SEX,RACEW,ISHAK))
X<-X[!X$id%in%c(82,83),]
H.t.arr<-array(0,c(ncol(beta.cum.int.sm),ncol(beta.cum.int.sm),168))
A.beta.curl.new<-array(0,c(401,ncol(beta.cum.int.sm),168))
At<-array(0,c(401,ncol(beta.cum.int.sm),168))
for (t in 1:168){
  ZZ.t<-model.matrix(~svr+id+n.cum*n2.cum+SEX+RACEW+ISHAK ,data=X[X$Day==t,])
  R<-as.numeric(rowSums(is.na(X[X$Day==t,]))==0)
  Z.t<-model.matrix(~n.cum*n2.cum+SEX+RACEW+ISHAK ,data=X[X$Day==t,])
  D.beta.t<-apply(Z.t, 2,function(x){x*exp(Z.t%*%beta.cum.int[t,])/((1+exp(Z.t%*%beta.cum.int[t,]))^2)})
  Y.bar.t<- X$svr[X$Day==t]-exp(Z.t%*%beta.cum.int[t,])/((1+exp(Z.t%*%beta.cum.int[t,])))  
  V.beta.t<-(1+exp(Z.t%*%beta.cum.int[t,]))^2/exp(Z.t%*%beta.cum.int[t,])
  mat<-apply(D.beta.t,1,function(z){z%*%t(z)})
  vec<-V.beta.t*R
  Ht<-apply(mat,1,function(x){x*vec})
  H.t.arr[,,t]<-matrix(colSums(Ht),ncol(beta.cum.int.sm),ncol(beta.cum.int.sm))
  vec1<-V.beta.t*R*Y.bar.t
  At[ZZ.t[,3],,t]<-apply(D.beta.t,2,function(x){x*vec1})
  A.beta.curl.new[,,t]<-At[,,t]%*%solve(H.t.arr[,,t])
}

################## Confidence Bands ######################################

Var.curl.boot<-function(t,A){return(t(A[,,t])%*%A[,,t])}
Var.boot<-function(t,A,H){return(solve(H[t,,])%*%t(A[,,t])%*%A[,,t]%*%solve(H[t,,]))}

Var.curl.sm<-function(t,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:dim(A)[1]){
    A.beta.curl.sm.t<-t(A[i,,t])
    G.i<-t(A.beta.curl.sm.t)%*%A.beta.curl.sm.t
    G<-G+G.i
  }
  return(G)
}


A.beta.curl.sm<-apply(A.beta.curl.new,c(1,2),function(x){tri_smoothed(1:168,x,3)})
A.beta.curl.sm1<-aperm(A.beta.curl.sm,c(2,3,1))

EV.curl.sm<-vector("list",168)
for (t in 1:168) EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm1))))

B.quant=matrix(0,ncol(beta.cum.int.sm),168)
Quant<-matrix(0,ncol(beta.cum.int.sm),5000)
for (j in 1:5000){
  z<-rnorm(401,0,1)
  ss<-apply(A.beta.curl.sm,1,function(x){Sum<-colSums(401*z*x,2)})
  for (t in 1:168){
    B.quant[,t]=(solve(EV.curl.sm[[t]]))%*%ss[,t]/sqrt(401)
  }
  Quant[,j]<-apply(B.quant,1,function(x){max(abs(x))})
}


Quant.alpha.sm<-apply(Quant,1,function(x){quantile(x, probs=0.95)})

ucl<-vector("list",168)
lcl<-vector("list",168)

for (t in 1:168){
  
  ucl[[t]]<-beta.cum.int.sm[t,]+EV.curl.sm[[t]]%*%Quant.alpha.sm/sqrt(401)
  lcl[[t]]<-beta.cum.int.sm[t,]-EV.curl.sm[[t]]%*%Quant.alpha.sm/sqrt(401)
}

time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  y1[i]<-beta.cum.int.sm[i,7]
  y2[i]<-lcl[[i]][7]
  y3[i]<-ucl[[i]][7]
  
}

plot(time, y1, type='l',ylim=c(-100,300),axes=F,lwd=2,lty=1,xlab="Days",ylab="Cumulative Adherence",main="Cumulative adherence (Ribavirin analysis)")
par(new=T)
plot(time, y2, type='l',ylim=c(-100,300),axes=F,lwd=2,lty=2,xlab="Days",ylab="Cumulative Adherence")
par(new=T)
plot(time, y3, type='l',ylim=c(-100,300),axes=F,lwd=2,lty=6,xlab="Days",ylab="Cumulative Adherence")
Axis(side=2)
Axis(side=1,at=c(7,42,84,126,168))
lines(time,rep(0,length(time)))


T2<-function(C,c,W,beta,A){
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
    A.t<-A[,,t]
    Int.t.sm<-W.t*A.t%*%t(C.t)
    Z.sm<-Z.sm+if(t==1 || t==168){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,1,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  t2<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
  return(pval)
}


contrast<-rep(0,7)
contrast[2]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",168)
W<-lapply(W,function(x){x<-1})




Var.Ct<-function(t,A){
  G.Ct<-matrix(0,ncol(A),ncol(A))
  
  for (i in 1:nrow(A)){
    A.Ct.beta.t<-A[i,,t]
    G.i<-A.Ct.beta.t%*%t(A.Ct.beta.t)
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}



A.Ct.beta<-array(0,c(nrow(A.beta.curl.sm1),1,168))
for (t in 1:168){
  A.Ct.beta[,,t]<-A.beta.curl.sm1[,,t]%*%t(C[[t]])
}


EV.Ct<-vector("list",168)

for (t in 1:168){
  EV.Ct[[t]]<-sqrt(diag(diag(Var.Ct(t,A.Ct.beta))))
}



T2(C,c,W,beta.cum.sm,A.beta.curl.sm1)


T3<-function(C,c,beta,A){
  Stat<-rep(0,168)
  T.t<-rep(0,168)
  for (t in 1:168){
    C.t<-C[[t]]
    beta.t<-beta[t,]
    c.t<-c[[t]]
    T.t[t]<-(C.t%*%beta.t - c.t)    
    Cov.t<-C.t%*%Var.curl.sm(t,A)%*%t(C.t)
    Stat[t]<-t(T.t[t])%*%solve(Cov.t)%*%T.t[t]
  }
  T3<-max(Stat)
  return(T3)
}



Var.curl.z<-function(t,z,A){
  G<-matrix(0,ncol(A),ncol(A))
  for (i in 1:nrow(A)){
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
  ss<-apply(A.Ct.beta,3,function(x){Sum<-colSums(z[j,]*x,2)})
  T.Ct=ss
  for (t in 1:168){
    Cov.t<-C[[t]]%*%Var.curl.z(t,z[j,],A.beta.curl.sm1)%*%t(C[[t]])
    B.Ct[t]=t(T.Ct[t])%*%solve(Cov.t)%*%T.Ct[t]
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3(C,c,beta.cum.sm,A.beta.curl.sm1)
p.val<-mean(Quant.Ct>=t.stat)







########## cross sectional with total adherence ######### 

n2.cumsum<-rep(0,401)
n.cumsum<-rep(0,401)
for (i in 1:401){
  n2.cumsum[i]<-sum(data.to.fit$n2[data.to.fit$id==i])
  n.cumsum[i]<-sum(data.to.fit2$n[data.to.fit2$id==i])
}

demo$n2.cumsum<-n2.cumsum
demo$n.cumsum<-n.cumsum

fit.peg.riba.int<-glm(svr~n2.cumsum*n.cumsum+SEX+RACEW+ISHAK, data=demo,family = binomial(link = "logit"),na.action="na.omit")
fit.peg.riba<-glm(svr~n2.cumsum+n.cumsum+SEX+RACEW+ISHAK, data=demo,family = binomial(link = "logit"),na.action="na.omit")
fit.riba<-glm(svr~n2.cumsum+SEX+RACEW+ISHAK, data=demo,family = binomial(link = "logit"),na.action="na.omit")
fit.peg<-glm(svr~n.cumsum+SEX+RACEW+ISHAK, data=demo,family = binomial(link = "logit"),na.action="na.omit")


########## gamm() ################
library(mgcv)
library(gamm4)

gamm.fit1<-gamm4(svr~s(n2.cum)+SEX+RACEW+ISHAK,random = ~ Day|id, data=data.to.fit,family = binomial(link = "logit"),na.action="na.omit")
gamm.fit2<-gamm4(svr~n2+SEX+RACEW+ISHAK,random = ~ Day|id, data=data.to.fit,family = binomial(link = "logit"),na.action="na.omit")


gamm.fit3<-gamm4(svr~s(n.cum)+SEX+RACEW+ISHAK,random = ~ week|id, data=data.to.fit2,family = binomial(link = "logit"),na.action="na.omit")
gamm.fit4<-gamm4(svr~n+SEX+RACEW+ISHAK,random = ~ week|id, data=data.to.fit2,family = binomial(link = "logit"),na.action="na.omit")

####################################

