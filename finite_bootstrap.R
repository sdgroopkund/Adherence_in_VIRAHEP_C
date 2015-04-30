homdir<-'C:/material/coursework old/coursework/coursework/tracs/data'
setwd(homdir)

demo<-read.table('demo.txt',header=T,sep='\t')
pegifn<-read.table('pegifn.txt',header=T,sep='\t')
rib<-read.table('rib.txt',header=T,sep='\t')
vload<-read.table('vload.txt',header=T,sep='\t')

demo[demo==""]<-NA
demo$ISHAK[demo$ISHAK=="B"]<-NA
demo$ISHAK <- as.numeric(as.character(demo$ISHAK))
write.table(demo, file = "demo1.txt",col.names = TRUE,sep='\t')
demo<-read.table('demo1.txt',header=T,sep='\t')

rib<-rib[-c(26881,72914), ]

rib<-cbind(rib,n1=(rib$n==1)+(rib$n==2),n2=as.numeric(rib$n==2))


######################
## removing dropouts #
######################

rib<-rib[(rib$dropout==0 && rib$dropout2==0),]

#######################

L<-data.frame(cbind(vhcid=demo$vhcid,id=seq(1,length(demo$vhcid))))
rib_demo<-merge(rib,merge(L,demo,by='vhcid'),by='vhcid')


rib_demo1<-rib_demo[rib_demo$Day<=168,]

rib_demo21<-rib_demo[rib_demo$Day>168,]
rib_demo2<-rib_demo21[rib_demo21$response24=='Responder',]

rib_demo3<-rib_demo[rib_demo$response24=='Responder',]


################### vload baseline #######################################

vload.base<-vload[vload$time==0,]
vload.base.keep<-subset(vload.base,select=c(vhcid,vload_itt))
vload_b<-vload$vload_itt[which(vload$time==0)]

un.vb<-unique(vload.base$vhcid)
un.v<-unique(vload$vhcid)

demo<-merge(demo,vload.base.keep,by='vhcid',all.x=T)

############################################################################


peg1<-pegifn[-8449,]


######################
## removing dropouts #
######################

peg2<-peg1[is.na(peg1$do_not_include)==1,]
peg3<-peg2[is.na(peg2$do_not_include2)==1,]

#######################





L<-data.frame(cbind(vhcid=demo$vhcid,id=seq(1,length(demo$vhcid))))
peg_demo<-merge(peg3,merge(L,demo,by='vhcid'),by='vhcid')


peg_demo1<-peg_demo[peg_demo$week<=24,]

peg_demo21<-peg_demo[peg_demo$week>24,]
peg_demo2<-peg_demo21[peg_demo21$response24=='Responder',]

peg_demo3<-peg_demo[peg_demo$response24=='Responder',]






####################################
### Bootstrap in Ribavirin #########
####################################


xmt<-demo[,c(1:7,9:16)]

data.to.fit<-subset(rib_demo1,select = c(id,Day,svr,n2,SEX,RACEW,ISHAK))
data.to.fit<-cbind(data.to.fit[,1:3],offset=1,data.to.fit[,4:dim(data.to.fit)[2]])

data.to.fit2<-subset(peg_demo1,select = c(id,week,svr,n,SEX,RACEW,ISHAK))
data.to.fit2<-cbind(data.to.fit2[,1:3],offset=1,data.to.fit2[,4:dim(data.to.fit2)[2]])



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

n2.cumsum<-rep(0,401)
n.cumsum<-rep(0,401)
for (i in 1:401){
  n2.cumsum[i]<-sum(data.to.fit$n2[data.to.fit$id==i])
  n.cumsum[i]<-sum(data.to.fit2$n[data.to.fit2$id==i])
}

demo$n2.cumsum<-n2.cumsum
demo$n.cumsum<-n.cumsum


#xmt<-model.matrix(~SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=demo)


beta<-vector("list",168)
var<-vector("list",168)
beta_rob<-vector("list",168)
var_rob<-vector("list",168)

for (i in 1:168){
  beta[[i]]<-glm(svr~n2.cum+SEX+RACEW+ISHAK, data=data.to.fit[data.to.fit$Day==i,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
  #beta[[i]]<-glm(svr~n2+SEX+RACEW+ISHAK, data=data.to.fit[data.to.fit$Day==i,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
  
}


betas<-matrix(unlist(beta),168,length(beta[[1]]),byrow=T)
betas.sm<-apply(betas,2,function(x){tri_smoothed(time,x,3)})
beta.sm<-lapply(1:nrow(betas.sm), function(i) betas.sm[i,])


###############################################
### Model to test bootstrap validity ##########
###############################################


Var.curl.boot<-function(t,A){return(t(A[,,t])%*%A[,,t])}
Var.boot<-function(t,A,H){return(solve(H[t,,])%*%t(A[,,t])%*%A[,,t]%*%solve(H[t,,]))}


Var.curl.sm<-function(t,A){
  G<-matrix(0,dim(A)[3],dim(A)[3])
  for (i in 1:ncol(A)){
    A.beta.curl.sm.t<-t(A[t,i,])
    G.i<-t(A.beta.curl.sm.t)%*%A.beta.curl.sm.t
    G<-G+G.i
  }
  return(G)
}

repl<-50
expit<-function(x) 1/(1+exp(-x))
rM<-matrix(0,repl,5)
for (k in c(1:repl)){
  mu<-matrix(0,nrow=399,ncol=168)
  Y.t<-matrix(0,nrow=399,ncol=168)
  for (i in 1:168){
    xvmt<-data.to.fit[data.to.fit$Day==i,]
    #xmt1<-merge(xvmt,xmt,by="vhcid")
    X.t<-model.matrix(~n2.cum+SEX+RACEW+ISHAK ,data=xvmt)
    #X.t<-model.matrix(~n2+SEX+RACEW+ISHAK ,data=xvmt)
    Z.t<-model.matrix(~id+n2+SEX+RACEW+ISHAK ,data=xvmt)
    id.ob=Z.t[,2]
    mu[,i]=expit(X.t%*%beta.sm[[i]])
    Y.t[,i]=rbinom(399,1,mu[,i])
  }
  Y1<-data.frame(cbind(Y1=rbinom(n=nrow(mu),1,rowMeans(mu)),id=id.ob))
  Y2<-cbind(Y2=as.numeric(apply(Y.t,1,function(x){(sum(x==1)>sum(x==0))})),id=id.ob)
  new.data<-merge(data.to.fit,merge(Y1,Y2,by='id'),by='id')
  new.data$idn<-new.data$id*(new.data$id<82)+(new.data$id-2)*(new.data$id>=82)
  
  beta.y1<-matrix(0,168,5)
  beta.y2<-vector("list",168)
  for (i in 1:168){
    xvmt<-new.data[new.data$Day==i,]
    beta.y1[i,]<-glm(Y1~n2.cum+SEX+RACEW+ISHAK, data=xvmt,family = binomial(link = "logit"),na.action="na.omit")$coefficients
    #beta.y1[i,]<-glm(Y1~n2+SEX+RACEW+ISHAK, data=xvmt,family = binomial(link = "logit"),na.action="na.omit")$coefficients
    
  }
  
  
  betas.y1.sm<-apply(beta.y1,2,function(x){tri_smoothed(time,x,3)})
  beta.y1.sm<-lapply(1:nrow(betas.y1.sm), function(i) betas.y1.sm[i,])
  
#   plot(betas.sm[,2],type='l')
#   lines(1:168,betas.y1.sm[,2],type='l')
#   
  
  X<-subset(new.data,select = c(Y1,id,Day,n2,n2.cum,SEX,RACEW,ISHAK))
  H.t.arr<-array(0,c(ncol(beta.y1),ncol(beta.y1),168))
  A.beta.curl.new<-array(0,c(399,5,168))
  At<-array(0,c(399,5,168))

  for (t in 1:168){
    R<-as.numeric(rowSums(is.na(X[X$Day==t,]))==0)
    Z.t<-model.matrix(~n2.cum+SEX+RACEW+ISHAK ,data=X[X$Day==t,])
    #Z.t<-model.matrix(~n2+SEX+RACEW+ISHAK ,data=X[X$Day==t,])
    D.beta.t<-apply(Z.t, 2,function(x){x*exp(Z.t%*%beta.y1[t,])/((1+exp(Z.t%*%beta.y1[t,]))^2)})
    V.beta.t<-(1+exp(Z.t%*%beta.y1[t,]))^2/exp(Z.t%*%beta.y1[t,])
    mat<-apply(D.beta.t,1,function(z){z%*%t(z)})
    vec<-V.beta.t*R
    Ht<-apply(mat,1,function(x){x*vec})
    H.t.arr[,,t]<-matrix(colSums(Ht),ncol(beta.y1),ncol(beta.y1))
    vec1<-V.beta.t*R*Y.bar.t
    At[,,t]<-apply(D.beta.t,2,function(x){x*vec1})
    A.beta.curl.new[,,t]<-At[,,t]%*%solve(H.t.arr[,,t])
    
  }
  
  A.beta.curl.sm<-apply(A.beta.curl.new,c(1,2),function(x){tri_smoothed(time,x,3)})
  #A.beta.curl.sm1<-aperm(A.beta.curl.sm,c(2,3,1))
  
  EV.curl.sm<-vector("list",168)
  for (t in 1:168) EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm))))
  B.quant=matrix(0,ncol(beta.y1),168)
  Quant<-matrix(0,ncol(beta.y1),1000)
  for (j in 1:1000){
    z<-rnorm(399,0,1)
    ss<-apply(A.beta.curl.sm,1,function(x){Sum<-colSums(399*z*x,2)})
    for (t in 1:168){
      B.quant[,t]=(solve(EV.curl.sm[[t]]))%*%ss[,t]/sqrt(399)
    }
    Quant[,j]<-apply(B.quant,1,function(x){max(abs(x))})
  }
  Quant.alpha<-apply(Quant,1,function(x){quantile(x, probs=0.95)})
  
  ucl.ci.sm<-vector("list",168)
  lcl.ci.sm<-vector("list",168)
  for (t in 1:168){
    lcl.ci.sm[[t]]<-beta.y1.sm[[t]]-EV.curl.sm[[t]]%*%Quant.alpha/sqrt(399)
    ucl.ci.sm[[t]]<-beta.y1.sm[[t]]+EV.curl.sm[[t]]%*%Quant.alpha/sqrt(399)
  }
  rM[k,]<-rowMeans(mapply(function(x,y,z){x<=y & y<=z},lcl.ci.sm,beta.sm,ucl.ci.sm))
  
}

y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)

for (i in 1:168){
  y1[i]<-beta.sm[[i]][2]
  y2[i]<-lcl.ci.sm[[i]][2]
  y3[i]<-ucl.ci.sm[[i]][2]
  y4[i]<-beta.y1.sm[[i]][2]
}

plot(time, y1, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=1,xlab="Days",ylab="Cumulative Adherence",main="Validity of Confidence Bands (Ribavirin analysis)")
par(new=T)
plot(time, y2, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=2,xlab="Days",ylab="Cumulative Adherence")
par(new=T)
plot(time, y3, type='l',ylim=c(-20,30),axes=F,lwd=2,lty=6,xlab="Days",ylab="Cumulative Adherence")
par(new=T)
plot(time, y4, type='l',ylim=c(-20,30),axes=F,lwd=3,lty=5,xlab="Days",ylab="Cumulative Adherence")
legend(95,30,lty=c(1,6,2,5),lwd=c(2,2,2,3), c('True beta','Lower CB','Upper CB','Estimated beta'))
Axis(side=2)
Axis(side=1,at=c(7,42,84,126,168))
lines(time,rep(0,length(time)))

######################### additional formulations if needed ########
beta<-vector("list",167)
for (i in 1:167){
  beta[[i]]<-glm(svr~n2+n2.cum+SEX+RACEW+ISHAK, data=data.to.fit[data.to.fit$Day==i+1,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
}
betas<-matrix(unlist(beta),167,length(beta[[1]]),byrow=T)
betas.sm<-apply(betas,2,function(x){tri_smoothed(1:167,x,3)})
beta.sm<-lapply(1:nrow(betas.sm), function(i) betas.sm[i,])

mu<-matrix(0,nrow=399,ncol=167)
Y.t<-matrix(0,nrow=399,ncol=167)
for (i in 1:167){
  xvmt<-data.to.fit[data.to.fit$Day==i+1,]
  #xmt1<-merge(xvmt,xmt,by="vhcid")
  X.t<-model.matrix(~n2+n2.cum+SEX+RACEW+ISHAK ,data=xvmt)
  Z.t<-model.matrix(~id+n2+SEX+RACEW+ISHAK ,data=xvmt)
  id.ob=Z.t[,2]
  mu[,i]=expit(X.t%*%beta.sm[[i]])
  Y.t[,i]=rbinom(399,1,mu[,i])
}
Y1<-data.frame(cbind(Y1=rbinom(n=nrow(mu),1,rowMeans(mu)),id=id.ob))
Y2<-cbind(Y2=as.numeric(apply(Y.t,1,function(x){(sum(x==1)>sum(x==0))})),id=id.ob)
new.data<-merge(data.to.fit[data.to.fit$Day!=1,],merge(Y1,Y2,by='id'),by='id')
new.data$idn<-new.data$id*(new.data$id<82)+(new.data$id-2)*(new.data$id>=82)

beta.y1<-matrix(0,167,6)
beta.y2<-vector("list",167)
for (i in 1:167){
  xvmt<-new.data[new.data$Day==i+1,]
  beta.y1[i,]<-glm(Y1~n2+n2.cum+SEX+RACEW+ISHAK, data=xvmt,family = binomial(link = "logit"),na.action="na.omit")$coefficients
}


betas.y1.sm<-apply(beta.y1,2,function(x){tri_smoothed(1:167,x,3)})
beta.y1.sm<-lapply(1:nrow(betas.y1.sm), function(i) betas.y1.sm[i,])

plot(betas.sm[,2],type='l')
lines(1:167,betas.y1.sm[,2],type='l')

#################################################################################






####################################
### Bootstrap in Peginterferon #####
####################################

beta.peg<-vector("list",24)

for (i in 1:24){
  beta.peg[[i]]<-glm(svr~n.cum+SEX+RACEW+ISHAK, data=data.to.fit2[data.to.fit2$week==i,],family = binomial(link = "logit"),na.action="na.omit")$coefficients
}


betas.peg<-matrix(unlist(beta.peg),24,length(beta.peg[[1]]),byrow=T)
betas.peg.sm<-apply(betas.peg,2,function(x){tri_smoothed(1:24,x,3)})
beta.peg.sm<-lapply(1:nrow(betas.peg.sm), function(i) betas.peg.sm[i,])


###############################################
### Model to test bootstrap validity ##########
###############################################


Var.curl.boot<-function(t,A){return(t(A[,,t])%*%A[,,t])}
Var.boot<-function(t,A,H){return(solve(H[t,,])%*%t(A[,,t])%*%A[,,t]%*%solve(H[t,,]))}

options(warn=F)
repl<-50
expit<-function(x) 1/(1+exp(-x))
rM.peg<-matrix(0,repl,5)
for (k in c(1:repl)){
  mu.peg<-vector('list',24)
  Y.t.peg<-vector('list',24)
  for (i in 1:24){
    xvmt<-data.to.fit2[data.to.fit2$week==i,]
    #xmt1<-merge(xvmt,xmt,by="vhcid")
    X.t<-model.matrix(~n.cum+SEX+RACEW+ISHAK ,data=xvmt)
    Z.t<-model.matrix(~id+n.cum+SEX+RACEW+ISHAK ,data=xvmt)
    id.ob=Z.t[,2]
    mu.peg[[i]]=cbind(expit(X.t%*%beta.peg.sm[[i]]),id.ob)
    Y.t.peg[[i]]=cbind(rbinom(nrow(mu.peg[[i]]),1,mu.peg[[i]][1,]),id.ob)
  }
  yy<-matrix(0,401,24)
  for (i in 1:24){
    yy[mu.peg[[i]][,2],i]=mu.peg[[i]][,1]
  }
  marg.means<-rowSums(yy)/apply(yy,1,function(x){length(x[x!=0])})
  Y1.peg<-data.frame(cbind(Y1=rbinom(401,1,marg.means),id=1:401))
  #Y2.peg<-cbind(Y2=unlist(lapply(Y.t.peg,function(x){(sum(x==1)>sum(x==0))})),id=id.ob)
  new.data.peg<-merge(data.to.fit2,Y1.peg,by='id')
  new.data.peg<-new.data.peg[new.data.peg$id!=82 & new.data.peg$id!=83,]
  LL<-cbind(id=unique(new.data.peg$id),idn=1:length(unique(new.data.peg$id)))
  new.data.peg2=merge(new.data.peg,LL,by='id')
  
  beta.peg.y1<-matrix(0,24,5)
  beta.peg.y2<-vector("list",24)
  for (i in 1:24){
    xvmt<-new.data.peg2[new.data.peg2$week==i,]
    beta.peg.y1[i,]<-glm(Y1~n.cum+SEX+RACEW+ISHAK, data=xvmt,family = binomial(link = "logit"),na.action="na.omit")$coefficients
  }
  
  
  betas.peg.y1.sm<-apply(beta.peg.y1,2,function(x){tri_smoothed(1:24,x,3)})
  beta.peg.y1.sm<-lapply(1:nrow(betas.peg.y1.sm), function(i) betas.peg.y1.sm[i,])
  
  #   plot(betas.sm[,2],type='l')
  #   lines(1:168,betas.y1.sm[,2],type='l')
  #   
  
  X<-subset(new.data.peg2,select = c(Y1,id,idn,week,n,n.cum,SEX,RACEW,ISHAK))
  H.t.arr<-array(0,c(ncol(beta.peg.y1),ncol(beta.peg.y1),24))
  A.beta.curl.new<-array(0,c(397,5,24))
  At<-array(0,c(397,5,24))
  for (t in 1:24){
    R<-as.numeric(rowSums(is.na(X[X$week==t,]))==0)
    Z.t<-model.matrix(~n.cum+SEX+RACEW+ISHAK ,data=X[X$week==t,])
    ZZ.t<-model.matrix(~idn+n.cum+SEX+RACEW+ISHAK ,data=X[X$week==t,])
    D.beta.t<-apply(Z.t, 2,function(x){x*exp(Z.t%*%beta.peg.y1[t,])/((1+exp(Z.t%*%beta.peg.y1[t,]))^2)})
    Y.bar.t<- X$Y1[X$week==t]-exp(Z.t%*%beta.peg.y1[t,])/((1+exp(Z.t%*%beta.peg.y1[t,])))  
    V.beta.t<-(1+exp(Z.t%*%beta.peg.y1[t,]))^2/exp(Z.t%*%beta.peg.y1[t,])
    mat<-apply(D.beta.t,1,function(z){z%*%t(z)})
    vec<-V.beta.t*R
    Ht<-apply(mat,1,function(x){x*vec})
    H.t.arr[,,t]<-matrix(colSums(Ht),ncol(beta.peg.y1),ncol(beta.peg.y1))
    vec1<-V.beta.t*R*Y.bar.t
    At[ZZ.t[,2],,t]<-apply(D.beta.t,2,function(x){x*vec1})
    A.beta.curl.new[,,t]<-At[,,t]%*%solve(H.t.arr[,,t])
  }
  
  A.beta.curl.sm<-apply(A.beta.curl.new,c(1,2),function(x){tri_smoothed(1:24,x,3)})
  #A.beta.curl.sm1<-aperm(A.beta.curl.sm,c(2,3,1))
  EV.curl.sm<-vector("list",24)
  for (t in 1:24) EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm))))
  B.quant=matrix(0,ncol(beta.peg.y1),24)
  Quant<-matrix(0,ncol(beta.peg.y1),1000)
  for (j in 1:1000){
    z<-rnorm(397,0,1)
    ss<-apply(A.beta.curl.sm,1,function(x){Sum<-colSums(397*z*x,2)})
    for (t in 1:24){
      B.quant[,t]=(solve(EV.curl.sm[[t]]))%*%ss[,t]/sqrt(397)
    }
    Quant[,j]<-apply(B.quant,1,function(x){max(abs(x))})
  }
  Quant.alpha<-apply(Quant,1,function(x){quantile(x, probs=0.95)})
  ucl.ci.sm<-vector("list",24)
  lcl.ci.sm<-vector("list",24)
  for (t in 1:24){
    lcl.ci.sm[[t]]<-beta.peg.y1.sm[[t]]-EV.curl.sm[[t]]%*%Quant.alpha/sqrt(397)
    ucl.ci.sm[[t]]<-beta.peg.y1.sm[[t]]+EV.curl.sm[[t]]%*%Quant.alpha/sqrt(397)
  }
  rM.peg[k,]<-rowMeans(mapply(function(x,y,z){x<=y & y<=z},lcl.ci.sm.peg,beta.peg.sm,ucl.ci.sm.peg))
  
}

options(warn=T)



z1<-rep(0,24)
z2<-rep(0,24)
z3<-rep(0,24)
z4<-rep(0,24)

for (i in 1:24){
  z1[i]<-beta.peg.y1.sm[[i]][2]
  z2[i]<-lcl.ci.sm.peg[[i]][2]
  z3[i]<-ucl.ci.sm.peg[[i]][2]
  z4[i]<-beta.peg.sm[[i]][2]
}
time2<-1:24

plot(time2*7, z1, type='l',ylim=c(-5,30),axes=F,lwd=2,lty=1,xlab="Days",ylab="Cumulative Adherence",main="Validity of Confidence Bands (Peginterferon analysis)")
par(new=T)
plot(time2*7, z2, type='l',ylim=c(-5,30),axes=F,lwd=2,lty=2,xlab="Days",ylab="Cumulative Adherence")
par(new=T)
plot(time2*7, z3, type='l',ylim=c(-5,30),axes=F,lwd=2,lty=6,xlab="Days",ylab="Cumulative Adherence")
par(new=T)
plot(time2*7, z4, type='l',ylim=c(-5,30),axes=F,lwd=3,lty=5,xlab="Days",ylab="Cumulative Adherence")
legend(95,30,lty=c(1,6,2,5),lwd=c(2,2,2,3), c('True beta','Lower CB','Upper CB','Estimated beta'))
Axis(side=2)
Axis(side=1,at=c(7,42,84,126,168))
lines(time,rep(0,length(time)))










time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)

for (i in 1:168){
  y1[i]<-beta.sm[[i]][2]
  y2[i]<-lcl.ci.sm[[i]][2]
  y3[i]<-ucl.ci.sm[[i]][2]
}

plot(time, y1, type='l',ylim=c(-2,3),ylab="Race")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),col="blue",ylab="Race")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),col=2,ylab="Race")


plot(time, y1, type='l',ylim=c(-2,4),axes=F,ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin (Smoothed)")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,4),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin (Smoothed)")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,4),axes=F,col=2,ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin (Smoothed)")
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)

















time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)

for (i in 1:168){
  y1[i]<-beta.sm[[i]][2]

}

plot(time, y1, type='l',ylim=c(-2,3),ylab="change in log-odds for svr",xlab="Days",main="Compliance for Ribavirin")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),col="blue",ylab="change in log-odds for svr",xlab="Days",main="Compliance for Ribavirin")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),col=2,ylab="change in log-odds for svr",xlab="Days",main="Compliance for Ribavirin")



plot(time, y1, type='l',ylim=c(-2,3),axes=F,ylab="change in log-odds for svr",xlab="Days",main="Adherence to Ribavirin")



plot(time, y1, type='l',ylim=c(-2,3),axes=F,ylab="change in log-odds for svr",xlab="Days",main="Adherence to Ribavirin")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Days")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),axes=F,col=2,ylab="change in log-odds for svr",xlab="Days")
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)



##########################################################################

Xi<-X[X$id==i & X$Day==t,]
#Z.t<-c(1,as.numeric(Xi$n2==1),as.numeric(Xi$SEX=='Male'),as.numeric(Xi$RACEW=='Yes'),
as.numeric(Xi$MXAD=='Yes'),Xi$age,Xi$ISHAK,as.numeric(Xi$infect=='Incidental Exposure'),
as.numeric(Xi$infect=='Other'),as.numeric(Xi$infect=='Unknown'),as.numeric(Xi$education=='>= high school'),
as.numeric(Xi$insurance=='Public'),as.numeric(Xi$insurance=='Uninsured'),
as.numeric(Xi$employ=='Unemployed'),as.numeric(Xi$marital=='Married/Relationship'),
as.numeric(Xi$marital=='Never Married'),as.numeric(Xi$alcohol=='> 2 drinks/week'))


############################################################################



cov<-function(s,t){
  
  X<-subset(rib_demo1,select = c(id,Day,n2,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))
  
  time<-seq(1,168)
  V<-seq(1,168)
  G<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  H.s<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  H.t<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  A.beta.t.mean<-rep(0,length(beta[[1]]))
  A.beta.s.mean<-rep(0,length(beta[[1]]))
  kkk.t<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  
  
  
  for (i in 1:401){
    R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
    
    
    Z.t<-model.matrix(~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=X[X$id==i & X$Day==t,])
    D.beta.t<-Z.t*exp(sum(beta[[t]]*Z.t))/((1+exp(sum(beta[[t]]*Z.t)))^2)
    Y.bar.t<- rib_demo1$svr[rib_demo1$id==i & rib_demo1$Day==t]-exp(sum(beta[[t]]*Z.t))/(1+exp(sum(beta[[t]]*Z.t)))
    V.beta.t<-(1+exp(sum(beta[[t]]*Z.t)))^2/exp(sum(beta[[t]]*Z.t))
    A.beta.t<-R[t]*V.beta.t*D.beta.t*Y.bar.t
    
    
    Z.s<-model.matrix(~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=X[X$id==i & X$Day==s,])
    D.beta.s<-Z.s*exp(sum(beta[[s]]*Z.s))/(1+exp(sum(beta[[s]]*Z.s)))^2
    Y.bar.s<- rib_demo1$svr[rib_demo1$id==i & rib_demo1$Day==s]-exp(sum(beta[[s]]*Z.s))/(1+exp(sum(beta[[s]]*Z.s)))
    V.beta.s<-(1+exp(sum(beta[[s]]*Z.s)))^2/exp(sum(beta[[s]]*Z.s))
    A.beta.s<-R[s]*V.beta.s*D.beta.s*Y.bar.s
    
    
    G.i<-t(A.beta.s)%*%A.beta.t
    H.s.i<-as.matrix(R[s]*t(D.beta.s)%*%V.beta.s%*%D.beta.s)
    H.t.i<-R[t]*t(D.beta.t)%*%V.beta.t%*%D.beta.t
    
    G<-G+G.i
    H.s<-H.s+H.s.i
    H.t<-H.t+H.t.i
    
  }
  
  return(solve(H.s)%*%G%*%solve(H.t))
}

r<-cov(5,7)

#Hypothesis Testing


Cov<-array(0,c(length(beta[[1]]),length(beta[[1]]),168,168))
for (i in 1:168){
  for (j in i:168){
    Cov[,,i,j]<-cov(i,j)
  }
  
}



#T1

T1<-function(k,C,c){
  
  C.star=as.matrix(.bdiag(C[k]))
  beta.star=as.vector(unlist(beta[k]))
  c.star=unlist(c[k])
  
  #for (i in k){
  #	Sigma.j=0
  #	for (j in k){
  #		Sigma.j<-cbind(Sigma.j,(if(i<=j){Cov[,,i,j]}else{Cov[,,j,i]}))
  #	}
  #	Sigma.ij<-rbind(Sigma.ij,Sigma.j)
  #}
  
  
  Sigma.ij=0
  for (i in k){
    Sigma.j=0
    for (j in k){
      Sigma.j<-cbind(Sigma.j,cov(i,j))
    }
    Sigma.ij<-rbind(Sigma.ij,Sigma.j)
  }
  
  
  
  Sigma.star<-Sigma.ij[-1,-1]	
  
  t1<-(t(C.star%*%beta.star-c.star))%*%(solve(C.star%*%Sigma.star%*%t(C.star)))%*%(C.star%*%beta.star-c.star)
  pval=1-pchisq(t1,nrow(as.matrix(C[[1]]))*length(k))
  
  return(pval)
}

k<-seq(5,165,50)
contrast<-rep(0,17)
contrast[4]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-0})



T1(k,C,c)

############################################################################################


H<-vector("list",168)
X<-subset(rib_demo1,select = c(id,Day,n2,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))

for (t in 1:168){
  H.t<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  for (i in 1:401){
    if (nrow(X[X$id==i & X$Day==t,])!=0){
      R<-as.numeric(sum(is.na(X[X$id==i & X$Day==t,]))==0)
      #R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
      Z.t<-model.matrix(~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=X[X$id==i & X$Day==t,])
      D.beta.t<-Z.t*exp(sum(beta[[t]]*Z.t))/((1+exp(sum(beta[[t]]*Z.t)))^2)
      Y.bar.t<- rib_demo1$svr[rib_demo1$id==i & rib_demo1$Day==t]-exp(sum(beta[[t]]*Z.t))/(1+exp(sum(beta[[t]]*Z.t)))
      V.beta.t<-(1+exp(sum(beta[[t]]*Z.t)))^2/exp(sum(beta[[t]]*Z.t))
      H.t.i<-R*t(D.beta.t)%*%V.beta.t%*%D.beta.t
      H.t<-H.t+H.t.i
    }
  }		
  
  H[[t]]<-H.t
}

############### SANDWICH ESTIMATOR ###############################3

H_rob<-vector("list",168)
X<-subset(rib_demo1,select = c(id,Day,n2,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))

for (t in 1:168){
  H_rob.t<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  for (i in 1:401){
    if (nrow(X[X$id==i & X$Day==t,])!=0){
      R<-as.numeric(sum(is.na(X[X$id==i & X$Day==t,]))==0)
      #R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
      Z.t<-model.matrix(~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=X[X$id==i & X$Day==t,])
      D.beta.t<-Z.t*exp(sum(beta[[t]]*Z.t))/((1+exp(sum(beta[[t]]*Z.t)))^2)
      Y.bar.t<- rib_demo1$svr[rib_demo1$id==i & rib_demo1$Day==t]-exp(sum(beta[[t]]*Z.t))/(1+exp(sum(beta[[t]]*Z.t)))
      V.beta.t<-(1+exp(sum(beta[[t]]*Z.t)))^2/exp(sum(beta[[t]]*Z.t))
      zz<-R*t(D.beta.t)%*%V.beta.t_rob%*%D.beta.t
      zzz<-solve(zz)%*%
        H.t.i<-zz
      H_rob.t<-H_rob.t+H.t.i
    }
  }		
  
  H_rob[[t]]<-H_rob.t
}


##############################################################################



A.beta<-vector("list",168)
A.beta.curl<-vector("list",168)
X<-subset(rib_demo1,select = c(id,Day,n2,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))

for (i in 1:401)
{
  A.beta.i<-matrix(0,168,17)
  A.beta.i.curl<-matrix(0,168,17)
  for (t in 1:168)
  {
    if (nrow(X[X$id==i & X$Day==t,])!=0){
      R<-as.numeric(sum(is.na(X[X$id==i & X$week==t,]))==0)
      #R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
      Z.t<-model.matrix(~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=X[X$id==i & X$Day==t,])
      D.beta.t<-Z.t*exp(sum(beta[[t]]*Z.t))/((1+exp(sum(beta[[t]]*Z.t)))^2)
      Y.bar.t<- rib_demo1$svr[rib_demo1$id==i & rib_demo1$Day==t]-exp(sum(beta[[t]]*Z.t))/(1+exp(sum(beta[[t]]*Z.t)))
      V.beta.t<-(1+exp(sum(beta[[t]]*Z.t)))^2/exp(sum(beta[[t]]*Z.t))
      A<-R*V.beta.t*D.beta.t*Y.bar.t
      A.beta.i[t,]<-if (length(A)==0) {rep(0,17)} else {A}
      A.beta.i.curl[t,]<-solve(H[[t]])%*%A.beta.i[t,]
    }
  }
  A.beta[[i]]<-A.beta.i
  A.beta.curl[[i]]<-A.beta.i.curl
}







Var<-function(t){
  G<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  for (i in 1:401){
    A.beta.t<-t(A.beta[[i]][t,])
    G.i<-t(A.beta.t)%*%A.beta.t
    G<-G+G.i
  }
  return(solve(H[[t]])%*%G%*%solve(H[[t]]))
}

ucl.ci<-vector("list",168)
lcl.ci<-vector("list",168)
for (t in 1:168){
  lcl.ci[[t]]<-beta[[t]]-sqrt(diag(Var(t)))*1.96
  ucl.ci[[t]]<-beta[[t]]+sqrt(diag(Var(t)))*1.96
}




time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)

for (i in 1:168){
  y1[i]<-beta[[i]][2]
  y2[i]<-lcl.ci[[i]][2]
  y3[i]<-ucl.ci[[i]][2]
}

plot(time, y1, type='l',ylim=c(-2,3),ylab="Race")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),col="blue",ylab="Race")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),col=2,ylab="Race")


plot(time, y1, type='l',ylim=c(-2,3),axes=F,ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),axes=F,col=2,ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin")
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)




zz<-ksmooth(time, y1, kernel = c("box", "normal"), bandwidth = 0.5,
            range.x = range(time),
            n.points = max(100, length(time)), time)

lines(ksmooth(time, y1, kernel = "normal", bandwidth = 4),col=2)
pp<-ksmooth(time, y1, kernel = "normal", bandwidth = 4)




################ Function to induce Triangular Kernel ###########################

tri_smoothed<-function(time,y,bandwidth){
  f0<-rep(0,length(y))
  D<-function(t){
    return(ifelse(abs(t)>1,0,1-abs(t)))
  }
  K_b<-function(x,y){
    return(D((x-y)/bandwidth))
  }
  for (t in time){
    f0[t]<-sum(K_b(t,time)*y)/sum(K_b(t,time))
  }
  return(f0)
}



betas<-matrix(unlist(beta),168,17,byrow=T)
betas.sm<-apply(betas,2,function(x){tri_smoothed(time,x,3)})
beta.sm<-lapply(1:nrow(betas.sm), function(i) betas.sm[i,])


plot(time, y1, type='l',ylim=c(-2,4),ylab="Race")
lines(tri_smoothed(time,y1,3),col=2)
par(new=T)
plot(time, y2, type='l',ylim=c(-2,4),col="blue",ylab="Race")
lines(tri_smoothed(time,y2,3),col=3)
par(new=T)
plot(time, y3, type='l',ylim=c(-2,4),col=2,ylab="Race")
lines(tri_smoothed(time,y3,3),col=8)

time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  y1[i]<-beta[[i]][2]
  y2[i]<-lcl[[i]][2]
  y3[i]<-ucl[[i]][2]
  y4[i]<-cl_l[[i]][2]
  y5[i]<-cl_u[[i]][2] 
}


plot(time, y1, type='l',ylim=c(-2,4),axes=F,ylab="change in log-odds for svr",xlab="Days",main="Adherence to Ribavirin")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,4),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Days")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,4),axes=F,col=2,ylab="change in log-odds for svr",xlab="Days")
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)


plot(time, y1, type='l',ylim=c(-6,1),ylab="Intercept")
par(new=T)
plot(time, y2, type='l',ylim=c(-6,1),col="red",ylab="Intercept")
par(new=T)
plot(time, y3, type='l',ylim=c(-6,1),col="red",ylab="Intercept")
par(new=T)
plot(time, y4, type='l',ylim=c(-6,1),col="blue",ylab="Intercept")
par(new=T)
plot(time, y5, type='l',ylim=c(-6,1),col="blue",ylab="Intercept")





plot(time, y1, type='l',ylim=c(-6,1),ylab="Intercept")
lines(tri_smoothed(time,y1,3),col=2)
par(new=T)
plot(time, y4, type='l',ylim=c(-6,1),col="blue",ylab="Intercept")
lines(tri_smoothed(time,y4,3),col=2)
par(new=T)
plot(time, y5, type='l',ylim=c(-6,1),col="red",ylab="Intercept")
lines(tri_smoothed(time,y5,3),col=1)



plot(time, tri_smoothed(time,y1,3), type='l',ylim=c(-6,1),ylab="Intercept")
lines(tri_smoothed(time,y4,3),col=2)
lines(tri_smoothed(time,y5,3),col=3)

plot(time, y1, type='l',ylim=c(-2,3),ylab="Compliance")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),col="red",ylab="Compliance")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),col="red",ylab="Compliance")
par(new=T)
plot(time, y4, type='l',ylim=c(-2,3),col="blue",ylab="Compliance")
par(new=T)
plot(time, y5, type='l',ylim=c(-2,3),col="blue",ylab="Compliance")


plot(time, y1, type='l',ylim=c(-2,3),ylab="Compliance")
lines(tri_smoothed(time,y1,3),col=2)
par(new=T)
plot(time, y4, type='l',ylim=c(-2,3),col="blue",ylab="Compliance")
lines(tri_smoothed(time,y4,3),col=2)
par(new=T)
plot(time, y5, type='l',ylim=c(-2,3),col="blue",ylab="Compliance")
lines(tri_smoothed(time,y5,3),col=1)





plot(time, y1, type='l',ylim=c(-2,1),ylab="Sex")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,1),col="red",ylab="Sex")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,1),col="red",ylab="Sex")
par(new=T)
plot(time, y4, type='l',ylim=c(-2,1),col="blue",ylab="Sex")
par(new=T)
plot(time, y5, type='l',ylim=c(-2,1),col="blue",ylab="Sex")



plot(time, y1, type='l',ylim=c(-2,1),ylab="Sex")
lines(tri_smoothed(time,y1,3),col=2)
par(new=T)
plot(time, y4, type='l',ylim=c(-2,1),col="blue",ylab="Sex")
lines(tri_smoothed(time,y4,3),col=2)
par(new=T)
plot(time, y5, type='l',ylim=c(-2,1),col="blue",ylab="Sex")
lines(tri_smoothed(time,y5,3),col=1)




plot(time, y1, type='l',ylim=c(0,2),ylab="Race")
par(new=T)
plot(time, y2, type='l',ylim=c(0,2),col="red",ylab="Race")
par(new=T)
plot(time, y3, type='l',ylim=c(0,2),col="red",ylab="Race")
par(new=T)
plot(time, y4, type='l',ylim=c(0,2),col="blue",ylab="Race")
par(new=T)
plot(time, y5, type='l',ylim=c(0,2),col="blue",ylab="Race")




plot(time, y1, type='l',ylim=c(-2,4),ylab="Race")
lines(tri_smoothed(time,y1,3),col=2)
par(new=T)
plot(time, y2, type='l',ylim=c(-2,4),col="blue",ylab="Race")
lines(tri_smoothed(time,y2,3),col=2)
par(new=T)
plot(time, y3, type='l',ylim=c(-2,4),col="blue",ylab="Race")
lines(tri_smoothed(time,y3,3),col=2)



plot(time, y1, type='l',ylim=c(-0.5,1.5),ylab="Education")
par(new=T)
plot(time, y2, type='l',ylim=c(-0.5,1.5),col="red",ylab="Education")
par(new=T)
plot(time, y3, type='l',ylim=c(-0.5,1.5),col="red",ylab="Education")
par(new=T)
plot(time, y4, type='l',ylim=c(-0.5,1.5),col="blue",ylab="Education")
par(new=T)
plot(time, y5, type='l',ylim=c(-0.5,1.5),col="blue",ylab="Education")



plot(time, y1, type='l',ylim=c(-0.5,0),ylab="ISHAAK")
par(new=T)
plot(time, y2, type='l',ylim=c(-0.5,0),col="red",ylab="ISHAAK")
par(new=T)
plot(time, y3, type='l',ylim=c(-0.5,0),col="red",ylab="ISHAAK")
par(new=T)
plot(time, y4, type='l',ylim=c(-0.5,0),col="blue",ylab="ISHAAK")
par(new=T)
plot(time, y5, type='l',ylim=c(-0.5,0),col="blue",ylab="ISHAAK")


A.beta2<-A.beta
A.beta.curl2<-A.beta.curl

A.beta.sm<-lapply(A.beta2,function(x){apply(x,2,function(x){tri_smoothed(time,x,3)})})
A.beta.curl.sm<-lapply(A.beta.curl2,function(x){apply(x,2,function(x){tri_smoothed(time,x,3)})})

#################################################################################################################

Var.sm<-function(t,A){
  
  X<-subset(rib_demo1,select = c(id,Day,n2,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))
  
  time<-seq(1,168)
  G<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  H.t<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  
  
  for (i in 1:401){
    R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
    Z.t<-model.matrix(~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=X[X$id==i & X$Day==t,])
    D.beta.t<-Z.t*exp(sum(beta[[t]]*Z.t))/((1+exp(sum(beta[[t]]*Z.t)))^2)
    Y.bar.t<- rib_demo1$svr[rib_demo1$id==i & rib_demo1$Day==t]-exp(sum(beta[[t]]*Z.t))/(1+exp(sum(beta[[t]]*Z.t)))
    V.beta.t<-(1+exp(sum(beta[[t]]*Z.t)))^2/exp(sum(beta[[t]]*Z.t))
    
    
    A.beta.sm.t<-t(A[[i]][t,])
    
    G.i<-t(A.beta.sm.t)%*%A.beta.sm.t
    H.t.i<-R[t]*t(D.beta.t)%*%V.beta.t%*%D.beta.t
    G<-G+G.i
    H.t<-H.t+H.t.i
    
  }
  
  return(solve(H.t)%*%G%*%solve(H.t))
}


###############################################################################################################



Var.curl.sm<-function(t,A){
  X<-subset(rib_demo1,select = c(id,Day,n2,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))
  time<-seq(1,168)
  G<-matrix(0,length(beta[[1]]),length(beta[[1]]))
  
  for (i in 1:401){
    A.beta.curl.sm.t<-t(A[[i]][t,])
    G.i<-t(A.beta.curl.sm.t)%*%A.beta.curl.sm.t
    G<-G+G.i
  }
  return(G)
}



############################### smoothed ci ################################


ucl.ci.sm<-vector("list",168)
lcl.ci.sm<-vector("list",168)
for (t in 1:168){
  lcl.ci.sm[[t]]<-beta.sm[[t]]-sqrt(diag(Var.curl.sm(t,A.beta.curl.sm)))*1.96
  ucl.ci.sm[[t]]<-beta.sm[[t]]+sqrt(diag(Var.curl.sm(t,A.beta.curl.sm)))*1.96
}




time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)

for (i in 1:168){
  y1[i]<-beta.sm[[i]][2]
  y2[i]<-lcl.ci.sm[[i]][2]
  y3[i]<-ucl.ci.sm[[i]][2]
}

plot(time, y1, type='l',ylim=c(-2,3),ylab="Race")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),col="blue",ylab="Race")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),col=2,ylab="Race")


plot(time, y1, type='l',ylim=c(-2,4),axes=F,ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin (Smoothed)")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,4),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin (Smoothed)")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,4),axes=F,col=2,ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin (Smoothed)")
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)


################################################################################





EV.curl.sm<-vector("list",168)
for (t in 1:168){
  EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm))))
  
}

B.curl.sm1=matrix(0,length(beta[[1]]),168)
Quant<-matrix(0,17,5000)
for (j in 1:5000){
  z<-rnorm(401,0,1)
  for (t in 1:168)
  {
    Sum.curl.sm=rep(0,length(beta[[1]]))
    for (i in 1:401)
    {
      A.beta.curl.sm.t<-t(A.beta.curl.sm[[i]][t,])
      Sum.curl.sm=Sum.curl.sm+ if (length(A.beta.curl.sm.t)==0) {0} else {401*z[i]*t(A.beta.curl.sm.t)}
    }
    #EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm))))
    B.curl.sm1[,t]=(solve(EV.curl.sm[[t]]))%*%Sum.curl.sm/sqrt(401)
  }
  Quant[,j]<-apply(B.curl.sm1,1,function(x){max(abs(x))})
}


Quant.alpha<-apply(Quant,1,function(x){quantile(x, probs=0.95)})


ucl.curl.sm<-vector("list",168)
lcl.curl.sm<-vector("list",168)


for (t in 1:168){
  
  ucl.curl.sm[[t]]<-beta.sm[[t]]+EV.curl.sm[[t]]%*%Quant.alpha/sqrt(401)
  lcl.curl.sm[[t]]<-beta.sm[[t]]-EV.curl.sm[[t]]%*%Quant.alpha/sqrt(401)
}

time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  y1[i]<-beta.sm[[i]][2]
  y2[i]<-lcl.curl.sm[[i]][2]
  y3[i]<-ucl.curl.sm[[i]][2]
  
}






plot(time, y1, type='l',ylim=c(-2,3),ylab="change in log-odds for svr",xlab="Days",main="Compliance for Ribavirin (Hypothesis Tests)")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),col="blue",ylab="change in log-odds for svr",xlab="Days")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),col="red",ylab="change in log-odds for svr",xlab="Days")
lines(rep(0,168),lwd=2)



########################################################################################################



EV<-vector("list",168)
for (t in 1:168){
  EV[[t]]<-sqrt(diag(diag(Var(t))))
}





B.quant=matrix(0,length(beta[[1]]),168)
Quant<-matrix(0,17,5000)
for (j in 1:5000){
  z<-rnorm(401,0,1)
  for (t in 1:168)
  {
    Sum=rep(0,length(beta[[1]]))
    for (i in 1:401)
    {
      A.beta.t<-t(A.beta[[i]][t,])
      Sum=Sum+ if (length(A.beta.t)==0) {0} else {z[i]*(401*var[[t]])%*%t(A.beta.t)}
    }
    #EV[[t]]<-sqrt(diag(diag(Var(t))))
    B.quant[,t]=(solve(EV[[t]]))%*%Sum/sqrt(401)
  }
  Quant[,j]<-apply(B.quant,1,function(x){max(abs(x))})
}


Quant.alpha.unsm<-apply(Quant,1,function(x){quantile(x, probs=0.95)})

ucl<-vector("list",168)
lcl<-vector("list",168)

for (t in 1:168){
  
  ucl[[t]]<-beta[[t]]+EV[[t]]%*%Quant.alpha.unsm/sqrt(401)
  lcl[[t]]<-beta[[t]]-EV[[t]]%*%Quant.alpha.unsm/sqrt(401)
}






time<-seq(1,168)
y<-rep(0,168)
y1<-rep(0,168)
y12<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  
  y1[i]<-beta[[i]][3]
  y4[i]<-lcl[[i]][3]
  y5[i]<-ucl[[i]][3] 
}




plot(time, y1, type='l',ylim=c(-2,4),ylab="change in log-odds for svr",xlab="Days",main="Compliance for Ribavirin (Confidence Bands)")
par(new=T)
plot(time, y4, type='l',ylim=c(-2,4),col="blue",ylab="change in log-odds for svr",xlab="Days")
par(new=T)
plot(time, y5, type='l',ylim=c(-2,4),col="red",ylab="change in log-odds for svr",xlab="Days")
lines(rep(0,168))



for (i in 1:168){
  y1[i]<-beta.sm[[i]][4]
  y2[i]<-lcl.curl.sm[[i]][4]
  y3[i]<-ucl.curl.sm[[i]][4]
}

plot(time, y1, type='l',ylim=c(-2,3),axes=F,ylab="change in log-odds for svr",xlab="Days",main="Adherence to Ribavirin (Hypothesis Test)")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Days")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),axes=F,col=2,ylab="change in log-odds for svr",xlab="Days")
lines(rep(0,168))
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)


plot(time, y1, type='l',ylim=c(-1,1),axes=F,ylab="change in log-odds for svr",xlab="Days",main="Effect of ISHAK on SVR in Ribavirin Analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-1,1),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Days")
par(new=T)
plot(time, y3, type='l',ylim=c(-1,1),axes=F,col=2,ylab="change in log-odds for svr",xlab="Days")
lines(rep(0,168))
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)

plot(time, y1, type='l',ylim=c(-2,2),axes=F,ylab="change in log-odds for males v females on svr",xlab="Days",main="Effect of Gender on SVR in Ribavirin Analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,2),axes=F,col="blue",ylab="change in log-odds for males v females on svr",xlab="Days")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,2),axes=F,col=2,ylab="change in log-odds for males v females on svr",xlab="Days")
lines(rep(0,168))
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)

plot(time, y1, type='l',ylim=c(-2,3),axes=F,ylab="change in log-odds for caucasians v non-caucasians on svr",xlab="Days",main="Effect of Race on SVR in Ribavirin Analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),axes=F,col="blue",ylab="change in log-odds for caucasians v non-caucasians on svr",xlab="Days")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),axes=F,col=2,ylab="change in log-odds for caucasians v non-caucasians on svr",xlab="Days")
lines(rep(0,168))
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)


save.image("G:\\inside\\coursework\\tracs\\data\\WS_rib_n2_f24")














#T2

T2<-function(C,c,W){
  
  int<-matrix(0,length(c[[1]]),168)
  
  for (t in 1:168){
    
    C.t<-C[[t]]
    beta.t<-beta[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    
    int[,t]<-(C.t%*%beta.t - c.t)*W.t
  }
  
  int1<- rowSums(int) - (int[,1]+int[,168])/2
  Z=0
  for (t in 1:168){
    H.t<-matrix(0,length(beta[[1]]),length(beta[[1]]))
    C.t<-C[[t]]
    beta.t<-beta[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    A.t<-matrix(0,length(beta[[1]]),401)
    for (i in 1:401){
      A.t[,i]<-A.beta[[i]][t,]
    }
    Int.t<-C.t%*%solve(H[[t]]/401)%*%A.t*W.t
    Z<-Z+if(t==1 || t==168){Int.t/2}else{Int.t}
  }
  sigma<-(1/(401^2))*matrix(rowSums(rbind(apply(Z,2,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  #Rsum<-rowSums(as.matrix(Z))
  #sigma<-(1/(401^2))*t(Rsum)%*%Rsum
  
  t2<-t(int1)%*%solve(sigma)%*%int1
  pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
  return(pval)
}

contrast<-rep(0,17)
contrast[6]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",168)
W<-lapply(W,function(x){x<-1})
T2(C,c,W)


contrast<-matrix(0,2,17)
contrast[1,15]<-1
contrast[2,16]<-1
contrast[3,10]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-contrast})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-rep(0,2)})
T2(C,c,W)








T2.sm<-function(C,c,W){
  
  int<-matrix(0,length(c[[1]]),168)
  for (t in 1:168){
    C.t<-C[[t]]
    beta.sm.t<-beta.sm[[t]]
    beta.t<-beta[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    int[,t]<-(C.t%*%beta.sm.t - c.t)*W.t
  }
  int1<- rowSums(int) - (int[,1]+int[,168])/2
  Z.sm=0
  for (t in 1:168){
    C.t<-C[[t]]
    c.t<-c[[t]]
    W.t<-W[[t]]
    A.t<-matrix(0,length(beta[[1]]),401)
    for (i in 1:401){
      A.t[,i]<-A.beta.curl[[i]][t,]
    }
    Int.t.sm<-C.t%*%A.t*W.t
    Z.sm<-Z.sm+if(t==1 || t==168){Int.t.sm/2}else{Int.t.sm}
  }
  sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,2,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
  t2<-t(int1)%*%solve(sigma.sm)%*%int1
  pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
  return(pval)
}

k<-seq(5,165,50)


contrast<-rep(0,17)
contrast[8]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",168)
W<-lapply(W,function(x){x<-1})

T2.sm(C,c,W)












T3<-function(C,c){
  Stat<-rep(0,168)
  for (t in 1:168){
    C.t<-C[[t]]
    beta.t<-beta[[t]]
    c.t<-c[[t]]
    T.t<-(C.t%*%beta.t - c.t)		
    Cov.t<-C.t%*%Cov[,,t,t]%*%t(C.t)
    Stat[t]<-t(T.t)%*%Cov.t%*%T.t
  }
  T3<-max(abs(Stat))
  return(T3)
}














k<-seq(5,165,50)
contrast<-rep(0,17)
contrast[2]<-1
C<-vector(mode="list",168)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",168)
c<-lapply(c,function(x){x<-0})



A.Ct.beta<-lapply(A.beta.curl2,function(x){apply(x,1,function(z){C[[t]]%*%z})})
A.Ct.beta.sm<-lapply(A.Ct.beta,function(x){apply(cbind(x),2,function(z){tri_smoothed(time,z,3)})})



Var.Ct.sm<-function(t,A){
  X<-subset(rib_demo1,select = c(id,Day,n2,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))
  time<-seq(1,168)
  G.Ct<-matrix(0,length(C[[1]]%*%beta[[1]]),length(C[[1]]%*%beta[[1]]))
  
  for (i in 1:401){
    A.Ct.beta.sm.t<-t(A[[i]][t,])
    G.i<-t(A.Ct.beta.sm.t)%*%A.Ct.beta.sm.t
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}

EV.Ct<-vector("list",168)

for (t in 1:168){
  EV.Ct[[t]]<-sqrt(diag(Var.Ct.sm(t,A.Ct.beta.sm)))
}


B.Ct=matrix(0,length(C[[1]]%*%beta[[1]]),168)

Quant.Ct<-matrix(0,length(C[[1]]%*%beta[[1]]),5000)
for (j in 1:5000){
  z<-rnorm(401,0,1)
  for (t in 1:168)
  {
    Sum.Ct=rep(0,length(C[[1]]%*%beta[[1]]))
    for (i in 1:401)
    {
      A.Ct.beta.t<-t(A.Ct.beta.sm[[i]][t,])
      Sum.Ct=Sum.Ct+ if (length(A.Ct.beta.t)==0) {0} else {401*z[i]*t(A.Ct.beta.t)}
    }
    #EV.Ct[[t]]<-sqrt(diag(diag(Var.Ct.sm(t,A.Ct.beta.sm))))
    B.Ct[,t]=(solve(EV.Ct[[t]]))%*%Sum.Ct/sqrt(401)
  }
  Quant.Ct[,j]<-apply(B.Ct,1,function(x){max(abs(x))})
}


Quant.Ct.alpha<-apply(Quant.Ct,1,function(x){quantile(x, probs=0.95)})


ucl.Ct.sm<-vector("list",168)
lcl.Ct.sm<-vector("list",168)


for (t in 1:168){
  
  ucl.Ct.sm[[t]]<-(C[[t]]%*%beta.sm[[t]]-c[[t]])+EV.Ct[[t]]%*%Quant.Ct.alpha[1]/sqrt(401)
  lcl.Ct.sm[[t]]<-(C[[t]]%*%beta.sm[[t]]-c[[t]])-EV.Ct[[t]]%*%Quant.Ct.alpha[1]/sqrt(401)
}




time<-seq(1,168)
y1<-rep(0,168)
y2<-rep(0,168)
y3<-rep(0,168)
y4<-rep(0,168)
y5<-rep(0,168)


for (i in 1:168){
  y1[i]<-beta.sm[[i]][2]
  y2[i]<-lcl.curl.sm[[i]][2]
  y3[i]<-ucl.curl.sm[[i]][2]
  y4[i]<-lcl.Ct.sm[[i]]
  y5[i]<-ucl.Ct.sm[[i]] 
}




plot(time, y1, type='l',ylim=c(-2,3),ylab="Compliance")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),col="red",ylab="Compliance")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),col="red",ylab="Compliance")
par(new=T)
plot(time, y4, type='l',ylim=c(-2,3),col="blue",ylab="Compliance")
par(new=T)
plot(time, y5, type='l',ylim=c(-2,3),col="blue",ylab="Compliance")





save.image("G:\\inside\\coursework\\tracs\\data\\anls6")



# Cross Validation



betas<-matrix(unlist(beta),168,17,byrow=T)


Mse<-matrix(0,100,168)
V<-matrix(0,100,168)
B<-matrix(0,100,168)
for (h in 1:100){
  betas.sm.h<-apply(betas,2,function(x){tri_smoothed(time,x,h)})
  beta.h<-lapply(1:nrow(betas.sm.h), function(i) betas.sm.h[i,])
  A.h<-lapply(A.beta.curl2,function(x){apply(x,2,function(x){tri_smoothed(time,x,h)})})
  for (t in 1:168){
    var.beta.h.t<-sum(diag(Var.curl.sm(t,A.h)))
    V[h,t]<-var.beta.h.t
    sq_bias.h.t<-sum((beta[[t]]-beta.h[[t]])^2)
    B[h,t]<-sq_bias.h.t
    Mse[h,t]<-var.beta.h.t+sq_bias.h.t
  }
}



betas<-matrix(unlist(beta),168,17,byrow=T)
betas.sm2<-apply(betas,2,function(x){tri_smoothed(time,x,2)})
beta.sm2<-lapply(1:nrow(betas.sm2), function(i) betas.sm2[i,])














A.Ct.beta<-lapply(A.beta.curl2,function(x){apply(x,1,function(z){C[[t]]%*%z})})


Var.Ct<-function(t,A){
  time<-seq(1,24)
  G.Ct<-matrix(0,length(C[[1]]%*%beta[[1]]),length(C[[1]]%*%beta[[1]]))
  
  for (i in 1:401){
    A.Ct.beta.t<-t(t(A[[i]])[t,])
    G.i<-t(A.Ct.beta.t)%*%A.Ct.beta.t
    G.Ct<-G.Ct+G.i
  }
  return(G.Ct)
}

EV.Ct<-vector("list",24)

for (t in 1:24){
  EV.Ct[[t]]<-Var.Ct(t,A.Ct.beta)
}




B.Ct=rep(0,24)
Quant.Ct<-rep(0,5000)
for (j in 1:5000){
  z<-rnorm(401,0,1)
  for (t in 1:168)
  {
    Sum.Ct=rep(0,length(C[[1]]%*%beta[[1]]))
    for (i in 1:401)
    {
      A.Ct.beta.t<-t(t(A.Ct.beta[[i]])[t,])
      Sum.Ct=Sum.Ct+ if (length(A.Ct.beta.t)==0) {0} else {401*z[i]*t(A.Ct.beta.t)}
    }
    T.Ct=Sum.Ct/401
    B.Ct[t]=t(T.Ct)%*%solve(EV.Ct[[t]])%*%T.Ct
  }
  Quant.Ct[j]<-max(B.Ct)
}
t.stat<-T3(C,c)
p.val<-mean(Quant.Ct>=t.stat)



library(matrixcalc)
 

A.beta<-vector("list",399)
A.beta.curl<-vector("list",399)
for (i in 1:399){
  A.beta.i<-matrix(0,168,5)
  A.beta.i.curl<-matrix(0,168,5)
  for (t in 1:168){
    if (nrow(X[X$id==i & X$Day==t,])!=0){
      R<-as.numeric(sum(is.na(X[X$id==i & X$Day==t,]))==0)
      #R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
      Z.t<-model.matrix(~n2+SEX+RACEW+ISHAK,data=X[X$id==i & X$Day==t,])
      D.beta.t<-Z.t*exp(sum(betas.y1.sm[t,]*Z.t))/((1+exp(sum(betas.y1.sm[t,]*Z.t)))^2)
      Y.bar.t<- new.data$Y1[new.data$id==i & new.data$Day==t]-exp(sum(betas.y1.sm[t,]*Z.t))/(1+exp(sum(betas.y1.sm[t,]*Z.t)))
      V.beta.t<-(1+exp(sum(betas.y1.sm[t,]*Z.t)))^2/exp(sum(betas.y1.sm[t,]*Z.t))
      A<-R*V.beta.t*D.beta.t*Y.bar.t
      A.beta.i[t,]<-if (length(A)==0) {rep(0,5)} else {A}
      A.beta.i.curl[t,]<-solve(H.t.arr[,,t])%*%A.beta.i[t,]
    }
  }
  A.beta[[i]]<-A.beta.i
  A.beta.curl[[i]]<-A.beta.i.curl
}


  



H<-vector("list",168)
X<-subset(new.data,select = c(Y1,id,idn,Day,n2,SEX,RACEW,ISHAK))
D.beta.t.i<-matrix(0,399,5)
Y.bar.t.i<-rep(0,399)
V.beta.t.i<-rep(0,399)
Z.t.i<-matrix(0,399,5)
for (t in 1:168){
  
  H.t<-matrix(0,ncol(beta.y1),ncol(beta.y1))
  mat<-array(0,c(ncol(beta.y1),ncol(beta.y1),399))
  vec1<-rep(0,399)
  for (i in 1:399){
    if (nrow(X[X$idn==i & X$Day==t,])!=0){
      R[i]<-as.numeric(sum(is.na(X[X$idn==i & X$Day==t,]))==0)
      Z.t.i[i,]<-model.matrix(~n2+SEX+RACEW+ISHAK ,data=X[X$idn==i & X$Day==t,])
      D.beta.t.i[i,]<-Z.t.i[i,]*exp(sum(betas.y1.sm[t,]*Z.t.i[i,]))/((1+exp(sum(betas.y1.sm[t,]*Z.t.i[i,])))^2)
      Y.new.i[i]<-new.data$Y1[new.data$idn==i & new.data$Day==t]
      Y.bar.t.i[i]<- new.data$Y1[new.data$idn==i & new.data$Day==t]-exp(sum(betas.y1.sm[t,]*Z.t.i[i,]))/(1+exp(sum(betas.y1.sm[t,]*Z.t.i[i,])))
      V.beta.t.i[i]<-(1+exp(sum(betas.y1.sm[t,]*Z.t.i[i,])))^2/exp(sum(betas.y1.sm[t,]*Z.t.i[i,]))
      mat[,,i]<-D.beta.t[i,]%*%t(D.beta.t[i,])
      vec1[i]<-R[i]*Y.bar.t.i[i]*V.beta.t.i[i]
      H.t.i<-R[i]*V.beta.t.i[i]*D.beta.t.i[i,]%*%t(D.beta.t.i[i,])
      H.t<-H.t+H.t.i
    }
  }  	
  
  H[[t]]<-H.t
}
