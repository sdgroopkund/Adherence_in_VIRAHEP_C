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


######################
## removing dropouts #
######################

rib<-rib[rib$dropout==0,]
rib<-rib[rib$dropout2==0,]

#######################

L<-data.frame(cbind(vhcid=demo$vhcid,id=seq(1,length(demo$vhcid))))
rib_demo<-merge(rib,merge(L,demo,by='vhcid'),by='vhcid')


rib_demo1<-rib_demo[rib_demo$Day<=168,]
rib_demo21<-rib_demo[rib_demo$Day>168,]
rib_demo2<-rib_demo21[rib_demo21$response24=='Responder',]
rib_demo3<-rib_demo[rib_demo$response24=='Responder',]



######################
## removing dropouts #
######################
 peg<-peg[-8449,]
 peg2<-peg[is.na(peg$do_not_include)==1,]
 peg3<-peg2[is.na(peg2$do_not_include2)==1,]

#######################





L<-data.frame(cbind(vhcid=demo$vhcid,id=seq(1,length(demo$vhcid))))
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


################### vload baseline #######################################

vload.base<-vload[vload$time==0,]
vload.base.keep<-subset(vload.base,select=c(vhcid,vload_itt))
vload_b<-vload$vload_itt[which(vload$time==0)]

un.vb<-unique(vload.base$vhcid)
un.v<-unique(vload$vhcid)

demo<-merge(demo,vload.base.keep,by='vhcid',all.x=T)
########################################################################



xmt<-demo[,c(1:7,9:16)]
#xmt<-model.matrix(~SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=demo)


beta<-vector("list",24)
var<-vector("list",24)
cl_u<-vector("list",24)
cl_l<-vector("list",24)

for (i in 1:24){
xvmt<-rib_peg[rib_peg$week==i,]
xmt1<-merge(xvmt,xmt,by="vhcid")
beta[[i]]<-glm(svr~n+n2.1+n2.2+n2.3+n2.4+n2.5+n2.6+n2.7+SEX+RACEW+ISHAK, data=xmt1,family = binomial(link = "logit"),na.action="na.omit")$coefficients
var[[i]]<-vcov(glm(svr~n+n2.1+n2.2+n2.3+n2.4+n2.5+n2.6+n2.7+SEX+RACEW+ISHAK, data=xmt1,family = binomial(link = "logit"),na.action="na.exclude"))
cl_l[[i]]<-beta[[i]]-sqrt(diag(var[[i]]))*1.96
cl_u[[i]]<-beta[[i]]+sqrt(diag(var[[i]]))*1.96
}

time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)

for (i in 1:24){
y1[i]<-beta[[i]][2]
y2[i]<-cl_l[[i]][2]
y3[i]<-cl_u[[i]][2]
 }

plot(time, y1, type='l',ylim=c(-3,8),ylab="Compliance")
par(new=T)
plot(time, y2, type='l',ylim=c(-3,8),col="blue",ylab="Compliance")
par(new=T)
plot(time, y3, type='l',ylim=c(-3,8),col=2,ylab="Compliance")




############################ Analysis as mean of effects ###############################



xmt<-demo[,c(1:7,9:16)]
#xmt<-model.matrix(~SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=demo)


beta.m<-vector("list",24)
var.m<-vector("list",24)
cl_u<-vector("list",24)
cl_l<-vector("list",24)

for (i in 1:24){
xvmt<-rib_peg[rib_peg$week==i,]
xmt1<-merge(xvmt,xmt,by="vhcid")
beta.m[[i]]<-glm(svr~n*n2+SEX+RACEW+ISHAK, data=xmt1,family = binomial(link = "logit"),na.action="na.omit")$coefficients
var.m[[i]]<-vcov(glm(svr~n*n2+SEX+RACEW+ISHAK, data=xmt1,family = binomial(link = "logit"),na.action="na.exclude"))
cl_l[[i]]<-beta.m[[i]]-sqrt(diag(var.m[[i]]))*1.96
cl_u[[i]]<-beta.m[[i]]+sqrt(diag(var.m[[i]]))*1.96
}

time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)
yy<-rep(0,24)

for (i in 1:24){
y1[i]<-beta.m[[i]][3]
y2[i]<-beta.m[[i]][2]
yy[i]<-beta.m[[i]][7]
}

plot(time, y1, type='l',lwd=2,ylim=c(-40,40),axes='F',ylab="change in log-odds for svr",xlab="Weeks",main="Adherence to Ribavirin/Peginterferon")
par(new=T)
plot(time, y2, type='l',lty=5,lwd=2,ylim=c(-40,40),axes='F',ylab="change in log-odds for svr",xlab='Weeks')
par(new=T)
plot(time, yy, type='l',lwd=2,lty=6,ylim=c(-40,40),axes='F',ylab="change in log-odds for svr",xlab='Weeks')
legend(15,40 ,lty=c(1,3,5),lwd=c(2,2,2), c('Ribavirin','Peginterferon','Interaction') )
Axis(side=2)
Axis(side=1,at=c(1,6,12,18,24))


par(new=T)
plot(time, y2, type='l',ylim=c(),col="blue",ylab="Compliance_peg")
par(new=T)
plot(time, y3, type='l',ylim=c(-50,80),col=2,ylab="Compliance_peg")
lines(rep(0,24))

##########################################################################################

xvmt<-rib_peg[rib_peg$week==i,]
xmt1<-merge(xvmt,xmt,by="vhcid")
xmt2<-subset(xmt1,select=c("n","n2","svr"))
int<-xmt2$n*xmt2$n2
xmt2<-cbind(xmt2,"n-n2"=int)


library(e1071)
library(rpart)

index <- 1:nrow(xmt2)
testindex <- sample(index, trunc(length(index)/3))
 testset <- xmt2[testindex, ]
 trainset <- xmt2[-testindex, ]

svm.model <- svm(svr ~ ., data = trainset, cost = 1,gamma=0.1,type="C",kernel="linear")
svm.pred <- predict(svm.model, testset[, -3],type="class")

table(pred = svm.pred, true = testset[, 3])

rpart.model <- rpart(svr ~ n+n2+n-n2, data = trainset)
rpart.pred <- predict(rpart.model, testset[, -3])


table(pred = rpart.pred, true = testset[, 3])
###################################################################################


data2 <- seq(1,10)

classes2 <- c('b','b','b','a','a','a','a','b','b','b')


model2 <- svm(data2,classes2,type='C',kernel='linear')

predict(model2,data2)
table(predict(model2,data2), classes2)

###########################################################################################

H.m<-vector("list",24)
X<-subset(peg_rib_demo,select = c(id,week,n,n2,n2.1,n2.2,n2.3,n2.4,n2.5,n2.6,n2.7,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))
for (t in 1:24){
	H.t<-matrix(0,length(beta.m[[1]]),length(beta.m[[1]]))
	D<-vector("list",401)
	for (i in 1:401){
		if (nrow(X[X$id==i & X$week==t,])!=0){
			R<-as.numeric(sum(is.na(X[X$id==i & X$week==t,]))==0)
			#R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
			Z.t<-model.matrix(~n+n2+SEX+RACEW+ISHAK ,data=X[X$id==i & X$week==t,])
			D.beta.t<-Z.t*exp(sum(beta.m[[t]]*Z.t))/((1+exp(sum(beta.m[[t]]*Z.t)))^2)
			#Y.bar.t<- peg_rib_demo$svr[peg_rib_demo$id==i & peg_rib_demo$week==t]-exp(sum(beta.m[[t]]*Z.t))/(1+exp(sum(beta.m[[t]]*Z.t)))
			V.beta.t<-(1+exp(sum(beta.m[[t]]*Z.t)))^2/exp(sum(beta.m[[t]]*Z.t))
			H.t.i<-R*t(D.beta.t)%*%V.beta.t%*%D.beta.t
			H.t<-H.t+H.t.i
			D[[i]]<-H.t.i
		}
	}		

H.m[[t]]<-H.t
}




A.beta.m<-vector("list",24)
A.beta.curl.m<-vector("list",24)
X<-subset(peg_rib_demo,select = c(id,week,n,n2,n2.1,n2.2,n2.3,n2.4,n2.5,n2.6,n2.7,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))

for (i in 1:401)
	{
	A.beta.i<-matrix(0,24,length(beta.m[[1]]))
	A.beta.i.curl<-matrix(0,24,length(beta.m[[1]]))
	for (t in 1:24){
		if (nrow(X[X$id==i & X$week==t,])!=0){
			R<-as.numeric(sum(is.na(X[X$id==i & X$week==t,]))==0)
			#R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
			Z.t<-model.matrix(~n+n2+SEX+RACEW+ISHAK ,data=X[X$id==i & X$week==t,])
			D.beta.t<-Z.t*exp(sum(beta.m[[t]]*Z.t))/((1+exp(sum(beta.m[[t]]*Z.t)))^2)
			Y.bar.t<- peg_rib_demo$svr[peg_rib_demo$id==i & peg_rib_demo$week==t]-exp(sum(beta.m[[t]]*Z.t))/(1+exp(sum(beta.m[[t]]*Z.t)))
			V.beta.t<-(1+exp(sum(beta.m[[t]]*Z.t)))^2/exp(sum(beta.m[[t]]*Z.t))
			A<-R*V.beta.t*D.beta.t*Y.bar.t
			A.beta.i[t,]<-if (length(A)==0) {rep(0,length(beta.m[[1]]))} else {A}
			A.beta.i.curl[t,]<-solve(H.m[[t]])%*%A.beta.i[t,]
		}
	}
	A.beta.m[[i]]<-A.beta.i
	A.beta.curl.m[[i]]<-A.beta.i.curl
}



Var.m<-function(t){
G<-matrix(0,length(beta.m[[1]]),length(beta.m[[1]]))
for (i in 1:401){
A.beta.t<-t(A.beta.m[[i]][t,])
G.i<-t(A.beta.t)%*%A.beta.t
G<-G+G.i
}
return(solve(H.m[[t]])%*%G%*%solve(H.m[[t]]))
}

ucl.ci<-vector("list",24)
lcl.ci<-vector("list",24)
for (t in 1:24){
lcl.ci[[t]]<-beta.m[[t]]-sqrt(diag(Var.m(t)))*1.96
ucl.ci[[t]]<-beta.m[[t]]+sqrt(diag(Var.m(t)))*1.96
}




time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)

for (i in 1:24){
y1[i]<-beta.m[[i]][3]
y2[i]<-lcl.ci[[i]][3]
y3[i]<-ucl.ci[[i]][3]
 }

plot(time, y1, type='l',ylim=c(-3,8),ylab="Compliance_rib")
par(new=T)
plot(time, y2, type='l',ylim=c(-3,8),col="blue",ylab="Compliance_rib")
par(new=T)
plot(time, y3, type='l',ylim=c(-3,8),col=2,ylab="Compliance_rib")








#################################################################################



H<-vector("list",24)
X<-subset(peg_rib_demo,select = c(id,week,n,n2,n2.1,n2.2,n2.3,n2.4,n2.5,n2.6,n2.7,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))
for (t in 1:24){
	H.t<-matrix(0,length(beta[[1]]),length(beta[[1]]))
	D<-vector("list",401)
	for (i in 1:401){
		if (nrow(X[X$id==i & X$week==t,])!=0){
			R<-as.numeric(sum(is.na(X[X$id==i & X$week==t,]))==0)
			#R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
			Z.t<-model.matrix(~n+n2.1+n2.2+n2.3+n2.4+n2.5+n2.6+n2.7+SEX+RACEW+ISHAK ,data=X[X$id==i & X$week==t,])
			D.beta.t<-Z.t*exp(sum(beta[[t]]*Z.t))/((1+exp(sum(beta[[t]]*Z.t)))^2)
			#Y.bar.t<- peg_rib_demo$svr[peg_rib_demo$id==i & peg_rib_demo$week==t]-exp(sum(beta[[t]]*Z.t))/(1+exp(sum(beta[[t]]*Z.t)))
			V.beta.t<-(1+exp(sum(beta[[t]]*Z.t)))^2/exp(sum(beta[[t]]*Z.t))
			H.t.i<-R*t(D.beta.t)%*%V.beta.t%*%D.beta.t
			H.t<-H.t+H.t.i
			D[[i]]<-H.t.i
		}
	}		

H[[t]]<-H.t
}




A.beta<-vector("list",24)
A.beta.curl<-vector("list",24)
X<-subset(peg_rib_demo,select = c(id,week,n,n2.1,n2.2,n2.3,n2.4,n2.5,n2.6,n2.7,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))

for (i in 1:401)
	{
	A.beta.i<-matrix(0,24,length(beta[[1]]))
	A.beta.i.curl<-matrix(0,24,length(beta[[1]]))
	for (t in 1:24){
		if (nrow(X[X$id==i & X$week==t,])!=0){
			R<-as.numeric(sum(is.na(X[X$id==i & X$week==t,]))==0)
			#R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
			Z.t<-model.matrix(~n+n2.1+n2.2+n2.3+n2.4+n2.5+n2.6+n2.7+SEX+RACEW+ISHAK ,data=X[X$id==i & X$week==t,])
			D.beta.t<-Z.t*exp(sum(beta[[t]]*Z.t))/((1+exp(sum(beta[[t]]*Z.t)))^2)
			Y.bar.t<- peg_rib_demo$svr[peg_rib_demo$id==i & peg_rib_demo$week==t]-exp(sum(beta[[t]]*Z.t))/(1+exp(sum(beta[[t]]*Z.t)))
			V.beta.t<-(1+exp(sum(beta[[t]]*Z.t)))^2/exp(sum(beta[[t]]*Z.t))
			A<-R*V.beta.t*D.beta.t*Y.bar.t
			A.beta.i[t,]<-if (length(A)==0) {rep(0,length(beta[[1]]))} else {A}
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

ucl.ci<-vector("list",24)
lcl.ci<-vector("list",24)
for (t in 1:24){
lcl.ci[[t]]<-beta[[t]]-sqrt(diag(Var(t)))*1.96
ucl.ci[[t]]<-beta[[t]]+sqrt(diag(Var(t)))*1.96
}




time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)

for (i in 1:24){
y1[i]<-beta[[i]][2]
y2[i]<-lcl.ci[[i]][2]
y3[i]<-ucl.ci[[i]][2]
 }

plot(time, y1, type='l',ylim=c(-3,8),ylab="Compliance_rib_D7")
par(new=T)
plot(time, y2, type='l',ylim=c(-3,8),col="blue",ylab="Compliance_rib_D7")
par(new=T)
plot(time, y3, type='l',ylim=c(-3,8),col=2,ylab="Compliance_rib_D7")







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



betas<-matrix(unlist(beta),24,length(beta[[1]]),byrow=T)
betas.sm<-apply(betas,2,function(x){tri_smoothed(time,x,3)})
beta.sm<-lapply(1:nrow(betas.sm), function(i) betas.sm[i,])

A.beta2<-A.beta
A.beta.curl2<-A.beta.curl

A.beta.sm<-lapply(A.beta2,function(x){apply(x,2,function(x){tri_smoothed(time,x,3)})})
A.beta.curl.sm<-lapply(A.beta.curl2,function(x){apply(x,2,function(x){tri_smoothed(time,x,3)})})



######################################################################################################



Var.curl.sm<-function(t,A){
G<-matrix(0,length(beta[[1]]),length(beta[[1]]))

for (i in 1:401){
A.beta.curl.sm.t<-t(A[[i]][t,])
G.i<-t(A.beta.curl.sm.t)%*%A.beta.curl.sm.t
G<-G+G.i
}
return(G)
}




EV.curl.sm<-vector("list",24)
for (t in 1:24){
EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm))))

}


Start<-Sys.time()


B.curl.sm1=matrix(0,length(beta[[1]]),24)
Quant<-matrix(0,length(beta[[1]]),5000)
for (j in 1:5000){
z<-rnorm(401,0,1)
for (t in 1:24)
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


ucl.curl.sm<-vector("list",24)
lcl.curl.sm<-vector("list",24)


for (t in 1:24){

ucl.curl.sm[[t]]<-beta.sm[[t]]+EV.curl.sm[[t]]%*%Quant.alpha/sqrt(401)
lcl.curl.sm[[t]]<-beta.sm[[t]]-EV.curl.sm[[t]]%*%Quant.alpha/sqrt(401)
}

time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)
y4<-rep(0,24)
y5<-rep(0,24)


for (i in 1:24){
y1[i]<-beta.sm[[i]][2]
y2[i]<-lcl.curl.sm[[i]][2]
y3[i]<-ucl.curl.sm[[i]][2]
}




plot(time, y1, type='l',ylim=c(-3,8),,ylab="Compliance")
par(new=T)
plot(time, y2, type='l',ylim=c(-3,8),,col="red",ylab="Compliance")
par(new=T)
plot(time, y3, type='l',ylim=c(-3,8),col="red",ylab="Compliance")









EV<-vector("list",24)
for (t in 1:24){
EV[[t]]<-sqrt(diag(diag(Var(t))))
}


B.quant=matrix(0,length(beta[[1]]),24)
Quant<-matrix(0,length(beta[[1]]),5000)
for (j in 1:5000){
z<-rnorm(401,0,1)
for (t in 1:24)
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

ucl<-vector("list",24)
lcl<-vector("list",24)

for (t in 1:24){

ucl[[t]]<-beta[[t]]+EV[[t]]%*%Quant.alpha.unsm/sqrt(401)
lcl[[t]]<-beta[[t]]-EV[[t]]%*%Quant.alpha.unsm/sqrt(401)
}






time<-seq(1,24)
y<-rep(0,24)
y1<-rep(0,24)
y12<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)
y4<-rep(0,24)
y5<-rep(0,24)


for (i in 1:24){
#y[i]<-beta[[i]][2]
y1[i]<-beta[[i]][2]
#y12[i]<-beta.sm2[[i]][2]
#y2[i]<-lcl.curl.sm[[i]][2]
#y3[i]<-ucl.curl.sm[[i]][2]
y4[i]<-lcl[[i]][2]
y5[i]<-ucl[[i]][2] 
}



#plot(time, y, type='l',ylim=c(-3,8),ylab="Compliance")
#par(new=T)
plot(time, y1, type='l',ylim=c(-3,8),ylab="Compliance",col='red')
par(new=T)
#plot(time, y12, type='l',ylim=c(-2,3),ylab="Compliance",col='blue')
#par(new=T)
#plot(time, y3, type='l',ylim=c(-2,3),col="red",ylab="Compliance")
#par(new=T)
plot(time, y4, type='l',ylim=c(-3,8),col="blue",ylab="Compliance")
par(new=T)
plot(time, y5, type='l',ylim=c(-3,8),col="blue",ylab="Compliance")























#T2

T2<-function(C,c,W){

int<-matrix(0,length(c[[1]]),24)

for (t in 1:24){

C.t<-C[[t]]
beta.t<-beta[[t]]
c.t<-c[[t]]
W.t<-W[[t]]

int[,t]<-(C.t%*%beta.t - c.t)*W.t
}

int1<- rowSums(int) - (int[,1]+int[,24])/2
Z=0
for (t in 1:24){
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
	Z<-Z+if(t==1 || t==24){Int.t/2}else{Int.t}
}
AA<-matrix(0,nrow(C[[1]]),nrow(C[[1]]))
for (i in 1:401){
	AA<-AA+Z[,i]%*%t(Z[,i])
	}
sigma<-(1/(401^2))*AA
#sigma<-(1/(401^2))*matrix(rowSums(rbind(apply(Z,2,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))


t2<-t(int1)%*%solve(sigma)%*%int1
pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
 return(pval)
}


contrast<-rep(0,length(beta[[1]]))
contrast[2]<-1
C<-vector(mode="list",24)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",24)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",24)
W<-lapply(W,function(x){x<-1})



c.mat<-matrix(0,6,length(beta[[1]]))
cc<-cbind(c(4,2,0,0,0,0),c(4,-1,1,0,0,0),c(4,-1,-1,0,0,0),c(-3,0,0,3,0,0),c(-3,0,0,-1,2,0),c(-3,0,0,-1,-1,1),c(-3,0,0,-1,-1,-1))
c.mat[,3:9]<-cc
c.mat<-apply(c.mat,1,function(x){x/sqrt(sum(x^2))})

contrast<-t(c.mat)
C<-vector(mode="list",24)
C<-lapply(C,function(x){x<-contrast})
c<-vector(mode="list",24)
c<-lapply(c,function(x){x<-rep(0,6)})






T2(C,c,W)







T2.sm<-function(C,c,W){

int<-matrix(0,length(c[[1]]),24)
for (t in 1:24){
C.t<-C[[t]]
beta.sm.t<-beta.sm[[t]]
beta.t<-beta[[t]]
c.t<-c[[t]]
W.t<-W[[t]]
int[,t]<-(C.t%*%beta.sm.t - c.t)*W.t
}
int1<- rowSums(int) - (int[,1]+int[,24])/2
Z.sm=0
for (t in 1:24){
	C.t<-C[[t]]
	c.t<-c[[t]]
	W.t<-W[[t]]
	A.t<-matrix(0,length(beta[[1]]),401)
	for (i in 1:401){
		A.t[,i]<-A.beta.curl[[i]][t,]
	}
	Int.t.sm<-C.t%*%A.t*W.t
	Z.sm<-Z.sm+if(t==1 || t==24){Int.t.sm/2}else{Int.t.sm}
}
sigma.sm<-matrix(rowSums(rbind(apply(Z.sm,2,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))
t2<-t(int1)%*%solve(sigma.sm)%*%int1
pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
return(pval)
}


contrast<-rep(0,length(beta[[1]]))
contrast[2]<-1
C<-vector(mode="list",24)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",24)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",24)
W<-lapply(W,function(x){x<-1})

T2.sm(C,c,W)












T3<-function(C,c){
	Stat<-rep(0,24)
	for (t in 1:24){
		C.t<-C[[t]]
		beta.t<-beta[[t]]
		c.t<-c[[t]]
		T.t<-(C.t%*%beta.t - c.t)		
		Cov.t<-C.t%*%Var(t)%*%t(C.t)
		Stat[t]<-t(T.t)%*%solve(Cov.t)%*%T.t
	}
	T3<-max(Stat)
return(T3)
}




c.mat<-matrix(0,6,length(beta[[1]]))
cc<-cbind(c(4,2,0,0,0,0),c(4,-1,1,0,0,0),c(4,-1,-1,0,0,0),c(-3,0,0,3,0,0),c(-3,0,0,-1,2,0),c(-3,0,0,-1,-1,1),c(-3,0,0,-1,-1,-1))
c.mat[,3:9]<-cc
c.mat<-apply(c.mat,1,function(x){x/sqrt(sum(x^2))})

contrast<-t(c.mat)
C<-vector(mode="list",24)
C<-lapply(C,function(x){x<-contrast})
c<-vector(mode="list",24)
c<-lapply(c,function(x){x<-rep(0,6)})


A.Ct.beta<-lapply(A.beta.curl,function(x){apply(x,1,function(z){C[[t]]%*%z})})


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

########################### Testing for equality of effects ######################################



B.Ct=rep(0,24)
Quant.Ct<-rep(0,5000)
for (j in 1:5000){
z<-rnorm(401,0,1)
for (t in 1:24)
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


###################################################################################################





B.Ct=matrix(0,length(C[[1]]%*%beta[[1]]),24)

Quant.Ct<-matrix(0,length(C[[1]]%*%beta[[1]]),5000)
for (j in 1:5000){
z<-rnorm(401,0,1)
for (t in 1:24)
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


ucl.Ct.sm<-vector("list",24)
lcl.Ct.sm<-vector("list",24)


for (t in 1:24){

ucl.Ct.sm[[t]]<-(C[[t]]%*%beta.sm[[t]]-c[[t]])+EV.Ct[[t]]%*%Quant.Ct.alpha[1]/sqrt(401)
lcl.Ct.sm[[t]]<-(C[[t]]%*%beta.sm[[t]]-c[[t]])-EV.Ct[[t]]%*%Quant.Ct.alpha[1]/sqrt(401)
}




time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)
y4<-rep(0,24)
y5<-rep(0,24)


for (i in 1:24){
y1[i]<-beta.sm[[i]][2]
y2[i]<-lcl.curl.sm[[i]][2]
y3[i]<-ucl.curl.sm[[i]][2]
y4[i]<-lcl.Ct.sm[[i]]
y5[i]<-ucl.Ct.sm[[i]] 
}




plot(time, y1, type='l',ylim=c(-3,8),ylab="Compliance")
par(new=T)
plot(time, y2, type='l',ylim=c(-3,8),col="red",ylab="Compliance")
par(new=T)
plot(time, y3, type='l',ylim=c(-3,8),col="red",ylab="Compliance")
par(new=T)
plot(time, y4, type='l',ylim=c(-3,8),col="blue",ylab="Compliance")
par(new=T)
plot(time, y5, type='l',ylim=c(-3,8),col="blue",ylab="Compliance")

















# Cross Validation



betas<-matrix(unlist(beta),24,length(beta[[1]]),byrow=T)


Mse<-matrix(0,100,24)
V<-matrix(0,100,24)
B<-matrix(0,100,24)
for (h in 1:100){
	betas.sm.h<-apply(betas,2,function(x){tri_smoothed(time,x,h)})
	beta.h<-lapply(1:nrow(betas.sm.h), function(i) betas.sm.h[i,])
	A.h<-lapply(A.beta.curl2,function(x){apply(x,2,function(x){tri_smoothed(time,x,h)})})
	for (t in 1:24){
		var.beta.h.t<-sum(diag(Var.curl.sm(t,A.h)))
		V[h,t]<-var.beta.h.t
		sq_bias.h.t<-sum((beta[[t]]-beta.h[[t]])^2)
		B[h,t]<-sq_bias.h.t
		Mse[h,t]<-var.beta.h.t+sq_bias.h.t
	}
}



betas<-matrix(unlist(beta),24,length(beta[[1]]),byrow=T)
betas.sm2<-apply(betas,2,function(x){tri_smoothed(time,x,2)})
beta.sm2<-lapply(1:nrow(betas.sm2), function(i) betas.sm2[i,])






i=3
xvmt<-rib_peg[rib_peg$week==i,]
xmt1<-merge(xvmt,xmt,by="vhcid")
xmt1$svr[intersect(which(xmt1$n==0),which(xmt1$n2<1))]
mean(xmt1$svr[which(xmt1$n==0)])
xmt1$svr[which(xmt1$n2<1)]







