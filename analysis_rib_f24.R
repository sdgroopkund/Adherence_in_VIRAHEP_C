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
########################################################################



xmt<-demo[,c(1:7,9:16)]

#xmt<-model.matrix(~SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=demo)


beta<-vector("list",168)
var<-vector("list",168)
beta_rob<-vector("list",168)
var_rob<-vector("list",168)
cl_u_rob<-vector("list",168)
cl_l_rob<-vector("list",168)
cl_u<-vector("list",168)
cl_l<-vector("list",168)

for (i in 1:168){
xvmt<-rib[rib$Day==i,]
xmt1<-merge(xvmt,xmt,by="vhcid")
beta[[i]]<-glm(svr~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol, data=xmt1,family = binomial(link = "logit"),na.action="na.omit")$coefficients
var[[i]]<-vcov(glm(svr~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol, data=xmt1,family = binomial(link = "logit"),na.action="na.exclude"))
cl_l[[i]]<-beta[[i]]-sqrt(diag(var[[i]]))*1.96
cl_u[[i]]<-beta[[i]]+sqrt(diag(var[[i]]))*1.96
}

for (i in 1:168){
xvmt<-rib[rib$Day==i,]
xmt1<-merge(xvmt,xmt,by="vhcid")
beta_rob[[i]]<-glmRob(svr~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol, data=xmt1,family = binomial(link = "logit"),na.action="na.omit")$coefficients
var_rob[[i]]<-vcov(glmRob(svr~n2+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol, data=xmt1,family = binomial(link = "logit"),na.action="na.exclude"))
cl_l_rob[[i]]<-beta[[i]]-sqrt(diag(var[[i]]))*1.96
cl_u_rob[[i]]<-beta[[i]]+sqrt(diag(var[[i]]))*1.96
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
betas.sm<-apply(betas,2,function(x){tri_smoothed(time,x,2)})
beta.sm<-lapply(1:nrow(betas.sm), function(i) betas.sm[i,])

A.beta2<-A.beta
A.beta.curl2<-A.beta.curl

A.beta.sm<-lapply(A.beta2,function(x){apply(x,2,function(x){tri_smoothed(time,x,2)})})
A.beta.curl.sm<-lapply(A.beta.curl2,function(x){apply(x,2,function(x){tri_smoothed(time,x,2)})})

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


plot(time, y1, type='l',ylim=c(-2,3),ylab="Race")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),col="blue",ylab="Race")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),col=2,ylab="Race")


for (i in 1:168){
  y1[i]<-beta.sm[[i]][7]
  y2[i]<-lcl.ci.sm[[i]][7]
  y3[i]<-ucl.ci.sm[[i]][7]
}

plot(time, y1, type='l',ylim=c(-2,4),axes=F,ylab="change in log-odds of SVR for caucasians v non-caucasians",xlab="Days",main="Adherence for Ribavirin (Smoothed)")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,4),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin (Smoothed)")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,4),axes=F,col=2,ylab="change in log-odds for svr",xlab="Days",main="Adherence for Ribavirin (Smoothed)")
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)

plot(time, y1, type='l',ylim=c(-2,3),axes=F,ylab="change in log-odds of SVR for caucasians v non-caucasians",xlab="Days",main="Effect of Race on SVR in Ribavirin analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),axes=F,col="blue",ylab="change in log-odds of SVR for caucasians v non-caucasians",xlab="Days",main="Effect of Race on SVR in Ribavirin analysis")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),axes=F,col=2,ylab="change in log-odds of SVR for caucasians v non-caucasians",xlab="Days",main="Effect of Race on SVR in Ribavirin analysis")
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)
lines(1:168,rep(0,168),type='l')

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
y1[i]<-beta.sm[[i]][7]
y2[i]<-lcl.curl.sm[[i]][7]
y3[i]<-ucl.curl.sm[[i]][7]

}




plot(time, y1, type='l',ylim=c(-2,2),axes=F,ylab="change in log-odds of SVR for males v females",xlab="Days",main="Effect of Gender on SVR in Ribavirin analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,2),axes=F,col="blue",ylab="change in log-odds of SVR for males v females",xlab="Days",main="Effect of Gender on SVR in Ribavirin analysis")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,2),axes=F,col=2,ylab="change in log-odds of SVR for males v females",xlab="Days",main="Effect of Gender on SVR in Ribavirin analysis")
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)
lines(1:168,rep(0,168),type='l')


plot(time, y1, type='l',ylim=c(-1,1),axes=F,ylab="change in log-odds of SVR for ISHAK score",xlab="Days",main="Effect of ISHAK on SVR in Ribavirin analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-1,1),axes=F,col="blue",ylab="change in log-odds of SVR for ISHAK score",xlab="Days",main="Effect of ISHAK on SVR in Ribavirin analysis")
par(new=T)
plot(time, y3, type='l',ylim=c(-1,1),axes=F,col=2,ylab="change in log-odds of SVR for ISHAK score",xlab="Days",main="Effect of ISHAK on SVR in Ribavirin analysis")
Axis(side=1,at=c(0,42,84,126,168))
Axis(side=2)
lines(1:168,rep(0,168),type='l')



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

y1[i]<-beta[[i]][2]
y2[i]<-lcl[[i]][2]
y3[i]<-ucl[[i]][2] 
}
  



plot(time, y1, type='l',ylim=c(-2,4),ylab="change in log-odds for svr",xlab="Days",main="Compliance for Ribavirin (Confidence Bands)")
par(new=T)
plot(time, y4, type='l',ylim=c(-2,4),col="blue",ylab="change in log-odds for svr",xlab="Days")
par(new=T)
plot(time, y5, type='l',ylim=c(-2,4),col="red",ylab="change in log-odds for svr",xlab="Days")




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
for (h in 4:10){
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




