demo<-read.table('demo.txt',header=T,sep='\t')
peg<-read.table('pegifn.txt',header=T,sep='\t')
rib<-read.table('rib.txt',header=T,sep='\t')
vload<-read.table('vload.txt',header=T,sep='\t')

demo[demo==""]<-NA
demo$ISHAK[demo$ISHAK=="B"]<-NA
demo$ISHAK <- as.numeric(as.character(demo$ISHAK))
write.table(demo, file = "demo1.txt",col.names = TRUE,sep='\t')
demo<-read.table('demo1.txt',header=T,sep='\t')

peg<-peg[-8449,]


######################
## removing dropouts #
######################

 peg2<-peg[is.na(peg$do_not_include)==1,]
 peg3<-peg2[is.na(peg2$do_not_include2)==1,]

#######################

L<-data.frame(cbind(vhcid=demo$vhcid,id=seq(1,length(demo$vhcid))))
peg_demo<-merge(peg3,merge(L,demo,by='vhcid'),by='vhcid')


peg_demo1<-peg_demo[peg_demo$week<=24,]

peg_demo21<-peg_demo[peg_demo$week>24,]
peg_demo2<-peg_demo21[peg_demo21$response24=='Responder',]

peg_demo3<-peg_demo[peg_demo$response24=='Responder',]


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
xvmt<-peg3[peg3$week==i,]
xmt1<-merge(xvmt,xmt,by="vhcid")
beta[[i]]<-glm(svr~n+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol, data=xmt1,family = binomial(link = "logit"),na.action="na.omit")$coefficients
var[[i]]<-vcov(glm(svr~n+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol, data=xmt1,family = binomial(link = "logit"),na.action="na.exclude"))
}


####################################################################################################################






H<-vector("list",24)
X<-subset(peg_demo1,select = c(id,week,n,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))
for (t in 1:24){
	H.t<-matrix(0,length(beta[[1]]),length(beta[[1]]))
	#D<-vector("list",401)
	for (i in 1:401){
		if (nrow(X[X$id==i & X$week==t,])!=0){
			R<-as.numeric(sum(is.na(X[X$id==i & X$week==t,]))==0)
			#R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
			Z.t<-model.matrix(~n+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=X[X$id==i & X$week==t,])
			D.beta.t<-Z.t*exp(sum(beta[[t]]*Z.t))/((1+exp(sum(beta[[t]]*Z.t)))^2)
			Y.bar.t<- peg_demo1$svr[peg_demo1$id==i & peg_demo1$week==t]-exp(sum(beta[[t]]*Z.t))/(1+exp(sum(beta[[t]]*Z.t)))
			V.beta.t<-(1+exp(sum(beta[[t]]*Z.t)))^2/exp(sum(beta[[t]]*Z.t))
			H.t.i<-R[t]*t(D.beta.t)%*%V.beta.t%*%D.beta.t
			H.t<-H.t+H.t.i
			#D[[i]]<-H.t.i
		}
	}		

H[[t]]<-H.t
}




A.beta<-vector("list",24)
A.beta.curl<-vector("list",24)
X<-subset(peg_demo1,select = c(id,week,n,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))

for (i in 1:401)
	{
	A.beta.i<-matrix(0,24,17)
	A.beta.i.curl<-matrix(0,24,17)
	for (t in 1:24){
		if (nrow(X[X$id==i & X$week==t,])!=0){
			R<-as.numeric(sum(is.na(X[X$id==i & X$week==t,]))==0)
			#R<-as.numeric(apply(X[X$id==i,],1,function(x){sum(is.na(x))==0}))
			Z.t<-model.matrix(~n+SEX+RACEW+MXAD+age+ISHAK+infect+education+insurance+employ+marital+alcohol ,data=X[X$id==i & X$week==t,])
			D.beta.t<-Z.t*exp(sum(beta[[t]]*Z.t))/((1+exp(sum(beta[[t]]*Z.t)))^2)
			Y.bar.t<- peg_demo1$svr[peg_demo1$id==i & peg_demo1$week==t]-exp(sum(beta[[t]]*Z.t))/(1+exp(sum(beta[[t]]*Z.t)))
			V.beta.t<-(1+exp(sum(beta[[t]]*Z.t)))^2/exp(sum(beta[[t]]*Z.t))
			A<-R[t]*V.beta.t*D.beta.t*Y.bar.t
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


plot(time, y1, type='l',ylim=c(-2,8),axes=F,ylab="change in log-odds for svr",xlab="Weeks",main="Adherence to Peginterferon")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,8),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Weeks")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,8),axes=F,col=2,ylab="change in log-odds for svr",xlab="Weeks")
lines(rep(0,24))
Axis(side=1,at=c(1,6,12,18,24))
Axis(side=2)


plot(time, y1, type='l',ylim=c(-2,2),axes=F,ylab="change in log-odds for males v females on svr",xlab="Days",main="Effect of Gender on SVR in Ribavirin Analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,2),axes=F,col="blue",ylab="change in log-odds for males v females on svr",xlab="Days")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,2),axes=F,col=2,ylab="change in log-odds for males v females on svr",xlab="Days")
lines(rep(0,168))
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



betas<-matrix(unlist(beta),24,17,byrow=T)
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

time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)
y4<-rep(0,24)
y5<-rep(0,24)


for (i in 1:24){
y1[i]<-beta[[i]][2]
y2[i]<-lcl[[i]][2]
y3[i]<-ucl[[i]][2]
y4[i]<-cl_l[[i]][1]
y5[i]<-cl_u[[i]][1] 
}


plot(time, y1, type='l',ylim=c(-2,8),axes=F,ylab="change in log-odds for svr",xlab="Weeks",main="Adherence to Peginterferon")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,8),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Weeks")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,8),axes=F,col=2,ylab="change in log-odds for svr",xlab="Weeks")

Axis(side=1,at=c(1,6,12,18,24))
Axis(side=2)



betas<-matrix(unlist(beta),24,17,byrow=T)
betas.sm<-apply(betas,2,function(x){tri_smoothed(time,x,58)})
beta.sm<-lapply(1:nrow(betas.sm), function(i) betas.sm[i,])

A.beta2<-A.beta
A.beta.curl2<-A.beta.curl

A.beta.sm<-lapply(A.beta2,function(x){apply(x,2,function(x){tri_smoothed(time,x,58)})})
A.beta.curl.sm<-lapply(A.beta.curl2,function(x){apply(x,2,function(x){tri_smoothed(time,x,58)})})

############################################################################################################

Var.sm<-function(t,A){
G<-matrix(0,length(beta[[1]]),length(beta[[1]]))
for (i in 1:401){
A.beta.sm.t<-t(A[[i]][t,])
G.i<-t(A.beta.sm.t)%*%A.beta.sm.t
}
return(solve(H[[t]])%*%G%*%solve(H[[t]]))
}

############################################################################################################



Var.curl.sm<-function(t,A){
G<-matrix(0,length(beta[[1]]),length(beta[[1]]))

for (i in 1:401){
A.beta.curl.sm.t<-t(A[[i]][t,])
G.i<-t(A.beta.curl.sm.t)%*%A.beta.curl.sm.t
G<-G+G.i
}
return(G)
}



betas<-matrix(unlist(beta),24,17,byrow=T)
betas.sm<-apply(betas,2,function(x){tri_smoothed(time,x,5)})
beta.sm<-lapply(1:nrow(betas.sm), function(i) betas.sm[i,])

A.beta2<-A.beta
A.beta.curl2<-A.beta.curl

A.beta.sm<-lapply(A.beta2,function(x){apply(x,2,function(x){tri_smoothed(time,x,5)})})
A.beta.curl.sm<-lapply(A.beta.curl2,function(x){apply(x,2,function(x){tri_smoothed(time,x,5)})})



ucl.ci.sm<-vector("list",24)
lcl.ci.sm<-vector("list",24)
for (t in 1:24){
  lcl.ci.sm[[t]]<-beta.sm[[t]]-sqrt(diag(Var.curl.sm(t,A.beta.curl.sm)))*1.96
  ucl.ci.sm[[t]]<-beta.sm[[t]]+sqrt(diag(Var.curl.sm(t,A.beta.curl.sm)))*1.96
}




time<-seq(1,24)
y1<-rep(0,24)
y2<-rep(0,24)
y3<-rep(0,24)

for (i in 1:24){
  y1[i]<-beta.sm[[i]][2]
  y2[i]<-lcl.ci.sm[[i]][2]
  y3[i]<-ucl.ci.sm[[i]][2]
}



plot(time, y1, type='l',ylim=c(-2,8),axes=F,ylab="change in log-odds for svr",xlab="Weeks",main="Adherence to Peginterferon (Hypothesis Test)")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,8),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Weeks")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,8),axes=F,col=2,ylab="change in log-odds for svr",xlab="Weeks")
lines(rep(0,24))
Axis(side=1,at=c(1,6,12,18,24))
Axis(side=2)



EV.curl.sm<-vector("list",24)
for (t in 1:24){
EV.curl.sm[[t]]<-sqrt(diag(diag(Var.curl.sm(t,A.beta.curl.sm))))

}


Start<-Sys.time()


B.curl.sm1=matrix(0,length(beta[[1]]),24)
Quant<-matrix(0,17,5000)
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


plot(time*7, y1, type='l',ylim=c(-1,1),axes=F,ylab="change in log-odds of SVR for ISHAK score",xlab="Days",main="Effect of ISHAK on SVR in Peginterferon Analysis")
par(new=T)
plot(time*7, y2, type='l',ylim=c(-1,1),axes=F,col="blue",ylab="change in log-odds of SVR for ISHAK score",xlab="Days")
par(new=T)
plot(time*7, y3, type='l',ylim=c(-1,1),axes=F,col=2,ylab="change in log-odds of SVR for ISHAK score",xlab="Days")
Axis(side=1,at=c(7,42,84,126,168))
Axis(side=2)
lines(time*7,rep(0,24))


plot(time*7, y1, type='l',ylim=c(-2,3),axes=F,ylab="change in log-odds of SVR for caucasians v non-caucasians",xlab="Days",main="Effect of Race on SVR in Peginterferon Analysis")
par(new=T)
plot(time*7, y2, type='l',ylim=c(-2,3),axes=F,col="blue",ylab="change in log-odds of SVR for caucasians v non-caucasians",xlab="Days")
par(new=T)
plot(time*7, y3, type='l',ylim=c(-2,3),axes=F,col=2,ylab="change in log-odds of SVR for caucasians v non-caucasians",xlab="Days")
Axis(side=1,at=c(7,42,84,126,168))
Axis(side=2)
lines(time*7,rep(0,24))


plot(time*7, y1, type='l',ylim=c(-2,8),axes=F,lwd=2,ylab="change in log-odds for svr",xlab="Days",main="Adherence to Peginterferon (Smoothed)")
par(new=T)
plot(time*7, y2, type='l',ylim=c(-2,8),axes=F,lwd=2,lty=6,ylab="change in log-odds for svr",xlab="Days")
par(new=T)
plot(time*7, y3, type='l',ylim=c(-2,8),axes=F,lwd=2,lty=6,ylab="change in log-odds for svr",xlab="Days")

Axis(side=1,at=c(7,42,84,126,168))
Axis(side=2)



plot(time*7, y1, type='l',ylim=c(-2,2),axes=F,ylab="change in log-odds of SVR for males v females",xlab="Days",main="Effect of Gender on SVR in Peginterferon Analysis")
par(new=T)
plot(time*7, y2, type='l',ylim=c(-2,2),axes=F,col="blue",ylab="change in log-odds of SVR for males v females",xlab="Days")
par(new=T)
plot(time*7, y3, type='l',ylim=c(-2,2),axes=F,col=2,ylab="change in log-odds of SVR for males v females",xlab="Days")
Axis(side=1,at=c(7,42,84,126,168))
Axis(side=2)
lines(time*7,rep(0,24))

plot(time, y1, type='l',ylim=c(-2,8),axes=F,lwd=2,ylab="change in log-odds for svr",xlab="Weeks",main="Adherence to Peginterferon (Hypothesis Test)")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,8),axes=F,lwd=2,lty=6,ylab="change in log-odds for svr",xlab="Weeks")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,8),axes=F,lwd=2,lty=6,ylab="change in log-odds for svr",xlab="Weeks")
lines(rep(0,24))
Axis(side=1,at=c(1,6,12,18,24))
Axis(side=2)




plot(time, y1, type='l',ylim=c(-1,1),axes=F,ylab="change in log-odds for svr",xlab="Weeks",main="Effect of ISHAK on SVR in Peginterferon Analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-1,1),axes=F,col="blue",ylab="change in log-odds for svr",xlab="Weeks")
par(new=T)
plot(time, y3, type='l',ylim=c(-1,1),axes=F,col=2,ylab="change in log-odds for svr",xlab="Weeks")
lines(rep(0,24))
Axis(side=1,at=c(1,6,12,18,24))
Axis(side=2)


plot(time, y1, type='l',ylim=c(-2,2),axes=F,ylab="change in log-odds for males v females on svr",xlab="Weeks",main="Effect of Gender on SVR in Peginterferon Analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,2),axes=F,col="blue",ylab="change in log-odds for males v females on svr",xlab="Weeks")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,2),axes=F,col=2,ylab="change in log-odds for males v females on svr",xlab="Weeks")
lines(rep(0,24))
Axis(side=1,at=c(1,6,12,18,24))
Axis(side=2)



plot(time, y1, type='l',ylim=c(-2,3),axes=F,ylab="change in log-odds for caucasians v non-caucasians on svr",xlab="Weeks",main="Effect of Race on SVR in Peginterferon Analysis")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,3),axes=F,col="blue",ylab="change in log-odds for caucasians v non-caucasians on svr",xlab="Weeks")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,3),axes=F,col=2,ylab="change in log-odds for caucasians v non-caucasians on svr",xlab="Weeks")
lines(rep(0,24))
Axis(side=1,at=c(1,6,12,18,24))
Axis(side=2)







plot(time, y1, type='l',ylim=c(-2,8),ylab="change in log-odds for svr",xlab="Weeks",main="Compliance for Peginterferon (Hypothesis Tests)")
par(new=T)
plot(time, y2, type='l',ylim=c(-2,8),col="blue",ylab="change in log-odds for svr",xlab="Weeks")
par(new=T)
plot(time, y3, type='l',ylim=c(-2,8),col="red",ylab="change in log-odds for svr",xlab="Weeks")

lines(rep(0,24),lwd=2)








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

y1<-rep(0,24)

y4<-rep(0,24)
y5<-rep(0,24)


for (i in 1:24){
y1[i]<-beta[[i]][2]
y2[i]<-lcl[[i]][2]
y3[i]<-ucl[[i]][2] 
}

plot(time, y1, type='l',ylim=c(-2,8),ylab="change in log-odds for svr",xlab="Weeks",main="Compliance for Peginterferon (Confidence Bands)")
par(new=T)
plot(time, y4, type='l',ylim=c(-2,8),col="blue",ylab="change in log-odds for svr",xlab="Weeks")
par(new=T)
plot(time, y5, type='l',ylim=c(-2,8),col="red",ylab="change in log-odds for svr",xlab="Weeks")
lines(rep(0,24))




save.image("G:\\inside\\coursework\\tracs\\data\\WS_peg_f24")
















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
sigma<-(1/(401^2))*matrix(rowSums(rbind(apply(Z,2,function(x){x%*%t(x)}))),nrow(C[[1]]),nrow(C[[1]]))


t2<-t(int1)%*%solve(sigma)%*%int1
pval=1-pchisq(t2,nrow(as.matrix(C[[1]])))
 return(pval)
}


contrast<-rep(0,17)
contrast[17]<-1
C<-vector(mode="list",24)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",24)
c<-lapply(c,function(x){x<-0})

W<-vector(mode="list",24)
W<-lapply(W,function(x){x<-1})

T2(C,c,W)






contrast<-matrix(0,3,17)
contrast[1,8]<-1
contrast[2,9]<-1
contrast[3,10]<-1
C<-vector(mode="list",24)
C<-lapply(C,function(x){x<-contrast})
c<-vector(mode="list",24)
c<-lapply(c,function(x){x<-rep(0,3)})
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


contrast<-rep(0,17)
contrast[8]<-1
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
		Cov.t<-C.t%*%Cov[,,t,t]%*%t(C.t)
		Stat[t]<-t(T.t)%*%Cov.t%*%T.t
	}
	T3<-max(abs(Stat))
return(T3)
}









contrast<-rep(0,17)
contrast[2]<-1
C<-vector(mode="list",24)
C<-lapply(C,function(x){x<-t(contrast)})
c<-vector(mode="list",24)
c<-lapply(c,function(x){x<-0})



A.Ct.beta<-lapply(A.beta.curl2,function(x){apply(x,1,function(z){C[[t]]%*%z})})
A.Ct.beta.sm<-lapply(A.Ct.beta,function(x){apply(cbind(x),2,function(z){tri_smoothed(time,z,3)})})



Var.Ct.sm<-function(t,A){
X<-subset(peg_demo1,select = c(id,week,n,SEX,RACEW,MXAD,age,ISHAK,infect,education,insurance,employ,marital,alcohol))
time<-seq(1,24)
G.Ct<-matrix(0,length(C[[1]]%*%beta[[1]]),length(C[[1]]%*%beta[[1]]))

for (i in 1:401){
A.Ct.beta.sm.t<-t(A[[i]][t,])
G.i<-t(A.Ct.beta.sm.t)%*%A.Ct.beta.sm.t
G.Ct<-G.Ct+G.i
}
return(G.Ct)
}

EV.Ct<-vector("list",24)

for (t in 1:24){
EV.Ct[[t]]<-sqrt(diag(Var.Ct.sm(t,A.Ct.beta.sm)))
}


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





save.image("G:\\inside\\coursework\\tracs\\data\\anls6")



# Cross Validation



betas<-matrix(unlist(beta),24,17,byrow=T)


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



betas<-matrix(unlist(beta),24,17,byrow=T)
betas.sm2<-apply(betas,2,function(x){tri_smoothed(time,x,2)})
beta.sm2<-lapply(1:nrow(betas.sm2), function(i) betas.sm2[i,])

