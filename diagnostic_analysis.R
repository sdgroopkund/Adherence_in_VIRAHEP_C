sa<-read.table('sa.txt',header=T,sep='\t')
le<-read.table('le.txt',header=T,sep='\t')
vload_1<-read.table('vload_1.txt',header=T,sep='\t')


i=3
xvmt<-rib_peg[rib_peg$week==i,]
xmt1<-merge(xvmt,xmt,by="vhcid")
xmt1$svr[intersect(which(xmt1$n==0),which(xmt1$n2<1))]
mean(xmt1$svr[which(xmt1$n==0)])
xmt1$svr[which(xmt1$n2<1)]


i=2
xvmt<-rib_peg[rib_peg$week==i,]
xmt2<-merge(xvmt,xmt,by="vhcid")

mean(xmt2$svr[which(xmt2$n==0)])
xmt2$svr[which(xmt2$n2<1)]


i=4
xvmt<-rib_peg[rib_peg$week==i,]
xmt4<-merge(xvmt,xmt,by="vhcid")

mean(xmt4$svr[which(xmt4$n==0)])
xmt4$svr[which(xmt4$n2<1)]


22  84 120 145 161 172 195 207 213 222 231 254 261 271 275 284 320 322 339 363 389 391


xmt1$svr[(xmt1$n==0)]
which((xmt1$n==0))

id_w3_n0<-xmt1$vhcid[(xmt1$n==0)]
id_w3_n1<-xmt1$vhcid[(xmt1$n==1)]

id_w3_n0_n2<-xmt1$vhcid[intersect(which(xmt1$n==0),which(xmt1$n2<1))]


##########################################################################################
############ analysis on SA symptom scores ####################################
#################################################################################


######### overall #############

sa$SAOV[sa$SAOV=="B"]<-NA
sa$SAOV[sa$SAOV=="C"]<-NA
sa$SAOV[sa$SAOV=="D"]<-NA
sa$SAOV<-as.numeric(as.character(sa$SAOV))

sa_int<-sa[sa$vhcid%in%id_w3_n0,]
sa_rest<-sa[sa$vhcid%in%id_w3_n1,]

saov.int<-sa_int$SAOV[which(sa_int$TMPT=='TMT Day 28')]
saov.rest<-sa_rest$SAOV[which(sa_rest$TMPT=='TMT Day 28')]
saov.int<-sa_int$SAOV[which(sa_int$TMPT=='TMT Day 14')]
saov.rest<-sa_rest$SAOV[which(sa_rest$TMPT=='TMT Day 14')]
saov.int<-sa_int$SAOV[which(sa_int$TMPT=='TMT Day 7')]
saov.rest<-sa_rest$SAOV[which(sa_rest$TMPT=='TMT Day 7')]
saov.int<-sa_int$SAOV[which(sa_int$TMPT=='TMT Week 12')]
saov.rest<-sa_rest$SAOV[which(sa_rest$TMPT=='TMT Week 12')]
saov.int<-sa_int$SAOV[which(sa_int$TMPT=='TMT Week 16')]
saov.rest<-sa_rest$SAOV[which(sa_rest$TMPT=='TMT Week 16')]
saov.int<-sa_int$SAOV[which(sa_int$TMPT=='TMT Week 20')]
saov.rest<-sa_rest$SAOV[which(sa_rest$TMPT=='TMT Week 20')]
xsa<-NULL
ysa<-NULL

x<-saov.int
y<-saov.rest
xsa<-c(xsa,x)
ysa<-c(ysa,y)


xsa1<-xsa[-which(is.na(xsa)==T)]
ysa1<-ysa[-which(is.na(ysa)==T)]


plot(ecdf(ysa1),axes='F', verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='Symptom Score (Overall)',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(xsa),axes='F',verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='Symptom Score (Overall)',ylab='Cumulative Distribution Function')
legend(6,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)

plot(ecdf(ysa1),main='CF',xlab='symptom scores overall',ylab='Fn')
lines(ecdf(xsa))




mean(saov.int,na.rm=T)
mean(saov.rest,na.rm=T)
t.test(saov.int,saov.rest)

x<-c(rep(1,length(saov.int)),rep(0,length(saov.rest)))
y<-c(saov.int,saov.rest)
mod<-lm(y~x)
summary(mod)


########### muscle ache #############################



sa$SAMAC[sa$SAMAC=="B"]<-NA
sa$SAMAC[sa$SAMAC=="C"]<-NA
sa$SAMAC[sa$SAMAC=="D"]<-NA
sa$SAMAC<-as.numeric(as.character(sa$SAMAC))

sa_int<-sa[sa$vhcid%in%id_w3_n0,]
sa_rest<-sa[sa$vhcid%in%id_w3_n1,]

SAMAC.int<-sa_int$SAMAC[which(sa_int$TMPT=='TMT Day 28')]
SAMAC.rest<-sa_rest$SAMAC[which(sa_rest$TMPT=='TMT Day 28')]

SAMAC.int<-sa_int$SAMAC[which(sa_int$TMPT=='TMT Day 14')]
SAMAC.rest<-sa_rest$SAMAC[which(sa_rest$TMPT=='TMT Day 14')]

SAMAC.int<-sa_int$SAMAC[which(sa_int$TMPT=='TMT Day 7')]
SAMAC.rest<-sa_rest$SAMAC[which(sa_rest$TMPT=='TMT Day 7')]

SAMAC.int<-sa_int$SAMAC[which(sa_int$TMPT=='TMT Week 12')]
SAMAC.rest<-sa_rest$SAMAC[which(sa_rest$TMPT=='TMT Week 12')]

SAMAC.int<-sa_int$SAMAC[which(sa_int$TMPT=='TMT Week 16')]
SAMAC.rest<-sa_rest$SAMAC[which(sa_rest$TMPT=='TMT Week 16')]

SAMAC.int<-sa_int$SAMAC[which(sa_int$TMPT=='TMT Week 20')]
SAMAC.rest<-sa_rest$SAMAC[which(sa_rest$TMPT=='TMT Week 20')]

xmac<-NULL
ymac<-NULL

x<-SAMAC.int
y<-SAMAC.rest
xmac<-c(xmac,x)
ymac<-c(ymac,y)


xmac1<-xmac[-which(is.na(xmac)==T)]
ymac1<-ymac[-which(is.na(ymac)==T)]

plot(ecdf(ymac1),main='Pooled Cumulative Distribution Function',xlab='muscle ache',ylab='Cumulative Distribution Function')
lines(ecdf(xmac))
Axis(side=1)
Axis(side=2)

require(Hmisc)
Ecdf(ymac1,axes=F)
Ecdf(xmac,add=T,axes=F)

plot(ecdf(ymac1),axes='F', verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='Muscle Ache',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(xmac),axes='F',verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='Muscle Ache',ylab='Cumulative Distribution Function')
legend(6,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)

mean(SAMAC.int,na.rm=T)
mean(SAMAC.rest,na.rm=T)
t.test(SAMAC.int,SAMAC.rest)


###################### irritability #############################




sa$SAIRR[sa$SAIRR=="B"]<-NA
sa$SAIRR[sa$SAIRR=="C"]<-NA
sa$SAIRR[sa$SAIRR=="D"]<-NA
sa$SAIRR<-as.numeric(as.character(sa$SAIRR))

sa_int<-sa[sa$vhcid%in%id_w3_n0,]
sa_rest<-sa[sa$vhcid%in%id_w3_n1,]

SAIRR.int<-sa_int$SAIRR[which(sa_int$TMPT=='TMT Day 28')]
SAIRR.rest<-sa_rest$SAIRR[which(sa_rest$TMPT=='TMT Day 28')]

SAIRR.int<-sa_int$SAIRR[which(sa_int$TMPT=='TMT Day 14')]
SAIRR.rest<-sa_rest$SAIRR[which(sa_rest$TMPT=='TMT Day 14')]

SAIRR.int<-sa_int$SAIRR[which(sa_int$TMPT=='TMT Day 7')]
SAIRR.rest<-sa_rest$SAIRR[which(sa_rest$TMPT=='TMT Day 7')]

SAIRR.int<-sa_int$SAIRR[which(sa_int$TMPT=='TMT Week 12')]
SAIRR.rest<-sa_rest$SAIRR[which(sa_rest$TMPT=='TMT Week 12')]

SAIRR.int<-sa_int$SAIRR[which(sa_int$TMPT=='TMT Week 16')]
SAIRR.rest<-sa_rest$SAIRR[which(sa_rest$TMPT=='TMT Week 16')]

SAIRR.int<-sa_int$SAIRR[which(sa_int$TMPT=='TMT Week 20')]
SAIRR.rest<-sa_rest$SAIRR[which(sa_rest$TMPT=='TMT Week 20')]

xirr<-NULL
yirr<-NULL

x<-SAIRR.int
y<-SAIRR.rest
xirr<-c(xirr,x)
yirr<-c(yirr,y)


xirr1<-xirr[-which(is.na(xirr)==T)]
yirr1<-yirr[-which(is.na(yirr)==T)]

plot(ecdf(yirr1),axes='F', verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='Irritability',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(xirr),axes='F',verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='Irritability',ylab='Cumulative Distribution Function')
legend(6,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)




plot(ecdf(yirr1),main='CF',xlab='symptom scores irritability',ylab='Fn')
lines(ecdf(xirr))



mean(xsa,na.rm=T)
mean(ysa1,na.rm=T)
t.test(xsa,ysa1)




########################### headache ################################


sa$SAHAC[sa$SAHAC=="B"]<-NA
sa$SAHAC[sa$SAHAC=="C"]<-NA
sa$SAHAC[sa$SAHAC=="D"]<-NA
sa$SAHAC<-as.numeric(as.character(sa$SAHAC))

sa_int<-sa[sa$vhcid%in%id_w3_n0,]
sa_rest<-sa[sa$vhcid%in%id_w3_n1,]

SAHAC.int<-sa_int$SAHAC[which(sa_int$TMPT=='TMT Day 28')]
SAHAC.rest<-sa_rest$SAHAC[which(sa_rest$TMPT=='TMT Day 28')]

SAHAC.int<-sa_int$SAHAC[which(sa_int$TMPT=='TMT Day 14')]
SAHAC.rest<-sa_rest$SAHAC[which(sa_rest$TMPT=='TMT Day 14')]

SAHAC.int<-sa_int$SAHAC[which(sa_int$TMPT=='TMT Day 7')]
SAHAC.rest<-sa_rest$SAHAC[which(sa_rest$TMPT=='TMT Day 7')]

SAHAC.int<-sa_int$SAHAC[which(sa_int$TMPT=='TMT Week 12')]
SAHAC.rest<-sa_rest$SAHAC[which(sa_rest$TMPT=='TMT Week 12')]

SAHAC.int<-sa_int$SAHAC[which(sa_int$TMPT=='TMT Week 16')]
SAHAC.rest<-sa_rest$SAHAC[which(sa_rest$TMPT=='TMT Week 16')]

SAHAC.int<-sa_int$SAHAC[which(sa_int$TMPT=='TMT Week 20')]
SAHAC.rest<-sa_rest$SAHAC[which(sa_rest$TMPT=='TMT Week 20')]

xhac<-NULL
yhac<-NULL

x<-SAHAC.int
y<-SAHAC.rest
xhac<-c(xhac,x)
yhac<-c(yhac,y)


xhac1<-xhac[-which(is.na(xhac)==T)]
yhac1<-yhac[-which(is.na(yhac)==T)]


plot(ecdf(yhac1),axes='F', verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='Headache',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(xhac),axes='F',verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='Headache',ylab='Cumulative Distribution Function')
legend(6,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)


plot(ecdf(yhac1),main='CF',xlab='symptom scores headache',ylab='Fn')
lines(ecdf(xhac))



mean(SAHAC.int,na.rm=T)
mean(SAHAC.rest,na.rm=T)
t.test(xsa,ysa1)



########################### fatigue ########################################




sa$SAFAT[sa$SAFAT=="B"]<-NA
sa$SAFAT[sa$SAFAT=="C"]<-NA
sa$SAFAT[sa$SAFAT=="D"]<-NA
sa$SAFAT<-as.numeric(as.character(sa$SAFAT))

sa_int<-sa[sa$vhcid%in%id_w3_n0,]
sa_rest<-sa[sa$vhcid%in%id_w3_n1,]

SAFAT.int<-sa_int$SAFAT[which(sa_int$TMPT=='TMT Day 28')]
SAFAT.rest<-sa_rest$SAFAT[which(sa_rest$TMPT=='TMT Day 28')]

SAFAT.int<-sa_int$SAFAT[which(sa_int$TMPT=='TMT Day 14')]
SAFAT.rest<-sa_rest$SAFAT[which(sa_rest$TMPT=='TMT Day 14')]

SAFAT.int<-sa_int$SAFAT[which(sa_int$TMPT=='TMT Day 7')]
SAFAT.rest<-sa_rest$SAFAT[which(sa_rest$TMPT=='TMT Day 7')]

SAFAT.int<-sa_int$SAFAT[which(sa_int$TMPT=='TMT Week 12')]
SAFAT.rest<-sa_rest$SAFAT[which(sa_rest$TMPT=='TMT Week 12')]

SAFAT.int<-sa_int$SAFAT[which(sa_int$TMPT=='TMT Week 16')]
SAFAT.rest<-sa_rest$SAFAT[which(sa_rest$TMPT=='TMT Week 16')]

SAFAT.int<-sa_int$SAFAT[which(sa_int$TMPT=='TMT Week 20')]
SAFAT.rest<-sa_rest$SAFAT[which(sa_rest$TMPT=='TMT Week 20')]

xfat<-NULL
yfat<-NULL

x<-SAFAT.int
y<-SAFAT.rest
xfat<-c(xfat,x)
yfat<-c(yfat,y)


xfat1<-xfat[-which(is.na(xfat)==T)]
yfat1<-yfat[-which(is.na(yfat)==T)]

plot(ecdf(yfat1),axes='F', verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='Fatigue',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(xfat),axes='F',verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='Fatigue',ylab='Cumulative Distribution Function')
legend(6,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)


plot(ecdf(yfat1),main='CF',xlab='symptom scores fatigue',ylab='Fn')
lines(ecdf(xfat))




mean(SAFAT.int,na.rm=T)
mean(SAFAT.rest,na.rm=T)
t.test(xsa,ysa1)




########################### depression ################################################



sa$SADP[sa$SADP=="B"]<-NA
sa$SADP[sa$SADP=="C"]<-NA
sa$SADP[sa$SADP=="D"]<-NA
sa$SADP<-as.numeric(as.character(sa$SADP))

sa_int<-sa[sa$vhcid%in%id_w3_n0,]
sa_rest<-sa[sa$vhcid%in%id_w3_n1,]

SADP.int<-sa_int$SADP[which(sa_int$TMPT=='TMT Day 28')]
SADP.rest<-sa_rest$SADP[which(sa_rest$TMPT=='TMT Day 28')]

SADP.int<-sa_int$SADP[which(sa_int$TMPT=='TMT Day 14')]
SADP.rest<-sa_rest$SADP[which(sa_rest$TMPT=='TMT Day 14')]

SADP.int<-sa_int$SADP[which(sa_int$TMPT=='TMT Day 7')]
SADP.rest<-sa_rest$SADP[which(sa_rest$TMPT=='TMT Day 7')]

SADP.int<-sa_int$SADP[which(sa_int$TMPT=='TMT Week 12')]
SADP.rest<-sa_rest$SADP[which(sa_rest$TMPT=='TMT Week 12')]

SADP.int<-sa_int$SADP[which(sa_int$TMPT=='TMT Week 16')]
SADP.rest<-sa_rest$SADP[which(sa_rest$TMPT=='TMT Week 16')]

SADP.int<-sa_int$SADP[which(sa_int$TMPT=='TMT Week 20')]
SADP.rest<-sa_rest$SADP[which(sa_rest$TMPT=='TMT Week 20')]

xdep<-NULL
ydep<-NULL

x<-SADP.int
y<-SADP.rest
xdep<-c(xdep,x)
ydep<-c(ydep,y)


xdep1<-xdep[-which(is.na(xdep)==T)]
ydep1<-ydep[-which(is.na(ydep)==T)]


plot(ecdf(ydep1),axes='F', verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='Depression',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(xdep),axes='F',verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='Depression',ylab='Cumulative Distribution Function')
legend(6,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)

plot(ecdf(ydep1),main='CF',xlab='symptom scores depression',ylab='Fn')
lines(ecdf(xdep))


mean(SADP.int,na.rm=T)
mean(SADP.rest,na.rm=T)
t.test(xsa,ysa1)



##########################################################################################
############ analysis on LE platelet counts ####################################
#################################################################################




######################## Platelets #########################################

le$PC[le$PC=="A"]<-NA
le$PC[le$PC=="B"]<-NA
le$PC[le$PC=="C"]<-NA
le$PC[le$PC=="D"]<-NA
le$PC[le$PC=="E"]<-NA
le$PC[le$PC=="F"]<-NA

le$PC<-as.numeric(as.character(le$PC))

le_int<-le[le$vhcid%in%id_w3_n0,]
le_rest<-le[le$vhcid%in%id_w3_n1,]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Day 28')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Day 28')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Day 14')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Day 14')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Day 7')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Day 7')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 6')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 6')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 8')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 8')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 12')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 12')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 16')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 16')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 20')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 20')]


mean(pc.int,na.rm=T)
mean(pc.rest,na.rm=T)
t.test(pc.int,pc.rest)

x<-c(rep(1,length(pc.int)),rep(0,length(pc.rest)))
y<-c(pc.int,pc.rest)
mod<-lm(y~x)
summary(mod)


######################### WBC #####################################################


le$WBC[le$WBC=="A"]<-NA
le$WBC[le$WBC=="B"]<-NA
le$WBC[le$WBC=="C"]<-NA
le$WBC[le$WBC=="D"]<-NA
le$WBC[le$WBC=="E"]<-NA
le$WBC[le$WBC=="F"]<-NA

le$WBC<-as.numeric(as.character(le$WBC))

le_int<-le[le$vhcid%in%id_w3_n0,]
le_rest<-le[le$vhcid%in%id_w3_n1,]

WBC.int<-le_int$WBC[which(le_int$tmpt=='TMT Day 28')]
WBC.rest<-le_rest$WBC[which(le_rest$tmpt=='TMT Day 28')]

WBC.int<-le_int$WBC[which(le_int$tmpt=='TMT Day 14')]
WBC.rest<-le_rest$WBC[which(le_rest$tmpt=='TMT Day 14')]

WBC.int<-le_int$WBC[which(le_int$tmpt=='TMT Day 7')]
WBC.rest<-le_rest$WBC[which(le_rest$tmpt=='TMT Day 7')]

WBC.int<-le_int$WBC[which(le_int$tmpt=='TMT Week 6')]
WBC.rest<-le_rest$WBC[which(le_rest$tmpt=='TMT Week 6')]

WBC.int<-le_int$WBC[which(le_int$tmpt=='TMT Week 8')]
WBC.rest<-le_rest$WBC[which(le_rest$tmpt=='TMT Week 8')]

WBC.int<-le_int$WBC[which(le_int$tmpt=='TMT Week 12')]
WBC.rest<-le_rest$WBC[which(le_rest$tmpt=='TMT Week 12')]

WBC.int<-le_int$WBC[which(le_int$tmpt=='TMT Week 16')]
WBC.rest<-le_rest$WBC[which(le_rest$tmpt=='TMT Week 16')]

WBC.int<-le_int$WBC[which(le_int$tmpt=='TMT Week 20')]
WBC.rest<-le_rest$WBC[which(le_rest$tmpt=='TMT Week 20')]


xwbc<-NULL
ywbc<-NULL

x<-WBC.int
y<-WBC.rest
xwbc<-c(xwbc,x)
ywbc<-c(ywbc,y)


xwbc1<-xwbc[-which(is.na(xwbc)==T)]
ywbc1<-ywbc[-which(is.na(ywbc)==T)]


plot(ecdf(ywbc1),axes='F',xlim=c(0,15), verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='WBC Count',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(xwbc),axes='F',xlim=c(0,15),verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='WBC Count',ylab='Cumulative Distribution Function')
#lines(ecdf(xwbc))
legend(6,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)


plot(ecdf(ywbc1),main='CF',xlab='WBC count',ylab='Fn')
lines(ecdf(xwbc))




mean(WBC.int,na.rm=T)
mean(WBC.rest,na.rm=T)
t.test(WBC.int,WBC.rest)

x<-c(rep(1,length(WBC.int)),rep(0,length(WBC.rest)))
y<-c(WBC.int,WBC.rest)
mod<-lm(y~x)
summary(mod)




################# NPC ##################################################

le$NPC[le$NPC=="A"]<-NA
le$NPC[le$NPC=="B"]<-NA
le$NPC[le$NPC=="C"]<-NA
le$NPC[le$NPC=="D"]<-NA
le$NPC[le$NPC=="E"]<-NA
le$NPC[le$NPC=="F"]<-NA

le$NPC<-as.numeric(as.character(le$NPC))

le_int<-le[le$vhcid%in%id_w3_n0,]
le_rest<-le[le$vhcid%in%id_w3_n1,]

NPC.int<-le_int$NPC[which(le_int$tmpt=='TMT Day 28')]
NPC.rest<-le_rest$NPC[which(le_rest$tmpt=='TMT Day 28')]

NPC.int<-le_int$NPC[which(le_int$tmpt=='TMT Day 14')]
NPC.rest<-le_rest$NPC[which(le_rest$tmpt=='TMT Day 14')]

NPC.int<-le_int$NPC[which(le_int$tmpt=='TMT Day 7')]
NPC.rest<-le_rest$NPC[which(le_rest$tmpt=='TMT Day 7')]

NPC.int<-le_int$NPC[which(le_int$tmpt=='TMT Week 6')]
NPC.rest<-le_rest$NPC[which(le_rest$tmpt=='TMT Week 6')]

NPC.int<-le_int$NPC[which(le_int$tmpt=='TMT Week 8')]
NPC.rest<-le_rest$NPC[which(le_rest$tmpt=='TMT Week 8')]

NPC.int<-le_int$NPC[which(le_int$tmpt=='TMT Week 12')]
NPC.rest<-le_rest$NPC[which(le_rest$tmpt=='TMT Week 12')]

NPC.int<-le_int$NPC[which(le_int$tmpt=='TMT Week 16')]
NPC.rest<-le_rest$NPC[which(le_rest$tmpt=='TMT Week 16')]

NPC.int<-le_int$NPC[which(le_int$tmpt=='TMT Week 20')]
NPC.rest<-le_rest$NPC[which(le_rest$tmpt=='TMT Week 20')]
xnpc<-NULL
ynpc<-NULL

x<-NPC.int
y<-NPC.rest
xnpc<-c(xnpc,x)
ynpc<-c(ynpc,y)


xnpc1<-xnpc[-which(is.na(xnpc)==T)]
ynpc1<-ynpc[-which(is.na(ynpc)==T)]

plot(ecdf(ynpc1),axes='F',xlim=c(0,100), verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='NPC Count',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(xnpc),axes='F',xlim=c(0,100),verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='NPC Count',ylab='Cumulative Distribution Function')
#lines(ecdf(xwbc))
legend(60,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)

plot(ecdf(ynpc1),main='CF',xlab='NPC count',ylab='Fn')
lines(ecdf(xnpc))





mean(NPC.int,na.rm=T)
mean(NPC.rest,na.rm=T)
t.test(NPC.int,NPC.rest)

x<-c(rep(1,length(NPC.int)),rep(0,length(NPC.rest)))
y<-c(NPC.int,NPC.rest)
mod<-lm(y~x)
summary(mod)



##########################################################################################
############ analysis on VLOAD_1 vload_itt ####################################
#################################################################################




vload_int<-vload_1[vload_1$vhcid%in%id_w3_n0,]
vload_rest<-vload_1[vload_1$vhcid%in%id_w3_n1,]


vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Day 21')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Day 21')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Day 1')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Day 1')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Day 28')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Day 28')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Day 7')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Day 7')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Day 14')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Day 14')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Week 8')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Week 8')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Week 12')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Week 12')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Week 24')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Week 24')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Week 48')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Week 48')]


mean(vl.int,na.rm=T)
mean(vl.rest,na.rm=T)
t.test(vl.int,vl.rest)

x<-c(rep(1,length(vl.int)),rep(0,length(vl.rest)))
y<-c(vl.int,vl.rest)
mod<-lm(y~x)
summary(mod)
prop.test(c(sum(vl.int>599),sum(vl.rest>599)),c(length(vl.int),length(vl.rest)))



###########################################################################################
#################### subsequent adherance ##############################################
######################################################################################

vload_int<-vload_1[vload_1$vhcid%in%id_w3_n0,]
vload_rest<-vload_1[vload_1$vhcid%in%id_w3_n1,]

peg3.f24<-peg3[peg3$week<=24,]
peg.sub<-peg3.f24[peg3.f24$vhcid%in%id_w3_n0,]

pp<-vector("list",22)
for (i in 1:22){
	pp[[i]]<-c(id_w3_n0[i],peg.sub$n[peg.sub$vhcid==id_w3_n0[i]])
}

lapply(pp,function(x)mean(x[2:length(x)]))
lapply(pp,function(x)sum(x[2:6]))


##########################################################



vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Day 1')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Day 1')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Day 28')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Day 28')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Day 7')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Day 7')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Day 14')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Day 14')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Week 8')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Week 8')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Week 12')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Week 12')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Week 24')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Week 24')]

vl.int<-vload_int$vload_itt[which(vload_int$tmpt=='TMT Week 48')]
vl.rest<-vload_rest$vload_itt[which(vload_rest$tmpt=='TMT Week 48')]

vl.int.l<-log(vl.int)
vl.rest.l<-log(vl.rest)

x_vl<-NULL
y_vl<-NULL

x<-vl.int
y<-vl.rest
x_vl<-c(x_vl,x)
y_vl<-c(y_vl,y)



x_vl1<-x_vl[-which(is.na(x_vl)==T)]
y_vl1<-y_vl[-which(is.na(y_vl)==T)]

vl.int.l<-log(x_vl1)
vl.rest.l<-log(y_vl1)

plot(ecdf(vl.rest.l),axes='F',xlim=c(6,18), verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='Viral Load Scores in Ln scale',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(vl.int.l),axes='F',xlim=c(6,18),verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='Viral Load Scores in Ln scale',ylab='Cumulative Distribution Function')
#lines(ecdf(xwbc))
legend(14,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)


plot(ecdf(vl.int.l))
lines(ecdf(vl.rest.l))

x1<-c(x1,x)
y1<-c(y1,y)

x1<-x1[is.na(x1)==F]
y1<-y1[is.na(y1)==F]

y<-y+runif(length(y),0,0.001)
x<-x+runif(length(x),0,0.001)



library(CvM2SL2Test)
cvm <- cvmts.test(x, y)

## compute the p-value for the test.
 pval <- cvmts.pval(cvm, length(x), length(y) )


plot(ecdf(x1),main='CF',xlab='Viral Load scores in Ln scale',ylab='Fn')
lines(ecdf(y1))

plot(y,main='vload day 1',xlab='sample',ylab='value')
par(new=T)
plot(x,main='vload day 1',xlab='sample',ylab='value',)

lines(x)





z<-c(x,y)
ind<-c(rep(1,length(x)),rep(0,length(y)))
plot(z,pch= ifelse(ind > 0, "x", "o"),col= ifelse(ind > 0, "red", "black"))

Z<-cbind(z,ind)


par(new=T)

plot(ecdf(y))


#####################################################################

x2<-NULL
y2<-NULL

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Day 28')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Day 28')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Day 14')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Day 14')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Day 7')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Day 7')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 6')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 6')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 8')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 8')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 12')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 12')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 16')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 16')]

pc.int<-le_int$PC[which(le_int$tmpt=='TMT Week 20')]
pc.rest<-le_rest$PC[which(le_rest$tmpt=='TMT Week 20')]

x<-pc.int
y<-pc.rest
x2<-c(x2,x)
y2<-c(y2,y)

x1<-x[-which(is.na(x2)==T)]
y1<-y[-which(is.na(y2)==T)]

xplat<-x2
yplat<-y1

plot(ecdf(y1),axes='F',xlim=c(0,600), verticals=TRUE, pch=46,main='Pooled Cumulative Distribution Function',xlab='Platelet Count',ylab='Cumulative Distribution Function')
par(new=T)
plot(ecdf(x2),axes='F',xlim=c(0,600),verticals=TRUE, pch=46,lty=2,col='red',main='Pooled Cumulative Distribution Function',xlab='Platelet Count',ylab='Cumulative Distribution Function')
#lines(ecdf(xwbc))
legend(300,0.2 ,lty=c(2,1),col=c(2,1), c('Group 1','Group 2'), text.col = c(2,1))
Axis(side=1)
Axis(side=2)

plot(ecdf(y1),main='CF',xlab='Platelet Counts: Week 20',ylab='Fn')
lines(ecdf(x))


x21<-x2[-which(is.na(x2)==T)]
y21<-y2[-which(is.na(y2)==T)]

plot(ecdf(y21),main='CF',xlab='Platelet Counts',ylab='Fn')
lines(ecdf(x21))

y<-y[is.na(y)==F]
x<-x[is.na(x)==F]

plot(ecdf(y),main='CF week 6',xlab='platelet counts',ylab='Fn')
lines(ecdf(x))





y<-y+runif(length(y),0,0.001)
x<-x+runif(length(x),0,0.001)


library(CvM2SL2Test)


cvm <- cvmts.test(x, y)

## compute the p-value for the test.
 pval <- cvmts.pval(cvm, length(x), length(y) )



plot(ecdf(x))
plot(ecdf(y))
par(new=T)
plot(ecdf(x))
plot(ecdf(x))
par(new=T)
plot(ecdf(y))



vl.int_d1<-vload_int[which(vload_int$tmpt=='TMT Day 1'),c(2,13)]
vl.int_d1<-cbind(vl.int_d1,rep(1,nrow(vl.int_d1)))
colnames(vl.int_d1)<-c("vhcid","vload_d1","ind_d1")
vl.rest_d1<-vload_rest[which(vload_rest$tmpt=='TMT Day 1'),c(2,13)]
vl.rest_d1<-cbind(vl.rest_d1,rep(1,nrow(vl.rest_d1)))
colnames(vl.rest_d1)<-c("vhcid","vload_d1","ind_d1")


vl.int_d7<-vload_int[which(vload_int$tmpt=='TMT Day 7'),c(2,13)]
vl.int_d7<-cbind(vl.int_d7,rep(1,nrow(vl.int_d7)))
colnames(vl.int_d7)<-c("vhcid","vload_d7","ind_d7")
vl.rest_d7<-vload_rest[which(vload_rest$tmpt=='TMT Day 7'),c(2,13)]
vl.rest_d7<-cbind(vl.rest_d7,rep(1,nrow(vl.rest_d7)))
colnames(vl.rest_d7)<-c("vhcid","vload_d7","ind_d7")


vl.int_d14<-vload_int[which(vload_int$tmpt=='TMT Day 14'),c(2,13)]
vl.int_d14<-cbind(vl.int_d14,rep(1,nrow(vl.int_d14)))
colnames(vl.int_d14)<-c("vhcid","vload_d14","ind_d14")
vl.rest_d14<-vload_rest[which(vload_rest$tmpt=='TMT Day 14'),c(2,13)]
vl.rest_d14<-cbind(vl.rest_d14,rep(1,nrow(vl.rest_d14)))
colnames(vl.rest_d14)<-c("vhcid","vload_d14","ind_d14")

vl.int_d28<-vload_int[which(vload_int$tmpt=='TMT Day 28'),c(2,13)]
vl.int_d28<-cbind(vl.int_d28,rep(1,nrow(vl.int_d28)))
colnames(vl.int_d28)<-c("vhcid","vload_d28","ind_d28")
vl.rest_d28<-vload_rest[which(vload_rest$tmpt=='TMT Day 28'),c(2,13)]
vl.rest_d28<-cbind(vl.rest_d28,rep(1,nrow(vl.rest_d28)))
colnames(vl.rest_d28)<-c("vhcid","vload_d28","ind_d28")


vl.int_w8<-vload_int[which(vload_int$tmpt=='TMT Week 8'),c(2,13)]
vl.int_w8<-cbind(vl.int_w8,rep(1,nrow(vl.int_w8)))
colnames(vl.int_w8)<-c("vhcid","vload_w8","ind_w8")
vl.rest_w8<-vload_rest[which(vload_rest$tmpt=='TMT Week 8'),c(2,13)]
vl.rest_w8<-cbind(vl.rest_w8,rep(1,nrow(vl.rest_w8)))
colnames(vl.rest_w8)<-c("vhcid","vload_w8","ind_w8")

vl.int_w12<-vload_int[which(vload_int$tmpt=='TMT Week 12'),c(2,13)]
vl.int_w12<-cbind(vl.int_w12,rep(1,nrow(vl.int_w12)))
colnames(vl.int_w12)<-c("vhcid","vload_w12","ind_w12")
vl.rest_w12<-vload_rest[which(vload_rest$tmpt=='TMT Week 12'),c(2,13)]
vl.rest_w12<-cbind(vl.rest_w12,rep(1,nrow(vl.rest_w12)))
colnames(vl.rest_w12)<-c("vhcid","vload_w12","ind_w12")

vl.int_w24<-vload_int[which(vload_int$tmpt=='TMT Week 24'),c(2,13)]
vl.int_w24<-cbind(vl.int_w24,rep(1,nrow(vl.int_w24)))
colnames(vl.int_w24)<-c("vhcid","vload_w24","ind_w24")
vl.rest_w24<-vload_rest[which(vload_rest$tmpt=='TMT Week 24'),c(2,13)]
vl.rest_w24<-cbind(vl.rest_w24,rep(1,nrow(vl.rest_w24)))
colnames(vl.rest_w24)<-c("vhcid","vload_w24","ind_w24")

vl.int_w48<-vload_int[which(vload_int$tmpt=='TMT Week 48'),c(2,13)]
vl.int_w48<-cbind(vl.int_w48,rep(1,nrow(vl.int_w48)))
colnames(vl.int_w48)<-c("vhcid","vload_w48","ind_w48")
vl.rest_w48<-vload_rest[which(vload_rest$tmpt=='TMT Week 48'),c(2,13)]
vl.rest_w48<-cbind(vl.rest_w48,rep(1,nrow(vl.rest_w48)))
colnames(vl.rest_w48)<-c("vhcid","vload_w48","ind_w48")

data<-merge(vl.int_d1,vl.int_d7,all=T,by="vhcid")
data<-merge(vl.int_d14,data,all=T,by="vhcid")
data<-merge(vl.int_d28,data,all=T,by="vhcid")
data<-merge(vl.int_w8,data,all=T,by="vhcid")
data<-merge(vl.int_w12,data,all=T,by="vhcid")
data<-merge(vl.int_w24,data,all=T,by="vhcid")

data1<-merge(vl.rest_d1,vl.rest_d7,all=T,by="vhcid")
data1<-merge(vl.rest_d14,data1,all=T,by="vhcid")
data1<-merge(vl.rest_d28,data1,all=T,by="vhcid")
data1<-merge(vl.rest_w8,data1,all=T,by="vhcid")
data1<-merge(vl.rest_w12,data1,all=T,by="vhcid")
data1<-merge(vl.rest_w24,data1,all=T,by="vhcid")



vload<-data[,c(2,4,6,8,10,12,14)]
vload<-as.vector(as.matrix(vload))
vl.int<-vload[-which(is.na(vload)==T)]
vload2<-data1[,c(2,4,6,8,10,12,14)]
vload2<-as.vector(as.matrix(vload2))
vl.rest<-vload2[-which(is.na(vload2)==T)]

v1<-log(vl.rest)
v2<-log(vl.int)

kk<-ecdf(v1)
rr<-ecdf(v2)
kr<-sort(union(knots(kk),knots(rr)))
prob.G<-kk(kr)
prob.F<-rr(kr)
jumps.F<-c(prob.F[1],prob.F[-1]-prob.F[-length(prob.F)])
jumps.G<-c(prob.G[1],prob.G[-1]-prob.G[-length(prob.G)])



bs.stat<-rep(0,5000)
bs.stat2<-rep(0,5000)
for (i in 1:5000){
bs1<-data[sample(nrow(data),nrow(data),replace=T),]
bs2<-data1[sample(nrow(data1),nrow(data1),replace=T),]
vl.bs1<-bs1[,c(2,4,6,8,10,12,14)]
vl.bs1<-as.vector(as.matrix(vl.bs1))
vl.int.bs1<-vl.bs1[-which(is.na(vl.bs1)==T)]
vl.bs2<-bs2[,c(2,4,6,8,10,12,14)]
vl.bs2<-as.vector(as.matrix(vl.bs2))
vl.rest.bs2<-vl.bs2[-which(is.na(vl.bs2)==T)]
v1.bs<-log(vl.rest.bs2)
v2.bs<-log(vl.int.bs1)
kk.bs<-ecdf(v1.bs)
rr.bs<-ecdf(v2.bs)
prob.Gh<-kk.bs(kr)
prob.Fh<-rr.bs(kr)
jumps.Fh<-c(prob.Fh[1],prob.Fh[-1]-prob.Fh[-length(prob.Fh)])
jumps.Gh<-c(prob.Gh[1],prob.Gh[-1]-prob.Gh[-length(prob.Gh)])
FG<-cbind(kr,prob.F,jumps.F,prob.Fh,jumps.Fh,prob.G,jumps.G,prob.Gh,jumps.Gh)
bs.stat[i]<-sum((FG[,4]+FG[,6]-FG[,2]-FG[,8])^2*((length(v2)*FG[,3]+length(v1)*FG[,7])/(length(v1)+length(v2))))
bs.stat2[i]<-sum((FG[,4]+FG[,6]-FG[,2]-FG[,8])^2*((FG[,3]+FG[,7])/2))
}
test.stat<-sum((FG[,2]-FG[,6])^2*((length(v2)*FG[,3]+length(v1)*FG[,7])/(length(v1)+length(v2))))
test.stat2<-sum((FG[,2]-FG[,6])^2*((FG[,3]+FG[,7])/2))

pval<-mean(bs.stat > test.stat)

plot(ecdf(v2),main='CF day 1',xlab='sample',ylab='Fn')
lines(ecdf(v1))
















################### Corrected Version #########################


vl.int_d1<-vload_int[which(vload_int$tmpt=='TMT Day 1'),c(2,13)]
vl.int_d1<-cbind(vl.int_d1,rep(1,nrow(vl.int_d1)))
colnames(vl.int_d1)<-c("vhcid","vload_d1","ind_d1")
vl.rest_d1<-vload_rest[which(vload_rest$tmpt=='TMT Day 1'),c(2,13)]
vl.rest_d1<-cbind(vl.rest_d1,rep(1,nrow(vl.rest_d1)))
colnames(vl.rest_d1)<-c("vhcid","vload_d1","ind_d1")


vl.int_d7<-vload_int[which(vload_int$tmpt=='TMT Day 7'),c(2,13)]
vl.int_d7<-cbind(vl.int_d7,rep(1,nrow(vl.int_d7)))
colnames(vl.int_d7)<-c("vhcid","vload_d7","ind_d7")
vl.rest_d7<-vload_rest[which(vload_rest$tmpt=='TMT Day 7'),c(2,13)]
vl.rest_d7<-cbind(vl.rest_d7,rep(1,nrow(vl.rest_d7)))
colnames(vl.rest_d7)<-c("vhcid","vload_d7","ind_d7")


vl.int_d14<-vload_int[which(vload_int$tmpt=='TMT Day 14'),c(2,13)]
vl.int_d14<-cbind(vl.int_d14,rep(1,nrow(vl.int_d14)))
colnames(vl.int_d14)<-c("vhcid","vload_d14","ind_d14")
vl.rest_d14<-vload_rest[which(vload_rest$tmpt=='TMT Day 14'),c(2,13)]
vl.rest_d14<-cbind(vl.rest_d14,rep(1,nrow(vl.rest_d14)))
colnames(vl.rest_d14)<-c("vhcid","vload_d14","ind_d14")

vl.int_d28<-vload_int[which(vload_int$tmpt=='TMT Day 28'),c(2,13)]
vl.int_d28<-cbind(vl.int_d28,rep(1,nrow(vl.int_d28)))
colnames(vl.int_d28)<-c("vhcid","vload_d28","ind_d28")
vl.rest_d28<-vload_rest[which(vload_rest$tmpt=='TMT Day 28'),c(2,13)]
vl.rest_d28<-cbind(vl.rest_d28,rep(1,nrow(vl.rest_d28)))
colnames(vl.rest_d28)<-c("vhcid","vload_d28","ind_d28")


vl.int_w8<-vload_int[which(vload_int$tmpt=='TMT Week 8'),c(2,13)]
vl.int_w8<-cbind(vl.int_w8,rep(1,nrow(vl.int_w8)))
colnames(vl.int_w8)<-c("vhcid","vload_w8","ind_w8")
vl.rest_w8<-vload_rest[which(vload_rest$tmpt=='TMT Week 8'),c(2,13)]
vl.rest_w8<-cbind(vl.rest_w8,rep(1,nrow(vl.rest_w8)))
colnames(vl.rest_w8)<-c("vhcid","vload_w8","ind_w8")

vl.int_w12<-vload_int[which(vload_int$tmpt=='TMT Week 12'),c(2,13)]
vl.int_w12<-cbind(vl.int_w12,rep(1,nrow(vl.int_w12)))
colnames(vl.int_w12)<-c("vhcid","vload_w12","ind_w12")
vl.rest_w12<-vload_rest[which(vload_rest$tmpt=='TMT Week 12'),c(2,13)]
vl.rest_w12<-cbind(vl.rest_w12,rep(1,nrow(vl.rest_w12)))
colnames(vl.rest_w12)<-c("vhcid","vload_w12","ind_w12")

vl.int_w24<-vload_int[which(vload_int$tmpt=='TMT Week 24'),c(2,13)]
vl.int_w24<-cbind(vl.int_w24,rep(1,nrow(vl.int_w24)))
colnames(vl.int_w24)<-c("vhcid","vload_w24","ind_w24")
vl.rest_w24<-vload_rest[which(vload_rest$tmpt=='TMT Week 24'),c(2,13)]
vl.rest_w24<-cbind(vl.rest_w24,rep(1,nrow(vl.rest_w24)))
colnames(vl.rest_w24)<-c("vhcid","vload_w24","ind_w24")
vl.rest_w24$ind_w24[is.na(vl.rest_w24$vload_w24)==T]<-NA
vl.int_w24$ind_w24[is.na(vl.int_w24$vload_w24)==T]<-NA

vl.int_w48<-vload_int[which(vload_int$tmpt=='TMT Week 48'),c(2,13)]
vl.int_w48<-cbind(vl.int_w48,rep(1,nrow(vl.int_w48)))
colnames(vl.int_w48)<-c("vhcid","vload_w48","ind_w48")
vl.rest_w48<-vload_rest[which(vload_rest$tmpt=='TMT Week 48'),c(2,13)]
vl.rest_w48<-cbind(vl.rest_w48,rep(1,nrow(vl.rest_w48)))
colnames(vl.rest_w48)<-c("vhcid","vload_w48","ind_w48")

data<-merge(vl.int_d1,vl.int_d7,all=T,by="vhcid")
data<-merge(vl.int_d14,data,all=T,by="vhcid")
data<-merge(vl.int_d28,data,all=T,by="vhcid")
data<-merge(vl.int_w8,data,all=T,by="vhcid")
data<-merge(vl.int_w12,data,all=T,by="vhcid")
data<-merge(vl.int_w24,data,all=T,by="vhcid")

data1<-merge(vl.rest_d1,vl.rest_d7,all=T,by="vhcid")
data1<-merge(vl.rest_d14,data1,all=T,by="vhcid")
data1<-merge(vl.rest_d28,data1,all=T,by="vhcid")
data1<-merge(vl.rest_w8,data1,all=T,by="vhcid")
data1<-merge(vl.rest_w12,data1,all=T,by="vhcid")
data1<-merge(vl.rest_w24,data1,all=T,by="vhcid")

vload<-data[,c(2,4,6,8,10,12,14)]
vload<-as.vector(as.matrix(vload))
vl.int<-vload[-which(is.na(vload)==T)]
vload2<-data1[,c(2,4,6,8,10,12,14)]
vload2<-as.vector(as.matrix(vload2))
vl.rest<-vload2[-which(is.na(vload2)==T)]


v1.int_d1<-log(vl.int_d1[,2])
v1.int_d7<-log(vl.int_d7[,2])
v1.int_d14<-log(vl.int_d14[,2])
v1.int_d28<-log(vl.int_d28[,2])
v1.int_w8<-log(vl.int_w8[,2])
v1.int_w12<-log(vl.int_w12[,2])
v1.int_w24<-log(as.vector(na.omit(vl.int_w24[,2])))

v1.rest_d1<-log(vl.rest_d1[,2])
v1.rest_d7<-log(vl.rest_d7[,2])
v1.rest_d14<-log(vl.rest_d14[,2])
v1.rest_d28<-log(vl.rest_d28[,2])
v1.rest_w8<-log(vl.rest_w8[,2])
v1.rest_w12<-log(vl.rest_w12[,2])
v1.rest_w24<-log(as.vector(na.omit(vl.rest_w24[,2])))


rr_d1<-ecdf(v1.int_d1)
rr_d7<-ecdf(v1.int_d7)
rr_d14<-ecdf(v1.int_d14)
rr_d28<-ecdf(v1.int_d28)
rr_w8<-ecdf(v1.int_w8)
rr_w12<-ecdf(v1.int_w12)
rr_w24<-ecdf(v1.int_w24)


kk_d1<-ecdf(v1.rest_d1)
kk_d7<-ecdf(v1.rest_d7)
kk_d14<-ecdf(v1.rest_d14)
kk_d28<-ecdf(v1.rest_d28)
kk_w8<-ecdf(v1.rest_w8)
kk_w12<-ecdf(v1.rest_w12)
kk_w24<-ecdf(v1.rest_w24)

kr.new<-sort(union(knots(kk_d1),union(knots(kk_d7),union(knots(kk_d14),union(knots(kk_d28),union(knots(kk_w8),union(knots(kk_w12),union(knots(kk_w24),union(knots(rr_d1),union(knots(rr_d7),union(knots(rr_d14),union(knots(rr_d28),union(knots(rr_w8),union(knots(rr_w12),knots(rr_w24)))))))))))))))

prob.G_d1<-kk_d1(kr.new)
prob.G_d7<-kk_d7(kr.new)
prob.G_d14<-kk_d14(kr.new)
prob.G_d28<-kk_d28(kr.new)
prob.G_w8<-kk_w8(kr.new)
prob.G_w12<-kk_w12(kr.new)
prob.G_w24<-kk_w24(kr.new)


prob.F_d1<-rr_d1(kr.new)
prob.F_d7<-rr_d7(kr.new)
prob.F_d14<-rr_d14(kr.new)
prob.F_d28<-rr_d28(kr.new)
prob.F_w8<-rr_w8(kr.new)
prob.F_w12<-rr_w12(kr.new)
prob.F_w24<-rr_w24(kr.new)


jumps.F_d1<-c(prob.F_d1[1],prob.F_d1[-1]-prob.F_d1[-length(prob.F_d1)])
jumps.F_d7<-c(prob.F_d7[1],prob.F_d7[-1]-prob.F_d7[-length(prob.F_d7)])
jumps.F_d14<-c(prob.F_d14[1],prob.F_d14[-1]-prob.F_d14[-length(prob.F_d14)])
jumps.F_d28<-c(prob.F_d28[1],prob.F_d28[-1]-prob.F_d28[-length(prob.F_d28)])
jumps.F_w8<-c(prob.F_w8[1],prob.F_w8[-1]-prob.F_w8[-length(prob.F_w8)])
jumps.F_w12<-c(prob.F_w12[1],prob.F_w12[-1]-prob.F_w12[-length(prob.F_w12)])
jumps.F_w24<-c(prob.F_w24[1],prob.F_w24[-1]-prob.F_w24[-length(prob.F_w24)])


jumps.G_d1<-c(prob.G_d1[1],prob.G_d1[-1]-prob.G_d1[-length(prob.G_d1)])
jumps.G_d7<-c(prob.G_d7[1],prob.G_d7[-1]-prob.G_d7[-length(prob.G_d7)])
jumps.G_d14<-c(prob.G_d14[1],prob.G_d14[-1]-prob.G_d14[-length(prob.G_d14)])
jumps.G_d28<-c(prob.G_d28[1],prob.G_d28[-1]-prob.G_d28[-length(prob.G_d28)])
jumps.G_w8<-c(prob.G_w8[1],prob.G_w8[-1]-prob.G_w8[-length(prob.G_w8)])
jumps.G_w12<-c(prob.G_w12[1],prob.G_w12[-1]-prob.G_w12[-length(prob.G_w12)])
jumps.G_w24<-c(prob.G_w24[1],prob.G_w24[-1]-prob.G_w24[-length(prob.G_w24)])

S<-Sys.time()




bs.stat<-rep(0,500)
bs.stat2<-rep(0,500)

for (i in 1:500){

eps1<-rexp(nrow(data),1)
eps2<-rexp(nrow(data1),1)

data_bs<-data
data_bs[,c(3,5,7,9,11,13,15)]<-data_bs[,c(3,5,7,9,11,13,15)]*eps1
data_bs[,c(2,4,6,8,10,12,14)]<-log(data_bs[,c(2,4,6,8,10,12,14)])

data1_bs<-data1
data1_bs[,c(3,5,7,9,11,13,15)]<-data1_bs[,c(3,5,7,9,11,13,15)]*eps2
data1_bs[,c(2,4,6,8,10,12,14)]<-log(data1_bs[,c(2,4,6,8,10,12,14)])


F_bs<-matrix(0,length(kr.new),7)
for (k in c(2,4,6,8,10,12,14)){ 
F_bs[,k/2]<-sapply(kr.new,function(x){sum((as.vector(na.omit(data_bs[,k]))<=x)*as.vector(na.omit(data_bs[,k+1])))/sum(as.vector(na.omit(data_bs[,k+1])))})
}
Prob_Fh<-rowSums(F_bs)

G_bs<-matrix(0,length(kr.new),7)
for (k in c(2,4,6,8,10,12,14)){ 
G_bs[,k/2]<-sapply(kr.new,function(x){sum((as.vector(na.omit(data1_bs[,k]))<=x)*as.vector(na.omit(data1_bs[,k+1])))/sum(as.vector(na.omit(data1_bs[,k+1])))})
}
Prob_Gh<-rowSums(G_bs)


FG<-cbind(kr.new,prob.F_d1,prob.F_d7,prob.F_d14,prob.F_d28,prob.F_w8,prob.F_w12,prob.F_w24,
	jumps.F_d1,jumps.F_d7,jumps.F_d14,jumps.F_d28,jumps.F_w8,jumps.F_w12,jumps.F_w24,
	prob.G_d1,prob.G_d7,prob.G_d14,prob.G_d28,prob.G_w8,prob.G_w12,prob.G_w24,
	jumps.G_d1,jumps.G_d7,jumps.G_d14,jumps.G_d28,jumps.G_w8,jumps.G_w12,jumps.G_w24)



bs_sum_jump_size<-Prob_Fh-rowSums(FG[,c(2:8)])-Prob_Gh+rowSums(FG[,c(16:22)])
sum_jump_size<-rowSums(FG[,c(2:8)])-rowSums(FG[,c(16:22)])
jump_size<-(length(v2)*rowSums(FG[,c(9:15)])+length(v1)*rowSums(FG[,c(23:29)]))/(length(v2)+length(v1))

bs.stat[i]<-sum(bs_sum_jump_size^2*jump_size)
}

test.stat<-sum(sum_jump_size^2*jump_size)
pval<-mean(bs.stat > test.stat)

E<-Sys.time()-S



###########################################################################################################








########################### platelets ####################################



pc.int_w1<-le_int[which(le_int$tmpt=='TMT Day 7'),c(2,13)]
pc.int_w1<-cbind(pc.int_w1,rep(1,nrow(pc.int_w1)))
colnames(pc.int_w1)<-c("vhcid","pc_w1","ind_w1")
pc.rest_w1<-le_rest[which(le_rest$tmpt=='TMT Day 7'),c(2,13)]
pc.rest_w1<-cbind(pc.rest_w1,rep(1,nrow(pc.rest_w1)))
colnames(pc.rest_w1)<-c("vhcid","pc_w1","ind_w1")
pc.rest_w1$ind_w1[is.na(pc.rest_w1$pc_w1)==T]<-NA
pc.int_w1$ind_w1[is.na(pc.int_w1$pc_w1)==T]<-NA


pc.int_w2<-le_int[which(le_int$tmpt=='TMT Day 14'),c(2,13)]
pc.int_w2<-cbind(pc.int_w2,rep(1,nrow(pc.int_w2)))
colnames(pc.int_w2)<-c("vhcid","pc_w2","ind_w2")
pc.rest_w2<-le_rest[which(le_rest$tmpt=='TMT Day 14'),c(2,13)]
pc.rest_w2<-cbind(pc.rest_w2,rep(1,nrow(pc.rest_w2)))
colnames(pc.rest_w2)<-c("vhcid","pc_w2","ind_w2")
pc.rest_w2$ind_w2[is.na(pc.rest_w2$pc_w2)==T]<-NA
pc.int_w2$ind_w2[is.na(pc.int_w2$pc_w2)==T]<-NA


pc.int_w4<-le_int[which(le_int$tmpt=='TMT Day 28'),c(2,13)]
pc.int_w4<-cbind(pc.int_w4,rep(1,nrow(pc.int_w4)))
colnames(pc.int_w4)<-c("vhcid","pc_w4","ind_w4")
pc.rest_w4<-le_rest[which(le_rest$tmpt=='TMT Day 28'),c(2,13)]
pc.rest_w4<-cbind(pc.rest_w4,rep(1,nrow(pc.rest_w4)))
colnames(pc.rest_w4)<-c("vhcid","pc_w4","ind_w4")
pc.rest_w4$ind_w4[is.na(pc.rest_w4$pc_w4)==T]<-NA
pc.int_w4$ind_w4[is.na(pc.int_w4$pc_w4)==T]<-NA


pc.int_w6<-le_int[which(le_int$tmpt=='TMT Week 6'),c(2,13)]
pc.int_w6<-cbind(pc.int_w6,rep(1,nrow(pc.int_w6)))
colnames(pc.int_w6)<-c("vhcid","pc_w6","ind_w6")
pc.rest_w6<-le_rest[which(le_rest$tmpt=='TMT Week 6'),c(2,13)]
pc.rest_w6<-cbind(pc.rest_w6,rep(1,nrow(pc.rest_w6)))
colnames(pc.rest_w6)<-c("vhcid","pc_w6","ind_w6")
pc.rest_w6$ind_w6[is.na(pc.rest_w6$pc_w6)==T]<-NA
pc.int_w6$ind_w6[is.na(pc.int_w6$pc_w6)==T]<-NA


pc.int_w8<-le_int[which(le_int$tmpt=='TMT Week 8'),c(2,13)]
pc.int_w8<-cbind(pc.int_w8,rep(1,nrow(pc.int_w8)))
colnames(pc.int_w8)<-c("vhcid","pc_w8","ind_w8")
pc.rest_w8<-le_rest[which(le_rest$tmpt=='TMT Week 8'),c(2,13)]
pc.rest_w8<-cbind(pc.rest_w8,rep(1,nrow(pc.rest_w8)))
colnames(pc.rest_w8)<-c("vhcid","pc_w8","ind_w8")
pc.rest_w8$ind_w8[is.na(pc.rest_w8$pc_w8)==T]<-NA
pc.int_w8$ind_w8[is.na(pc.int_w8$pc_w8)==T]<-NA


pc.int_w12<-le_int[which(le_int$tmpt=='TMT Week 12'),c(2,13)]
pc.int_w12<-cbind(pc.int_w12,rep(1,nrow(pc.int_w12)))
colnames(pc.int_w12)<-c("vhcid","pc_w12","ind_w12")
pc.rest_w12<-le_rest[which(le_rest$tmpt=='TMT Week 12'),c(2,13)]
pc.rest_w12<-cbind(pc.rest_w12,rep(1,nrow(pc.rest_w12)))
colnames(pc.rest_w12)<-c("vhcid","pc_w12","ind_w12")
pc.rest_w12$ind_w12[is.na(pc.rest_w12$pc_w12)==T]<-NA
pc.int_w12$ind_w12[is.na(pc.int_w12$pc_w12)==T]<-NA


pc.int_w16<-le_int[which(le_int$tmpt=='TMT Week 16'),c(2,13)]
pc.int_w16<-cbind(pc.int_w16,rep(1,nrow(pc.int_w16)))
colnames(pc.int_w16)<-c("vhcid","pc_w16","ind_w16")
pc.rest_w16<-le_rest[which(le_rest$tmpt=='TMT Week 16'),c(2,13)]
pc.rest_w16<-cbind(pc.rest_w16,rep(1,nrow(pc.rest_w16)))
colnames(pc.rest_w16)<-c("vhcid","pc_w16","ind_w16")
pc.rest_w16$ind_w16[is.na(pc.rest_w16$pc_w16)==T]<-NA
pc.int_w16$ind_w16[is.na(pc.int_w16$pc_w16)==T]<-NA


pc.int_w20<-le_int[which(le_int$tmpt=='TMT Week 20'),c(2,13)]
pc.int_w20<-cbind(pc.int_w20,rep(1,nrow(pc.int_w20)))
colnames(pc.int_w20)<-c("vhcid","pc_w20","ind_w20")
pc.rest_w20<-le_rest[which(le_rest$tmpt=='TMT Week 20'),c(2,13)]
pc.rest_w20<-cbind(pc.rest_w20,rep(1,nrow(pc.rest_w20)))
colnames(pc.rest_w20)<-c("vhcid","pc_w20","ind_w20")
pc.rest_w20$ind_w20[is.na(pc.rest_w20$pc_w20)==T]<-NA
pc.int_w20$ind_w20[is.na(pc.int_w20$pc_w20)==T]<-NA


pc.int_d1<-le_int[which(le_int$tmpt=='Baseline'),c(2,13)]
pc.int_d1<-cbind(pc.int_d1,rep(1,nrow(pc.int_d1)))
colnames(pc.int_d1)<-c("vhcid","pc_d1","ind_d1")
pc.rest_d1<-le_rest[which(le_rest$tmpt=='Baseline'),c(2,13)]
pc.rest_d1<-cbind(pc.rest_d1,rep(1,nrow(pc.rest_d1)))
colnames(pc.rest_d1)<-c("vhcid","pc_d1","ind_d1")
pc.rest_d1$ind_d1[is.na(pc.rest_d1$pc_d1)==T]<-NA
pc.int_d1$ind_d1[is.na(pc.int_d1$pc_d1)==T]<-NA



pc.int_w24<-le_int[which(le_int$tmpt=='TMT Week 24'),c(2,13)]
pc.int_w24<-cbind(pc.int_w24,rep(1,nrow(pc.int_w24)))
colnames(pc.int_w24)<-c("vhcid","pc_w24","ind_w24")
pc.rest_w24<-le_rest[which(le_rest$tmpt=='TMT Week 24'),c(2,13)]
pc.rest_w24<-cbind(pc.rest_w24,rep(1,nrow(pc.rest_w24)))
colnames(pc.rest_w24)<-c("vhcid","pc_w24","ind_w24")
pc.rest_w24$ind_w24[is.na(pc.rest_w24$pc_w24)==T]<-NA
pc.int_w24$ind_w24[is.na(pc.int_w24$pc_w24)==T]<-NA





data_pc<-merge(pc.int_d1,pc.int_w1,all=T,by="vhcid")
data_pc<-merge(pc.int_w2,data_pc,all=T,by="vhcid")
data_pc<-merge(pc.int_w4,data_pc,all=T,by="vhcid")
data_pc<-merge(pc.int_w6,data_pc,all=T,by="vhcid")
data_pc<-merge(pc.int_w8,data_pc,all=T,by="vhcid")
data_pc<-merge(pc.int_w12,data_pc,all=T,by="vhcid")
data_pc<-merge(pc.int_w16,data_pc,all=T,by="vhcid")
data_pc<-merge(pc.int_w20,data_pc,all=T,by="vhcid")
data_pc<-merge(pc.int_w24,data_pc,all=T,by="vhcid")

data1_pc<-merge(pc.rest_d1,pc.rest_w1,all=T,by="vhcid")
data1_pc<-merge(pc.rest_w2,data1_pc,all=T,by="vhcid")
data1_pc<-merge(pc.rest_w4,data1_pc,all=T,by="vhcid")
data1_pc<-merge(pc.rest_w6,data1_pc,all=T,by="vhcid")
data1_pc<-merge(pc.rest_w8,data1_pc,all=T,by="vhcid")
data1_pc<-merge(pc.rest_w12,data1_pc,all=T,by="vhcid")
data1_pc<-merge(pc.rest_w16,data1_pc,all=T,by="vhcid")
data1_pc<-merge(pc.rest_w20,data1_pc,all=T,by="vhcid")
data1_pc<-merge(pc.rest_w24,data1_pc,all=T,by="vhcid")

pc<-data_pc[,c(2,4,6,8,10,12,14,16,18,20)]
pc<-as.vector(as.matrix(pc))
pc.int<-pc[-which(is.na(pc)==T)]
pc2<-data1_pc[,c(2,4,6,8,10,12,14,16,18,20)]
pc2<-as.vector(as.matrix(pc2))
pc.rest<-pc2[-which(is.na(pc2)==T)]

v1_pc<-pc.rest
v2_pc<-pc.int


rr_d1<-ecdf(pc.int_d1[,2])
rr_w1<-ecdf(pc.int_w1[,2])
rr_w2<-ecdf(pc.int_w2[,2])
rr_w4<-ecdf(pc.int_w4[,2])
rr_w6<-ecdf(pc.int_w6[,2])
rr_w8<-ecdf(pc.int_w8[,2])
rr_w12<-ecdf(pc.int_w12[,2])
rr_w16<-ecdf(pc.int_w16[,2])
rr_w20<-ecdf(pc.int_w20[,2])
rr_w24<-ecdf(pc.int_w24[,2])


kk_d1<-ecdf(pc.rest_d1[,2])
kk_w1<-ecdf(pc.rest_w1[,2])
kk_w2<-ecdf(pc.rest_w2[,2])
kk_w4<-ecdf(pc.rest_w4[,2])
kk_w6<-ecdf(pc.rest_w6[,2])
kk_w8<-ecdf(pc.rest_w8[,2])
kk_w12<-ecdf(pc.rest_w12[,2])
kk_w16<-ecdf(pc.rest_w16[,2])
kk_w20<-ecdf(pc.rest_w20[,2])
kk_w24<-ecdf(pc.rest_w24[,2])



kr.new<-sort(union(knots(kk_d1),union(knots(kk_w1),union(knots(kk_w2),union(knots(kk_w4),union(knots(kk_w6),union(knots(kk_w8),union(knots(kk_w12),union(knots(kk_w16),union(knots(kk_w20),union(knots(kk_w24),union(knots(rr_d1),union(knots(rr_w1),union(knots(rr_w2),union(knots(rr_w4),union(knots(rr_w6),union(knots(rr_w8),union(knots(rr_w12),union(knots(rr_w16),union(knots(rr_w20),knots(rr_w24)))))))))))))))))))))

prob.G_d1<-kk_d1(kr.new)
prob.G_d7<-kk_w1(kr.new)
prob.G_d14<-kk_w2(kr.new)
prob.G_d28<-kk_w4(kr.new)
prob.G_w6<-kk_w6(kr.new)
prob.G_w8<-kk_w8(kr.new)
prob.G_w12<-kk_w12(kr.new)
prob.G_w16<-kk_w16(kr.new)
prob.G_w20<-kk_w20(kr.new)
prob.G_w24<-kk_w24(kr.new)


prob.F_d1<-rr_d1(kr.new)
prob.F_d7<-rr_w1(kr.new)
prob.F_d14<-rr_w2(kr.new)
prob.F_d28<-rr_w4(kr.new)
prob.F_w6<-rr_w6(kr.new)
prob.F_w8<-rr_w8(kr.new)
prob.F_w12<-rr_w12(kr.new)
prob.F_w16<-rr_w16(kr.new)
prob.F_w20<-rr_w20(kr.new)
prob.F_w24<-rr_w24(kr.new)


jumps.F_d1<-c(prob.F_d1[1],prob.F_d1[-1]-prob.F_d1[-length(prob.F_d1)])
jumps.F_d7<-c(prob.F_d7[1],prob.F_d7[-1]-prob.F_d7[-length(prob.F_d7)])
jumps.F_d14<-c(prob.F_d14[1],prob.F_d14[-1]-prob.F_d14[-length(prob.F_d14)])
jumps.F_d28<-c(prob.F_d28[1],prob.F_d28[-1]-prob.F_d28[-length(prob.F_d28)])
jumps.F_w6<-c(prob.F_w6[1],prob.F_w6[-1]-prob.F_w6[-length(prob.F_w6)])
jumps.F_w8<-c(prob.F_w8[1],prob.F_w8[-1]-prob.F_w8[-length(prob.F_w8)])
jumps.F_w12<-c(prob.F_w12[1],prob.F_w12[-1]-prob.F_w12[-length(prob.F_w12)])
jumps.F_w16<-c(prob.F_w16[1],prob.F_w16[-1]-prob.F_w16[-length(prob.F_w16)])
jumps.F_w20<-c(prob.F_w20[1],prob.F_w20[-1]-prob.F_w20[-length(prob.F_w20)])
jumps.F_w24<-c(prob.F_w24[1],prob.F_w24[-1]-prob.F_w24[-length(prob.F_w24)])


jumps.G_d1<-c(prob.G_d1[1],prob.G_d1[-1]-prob.G_d1[-length(prob.G_d1)])
jumps.G_d7<-c(prob.G_d7[1],prob.G_d7[-1]-prob.G_d7[-length(prob.G_d7)])
jumps.G_d14<-c(prob.G_d14[1],prob.G_d14[-1]-prob.G_d14[-length(prob.G_d14)])
jumps.G_d28<-c(prob.G_d28[1],prob.G_d28[-1]-prob.G_d28[-length(prob.G_d28)])
jumps.G_w6<-c(prob.G_w6[1],prob.G_w6[-1]-prob.G_w6[-length(prob.G_w6)])
jumps.G_w8<-c(prob.G_w8[1],prob.G_w8[-1]-prob.G_w8[-length(prob.G_w8)])
jumps.G_w12<-c(prob.G_w12[1],prob.G_w12[-1]-prob.G_w12[-length(prob.G_w12)])
jumps.G_w16<-c(prob.G_w16[1],prob.G_w16[-1]-prob.G_w16[-length(prob.G_w16)])
jumps.G_w20<-c(prob.G_w20[1],prob.G_w20[-1]-prob.G_w20[-length(prob.G_w20)])
jumps.G_w24<-c(prob.G_w24[1],prob.G_w24[-1]-prob.G_w24[-length(prob.G_w24)])

S<-Sys.time()




bs.stat_pc<-rep(0,500)
bs.stat2_pc<-rep(0,500)

for (i in 1:500){

eps1<-rexp(nrow(data_pc),1)
eps2<-rexp(nrow(data1_pc),1)

data_pc_bs<-data_pc
data_pc_bs[,c(3,5,7,9,11,13,15,17,19,21)]<-data_pc_bs[,c(3,5,7,9,11,13,15,17,19,21)]*eps1


data1_pc_bs<-data1_pc
data1_pc_bs[,c(3,5,7,9,11,13,15,17,19,21)]<-data1_pc_bs[,c(3,5,7,9,11,13,15,17,19,21)]*eps2


F_bs<-matrix(0,length(kr.new),10)
for (k in c(2,4,6,8,10,12,14,16,18,20)){ 
F_bs[,k/2]<-sapply(kr.new,function(x){sum((as.vector(na.omit(data_pc_bs[,k]))<=x)*as.vector(na.omit(data_pc_bs[,k+1])))/sum(as.vector(na.omit(data_pc_bs[,k+1])))})
}
Prob_Fh<-rowSums(F_bs)

G_bs<-matrix(0,length(kr.new),10)
for (k in c(2,4,6,8,10,12,14,16,18,20)){ 
G_bs[,k/2]<-sapply(kr.new,function(x){sum((as.vector(na.omit(data1_pc_bs[,k]))<=x)*as.vector(na.omit(data1_pc_bs[,k+1])))/sum(as.vector(na.omit(data1_pc_bs[,k+1])))})
}
Prob_Gh<-rowSums(G_bs)


FG<-cbind(kr.new,prob.F_d1,prob.F_d7,prob.F_d14,prob.F_d28,prob.F_w6,prob.F_w8,prob.F_w12,prob.F_w16,prob.F_w20,prob.F_w24,
	jumps.F_d1,jumps.F_d7,jumps.F_d14,jumps.F_d28,jumps.F_w6,jumps.F_w8,jumps.F_w12,jumps.F_w16,jumps.F_w20,jumps.F_w24,
	prob.G_d1,prob.G_d7,prob.G_d14,prob.G_d28,prob.G_w6,prob.G_w8,prob.G_w12,prob.G_w16,prob.G_w20,prob.G_w24,
	jumps.G_d1,jumps.G_d7,jumps.G_d14,jumps.G_d28,jumps.G_w6,jumps.G_w8,jumps.G_w12,jumps.G_w16,jumps.G_w20,jumps.G_w24)



bs_sum_jump_size<-Prob_Fh-rowSums(FG[,c(2:11)])-Prob_Gh+rowSums(FG[,c(22:31)])
sum_jump_size<-rowSums(FG[,c(2:11)])-rowSums(FG[,c(22:31)])
jump_size<-(length(v2_pc)*rowSums(FG[,c(12:21)])+length(v1_pc)*rowSums(FG[,c(32:41)]))/(length(v2_pc)+length(v1_pc))

bs.stat_pc[i]<-sum(bs_sum_jump_size^2*jump_size)
}

test.stat_pc<-sum(sum_jump_size^2*jump_size)
pval.pc<-mean(bs.stat_pc > test.stat_pc)

E<-Sys.time()-S



###########################################################################################################






########################### NPC ####################################################



npc.int_w1<-le_int[which(le_int$tmpt=='TMT Day 7'),c(2,14)]
npc.int_w1<-cbind(npc.int_w1,rep(1,nrow(npc.int_w1)))
colnames(npc.int_w1)<-c("vhcid","npc_w1","ind_w1")
npc.rest_w1<-le_rest[which(le_rest$tmpt=='TMT Day 7'),c(2,14)]
npc.rest_w1<-cbind(npc.rest_w1,rep(1,nrow(npc.rest_w1)))
colnames(npc.rest_w1)<-c("vhcid","npc_w1","ind_w1")
npc.rest_w1$ind_w1[is.na(npc.rest_w1$npc_w1)==T]<-NA
npc.int_w1$ind_w1[is.na(npc.int_w1$npc_w1)==T]<-NA


npc.int_w2<-le_int[which(le_int$tmpt=='TMT Day 14'),c(2,14)]
npc.int_w2<-cbind(npc.int_w2,rep(1,nrow(npc.int_w2)))
colnames(npc.int_w2)<-c("vhcid","npc_w2","ind_w2")
npc.rest_w2<-le_rest[which(le_rest$tmpt=='TMT Day 14'),c(2,14)]
npc.rest_w2<-cbind(npc.rest_w2,rep(1,nrow(npc.rest_w2)))
colnames(npc.rest_w2)<-c("vhcid","npc_w2","ind_w2")
npc.rest_w2$ind_w2[is.na(npc.rest_w2$npc_w2)==T]<-NA
npc.int_w2$ind_w2[is.na(npc.int_w2$npc_w2)==T]<-NA


npc.int_w4<-le_int[which(le_int$tmpt=='TMT Day 28'),c(2,14)]
npc.int_w4<-cbind(npc.int_w4,rep(1,nrow(npc.int_w4)))
colnames(npc.int_w4)<-c("vhcid","npc_w4","ind_w4")
npc.rest_w4<-le_rest[which(le_rest$tmpt=='TMT Day 28'),c(2,14)]
npc.rest_w4<-cbind(npc.rest_w4,rep(1,nrow(npc.rest_w4)))
colnames(npc.rest_w4)<-c("vhcid","npc_w4","ind_w4")
npc.rest_w4$ind_w4[is.na(npc.rest_w4$npc_w4)==T]<-NA
npc.int_w4$ind_w4[is.na(npc.int_w4$npc_w4)==T]<-NA


npc.int_w6<-le_int[which(le_int$tmpt=='TMT Week 6'),c(2,14)]
npc.int_w6<-cbind(npc.int_w6,rep(1,nrow(npc.int_w6)))
colnames(npc.int_w6)<-c("vhcid","npc_w6","ind_w6")
npc.rest_w6<-le_rest[which(le_rest$tmpt=='TMT Week 6'),c(2,14)]
npc.rest_w6<-cbind(npc.rest_w6,rep(1,nrow(npc.rest_w6)))
colnames(npc.rest_w6)<-c("vhcid","npc_w6","ind_w6")
npc.rest_w6$ind_w6[is.na(npc.rest_w6$npc_w6)==T]<-NA
npc.int_w6$ind_w6[is.na(npc.int_w6$npc_w6)==T]<-NA


npc.int_w8<-le_int[which(le_int$tmpt=='TMT Week 8'),c(2,14)]
npc.int_w8<-cbind(npc.int_w8,rep(1,nrow(npc.int_w8)))
colnames(npc.int_w8)<-c("vhcid","npc_w8","ind_w8")
npc.rest_w8<-le_rest[which(le_rest$tmpt=='TMT Week 8'),c(2,14)]
npc.rest_w8<-cbind(npc.rest_w8,rep(1,nrow(npc.rest_w8)))
colnames(npc.rest_w8)<-c("vhcid","npc_w8","ind_w8")
npc.rest_w8$ind_w8[is.na(npc.rest_w8$npc_w8)==T]<-NA
npc.int_w8$ind_w8[is.na(npc.int_w8$npc_w8)==T]<-NA


npc.int_w12<-le_int[which(le_int$tmpt=='TMT Week 12'),c(2,14)]
npc.int_w12<-cbind(npc.int_w12,rep(1,nrow(npc.int_w12)))
colnames(npc.int_w12)<-c("vhcid","npc_w12","ind_w12")
npc.rest_w12<-le_rest[which(le_rest$tmpt=='TMT Week 12'),c(2,14)]
npc.rest_w12<-cbind(npc.rest_w12,rep(1,nrow(npc.rest_w12)))
colnames(npc.rest_w12)<-c("vhcid","npc_w12","ind_w12")
npc.rest_w12$ind_w12[is.na(npc.rest_w12$npc_w12)==T]<-NA
npc.int_w12$ind_w12[is.na(npc.int_w12$npc_w12)==T]<-NA


npc.int_w16<-le_int[which(le_int$tmpt=='TMT Week 16'),c(2,14)]
npc.int_w16<-cbind(npc.int_w16,rep(1,nrow(npc.int_w16)))
colnames(npc.int_w16)<-c("vhcid","npc_w16","ind_w16")
npc.rest_w16<-le_rest[which(le_rest$tmpt=='TMT Week 16'),c(2,14)]
npc.rest_w16<-cbind(npc.rest_w16,rep(1,nrow(npc.rest_w16)))
colnames(npc.rest_w16)<-c("vhcid","npc_w16","ind_w16")
npc.rest_w16$ind_w16[is.na(npc.rest_w16$npc_w16)==T]<-NA
npc.int_w16$ind_w16[is.na(npc.int_w16$npc_w16)==T]<-NA


npc.int_w20<-le_int[which(le_int$tmpt=='TMT Week 20'),c(2,14)]
npc.int_w20<-cbind(npc.int_w20,rep(1,nrow(npc.int_w20)))
colnames(npc.int_w20)<-c("vhcid","npc_w20","ind_w20")
npc.rest_w20<-le_rest[which(le_rest$tmpt=='TMT Week 20'),c(2,14)]
npc.rest_w20<-cbind(npc.rest_w20,rep(1,nrow(npc.rest_w20)))
colnames(npc.rest_w20)<-c("vhcid","npc_w20","ind_w20")
npc.rest_w20$ind_w20[is.na(npc.rest_w20$npc_w20)==T]<-NA
npc.int_w20$ind_w20[is.na(npc.int_w20$npc_w20)==T]<-NA


npc.int_d1<-le_int[which(le_int$tmpt=='Baseline'),c(2,14)]
npc.int_d1<-cbind(npc.int_d1,rep(1,nrow(npc.int_d1)))
colnames(npc.int_d1)<-c("vhcid","npc_d1","ind_d1")
npc.rest_d1<-le_rest[which(le_rest$tmpt=='Baseline'),c(2,14)]
npc.rest_d1<-cbind(npc.rest_d1,rep(1,nrow(npc.rest_d1)))
colnames(npc.rest_d1)<-c("vhcid","npc_d1","ind_d1")
npc.rest_d1$ind_d1[is.na(npc.rest_d1$npc_d1)==T]<-NA
npc.int_d1$ind_d1[is.na(npc.int_d1$npc_d1)==T]<-NA



npc.int_w24<-le_int[which(le_int$tmpt=='TMT Week 24'),c(2,14)]
npc.int_w24<-cbind(npc.int_w24,rep(1,nrow(npc.int_w24)))
colnames(npc.int_w24)<-c("vhcid","npc_w24","ind_w24")
npc.rest_w24<-le_rest[which(le_rest$tmpt=='TMT Week 24'),c(2,14)]
npc.rest_w24<-cbind(npc.rest_w24,rep(1,nrow(npc.rest_w24)))
colnames(npc.rest_w24)<-c("vhcid","npc_w24","ind_w24")
npc.rest_w24$ind_w24[is.na(npc.rest_w24$npc_w24)==T]<-NA
npc.int_w24$ind_w24[is.na(npc.int_w24$npc_w24)==T]<-NA





data_npc<-merge(npc.int_d1,npc.int_w1,all=T,by="vhcid")
data_npc<-merge(npc.int_w2,data_npc,all=T,by="vhcid")
data_npc<-merge(npc.int_w4,data_npc,all=T,by="vhcid")
data_npc<-merge(npc.int_w6,data_npc,all=T,by="vhcid")
data_npc<-merge(npc.int_w8,data_npc,all=T,by="vhcid")
data_npc<-merge(npc.int_w12,data_npc,all=T,by="vhcid")
data_npc<-merge(npc.int_w16,data_npc,all=T,by="vhcid")
data_npc<-merge(npc.int_w20,data_npc,all=T,by="vhcid")
data_npc<-merge(npc.int_w24,data_npc,all=T,by="vhcid")

data1_npc<-merge(npc.rest_d1,npc.rest_w1,all=T,by="vhcid")
data1_npc<-merge(npc.rest_w2,data1_npc,all=T,by="vhcid")
data1_npc<-merge(npc.rest_w4,data1_npc,all=T,by="vhcid")
data1_npc<-merge(npc.rest_w6,data1_npc,all=T,by="vhcid")
data1_npc<-merge(npc.rest_w8,data1_npc,all=T,by="vhcid")
data1_npc<-merge(npc.rest_w12,data1_npc,all=T,by="vhcid")
data1_npc<-merge(npc.rest_w16,data1_npc,all=T,by="vhcid")
data1_npc<-merge(npc.rest_w20,data1_npc,all=T,by="vhcid")
data1_npc<-merge(npc.rest_w24,data1_npc,all=T,by="vhcid")

npc<-data_npc[,c(2,4,6,8,10,12,14,16,18,20)]
npc<-as.vector(as.matrix(npc))
npc.int<-npc[-which(is.na(npc)==T)]
npc2<-data1_npc[,c(2,4,6,8,10,12,14,16,18,20)]
npc2<-as.vector(as.matrix(npc2))
npc.rest<-npc2[-which(is.na(npc2)==T)]

v1_npc<-npc.rest
v2_npc<-npc.int


rr_d1<-ecdf(npc.int_d1[,2])
rr_w1<-ecdf(npc.int_w1[,2])
rr_w2<-ecdf(npc.int_w2[,2])
rr_w4<-ecdf(npc.int_w4[,2])
rr_w6<-ecdf(npc.int_w6[,2])
rr_w8<-ecdf(npc.int_w8[,2])
rr_w12<-ecdf(npc.int_w12[,2])
rr_w16<-ecdf(npc.int_w16[,2])
rr_w20<-ecdf(npc.int_w20[,2])
rr_w24<-ecdf(npc.int_w24[,2])


kk_d1<-ecdf(npc.rest_d1[,2])
kk_w1<-ecdf(npc.rest_w1[,2])
kk_w2<-ecdf(npc.rest_w2[,2])
kk_w4<-ecdf(npc.rest_w4[,2])
kk_w6<-ecdf(npc.rest_w6[,2])
kk_w8<-ecdf(npc.rest_w8[,2])
kk_w12<-ecdf(npc.rest_w12[,2])
kk_w16<-ecdf(npc.rest_w16[,2])
kk_w20<-ecdf(npc.rest_w20[,2])
kk_w24<-ecdf(npc.rest_w24[,2])



kr.new<-sort(union(knots(kk_d1),union(knots(kk_w1),union(knots(kk_w2),union(knots(kk_w4),union(knots(kk_w6),union(knots(kk_w8),union(knots(kk_w12),union(knots(kk_w16),union(knots(kk_w20),union(knots(kk_w24),union(knots(rr_d1),union(knots(rr_w1),union(knots(rr_w2),union(knots(rr_w4),union(knots(rr_w6),union(knots(rr_w8),union(knots(rr_w12),union(knots(rr_w16),union(knots(rr_w20),knots(rr_w24)))))))))))))))))))))

prob.G_d1<-kk_d1(kr.new)
prob.G_d7<-kk_w1(kr.new)
prob.G_d14<-kk_w2(kr.new)
prob.G_d28<-kk_w4(kr.new)
prob.G_w6<-kk_w6(kr.new)
prob.G_w8<-kk_w8(kr.new)
prob.G_w12<-kk_w12(kr.new)
prob.G_w16<-kk_w16(kr.new)
prob.G_w20<-kk_w20(kr.new)
prob.G_w24<-kk_w24(kr.new)


prob.F_d1<-rr_d1(kr.new)
prob.F_d7<-rr_w1(kr.new)
prob.F_d14<-rr_w2(kr.new)
prob.F_d28<-rr_w4(kr.new)
prob.F_w6<-rr_w6(kr.new)
prob.F_w8<-rr_w8(kr.new)
prob.F_w12<-rr_w12(kr.new)
prob.F_w16<-rr_w16(kr.new)
prob.F_w20<-rr_w20(kr.new)
prob.F_w24<-rr_w24(kr.new)


jumps.F_d1<-c(prob.F_d1[1],prob.F_d1[-1]-prob.F_d1[-length(prob.F_d1)])
jumps.F_d7<-c(prob.F_d7[1],prob.F_d7[-1]-prob.F_d7[-length(prob.F_d7)])
jumps.F_d14<-c(prob.F_d14[1],prob.F_d14[-1]-prob.F_d14[-length(prob.F_d14)])
jumps.F_d28<-c(prob.F_d28[1],prob.F_d28[-1]-prob.F_d28[-length(prob.F_d28)])
jumps.F_w6<-c(prob.F_w6[1],prob.F_w6[-1]-prob.F_w6[-length(prob.F_w6)])
jumps.F_w8<-c(prob.F_w8[1],prob.F_w8[-1]-prob.F_w8[-length(prob.F_w8)])
jumps.F_w12<-c(prob.F_w12[1],prob.F_w12[-1]-prob.F_w12[-length(prob.F_w12)])
jumps.F_w16<-c(prob.F_w16[1],prob.F_w16[-1]-prob.F_w16[-length(prob.F_w16)])
jumps.F_w20<-c(prob.F_w20[1],prob.F_w20[-1]-prob.F_w20[-length(prob.F_w20)])
jumps.F_w24<-c(prob.F_w24[1],prob.F_w24[-1]-prob.F_w24[-length(prob.F_w24)])


jumps.G_d1<-c(prob.G_d1[1],prob.G_d1[-1]-prob.G_d1[-length(prob.G_d1)])
jumps.G_d7<-c(prob.G_d7[1],prob.G_d7[-1]-prob.G_d7[-length(prob.G_d7)])
jumps.G_d14<-c(prob.G_d14[1],prob.G_d14[-1]-prob.G_d14[-length(prob.G_d14)])
jumps.G_d28<-c(prob.G_d28[1],prob.G_d28[-1]-prob.G_d28[-length(prob.G_d28)])
jumps.G_w6<-c(prob.G_w6[1],prob.G_w6[-1]-prob.G_w6[-length(prob.G_w6)])
jumps.G_w8<-c(prob.G_w8[1],prob.G_w8[-1]-prob.G_w8[-length(prob.G_w8)])
jumps.G_w12<-c(prob.G_w12[1],prob.G_w12[-1]-prob.G_w12[-length(prob.G_w12)])
jumps.G_w16<-c(prob.G_w16[1],prob.G_w16[-1]-prob.G_w16[-length(prob.G_w16)])
jumps.G_w20<-c(prob.G_w20[1],prob.G_w20[-1]-prob.G_w20[-length(prob.G_w20)])
jumps.G_w24<-c(prob.G_w24[1],prob.G_w24[-1]-prob.G_w24[-length(prob.G_w24)])

S<-Sys.time()




bs.stat_npc<-rep(0,500)
bs.stat2_npc<-rep(0,500)

for (i in 1:500){

eps1<-rexp(nrow(data_npc),1)
eps2<-rexp(nrow(data1_npc),1)

data_npc_bs<-data_npc
data_npc_bs[,c(3,5,7,9,11,13,15,17,19,21)]<-data_npc_bs[,c(3,5,7,9,11,13,15,17,19,21)]*eps1


data1_npc_bs<-data1_npc
data1_npc_bs[,c(3,5,7,9,11,13,15,17,19,21)]<-data1_npc_bs[,c(3,5,7,9,11,13,15,17,19,21)]*eps2


F_bs<-matrix(0,length(kr.new),10)
for (k in c(2,4,6,8,10,12,14,16,18,20)){ 
F_bs[,k/2]<-sapply(kr.new,function(x){sum((as.vector(na.omit(data_npc_bs[,k]))<=x)*as.vector(na.omit(data_npc_bs[,k+1])))/sum(as.vector(na.omit(data_npc_bs[,k+1])))})
}
Prob_Fh<-rowSums(F_bs)

G_bs<-matrix(0,length(kr.new),10)
for (k in c(2,4,6,8,10,12,14,16,18,20)){ 
G_bs[,k/2]<-sapply(kr.new,function(x){sum((as.vector(na.omit(data1_npc_bs[,k]))<=x)*as.vector(na.omit(data1_npc_bs[,k+1])))/sum(as.vector(na.omit(data1_npc_bs[,k+1])))})
}
Prob_Gh<-rowSums(G_bs)


FG<-cbind(kr.new,prob.F_d1,prob.F_d7,prob.F_d14,prob.F_d28,prob.F_w6,prob.F_w8,prob.F_w12,prob.F_w16,prob.F_w20,prob.F_w24,
	jumps.F_d1,jumps.F_d7,jumps.F_d14,jumps.F_d28,jumps.F_w6,jumps.F_w8,jumps.F_w12,jumps.F_w16,jumps.F_w20,jumps.F_w24,
	prob.G_d1,prob.G_d7,prob.G_d14,prob.G_d28,prob.G_w6,prob.G_w8,prob.G_w12,prob.G_w16,prob.G_w20,prob.G_w24,
	jumps.G_d1,jumps.G_d7,jumps.G_d14,jumps.G_d28,jumps.G_w6,jumps.G_w8,jumps.G_w12,jumps.G_w16,jumps.G_w20,jumps.G_w24)



bs_sum_jump_size<-Prob_Fh-rowSums(FG[,c(2:11)])-Prob_Gh+rowSums(FG[,c(22:31)])
sum_jump_size<-rowSums(FG[,c(2:11)])-rowSums(FG[,c(22:31)])
jump_size<-(length(v2_npc)*rowSums(FG[,c(12:21)])+length(v1_npc)*rowSums(FG[,c(32:41)]))/(length(v2_npc)+length(v1_npc))

bs.stat_npc[i]<-sum(bs_sum_jump_size^2*jump_size)
}

test.stat_npc<-sum(sum_jump_size^2*jump_size)
pval.npc<-mean(bs.stat_npc > test.stat_npc)

E<-Sys.time()-S



###########################################################################################################























bs1<-data[sample(nrow(data),nrow(data),replace=T),]
bs2<-data1[sample(nrow(data1),nrow(data1),replace=T),]


bs1_d1<-as.vector(na.omit(bs1[,12]))
bs1_d7<-as.vector(na.omit(bs1[,14]))
bs1_d14<-as.vector(na.omit(bs1[,10]))
bs1_d28<-as.vector(na.omit(bs1[,8]))
bs1_w8<-as.vector(na.omit(bs1[,6]))
bs1_w12<-as.vector(na.omit(bs1[,4]))
bs1_w24<-as.vector(na.omit(bs1[,2]))


bs2_d1<-as.vector(na.omit(bs2[,12]))
bs2_d7<-as.vector(na.omit(bs2[,14]))
bs2_d14<-as.vector(na.omit(bs2[,10]))
bs2_d28<-as.vector(na.omit(bs2[,8]))
bs2_w8<-as.vector(na.omit(bs2[,6]))
bs2_w12<-as.vector(na.omit(bs2[,4]))
bs2_w24<-as.vector(na.omit(bs2[,2]))


kk.bs_d1<-ecdf(bs2_d1)
kk.bs_d7<-ecdf(bs2_d7)
kk.bs_d14<-ecdf(bs2_d14)
kk.bs_d28<-ecdf(bs2_d28)
kk.bs_w8<-ecdf(bs2_w8)
kk.bs_w12<-ecdf(bs2_w12)
kk.bs_w24<-ecdf(bs2_w24)


rr.bs_d1<-ecdf(bs1_d1)
rr.bs_d7<-ecdf(bs1_d7)
rr.bs_d14<-ecdf(bs1_d14)
rr.bs_d28<-ecdf(bs1_d28)
rr.bs_w8<-ecdf(bs1_w8)
rr.bs_w12<-ecdf(bs1_w12)
rr.bs_w24<-ecdf(bs1_w24)


prob.Gh_d1<-kk.bs_d1(kr.new)
prob.Gh_d7<-kk.bs_d7(kr.new)
prob.Gh_d14<-kk.bs_d14(kr.new)
prob.Gh_d28<-kk.bs_d28(kr.new)
prob.Gh_w8<-kk.bs_w8(kr.new)
prob.Gh_w12<-kk.bs_w12(kr.new)
prob.Gh_w24<-kk.bs_w24(kr.new)


prob.Fh_d1<-rr.bs_d1(kr.new)
prob.Fh_d7<-rr.bs_d7(kr.new)
prob.Fh_d14<-rr.bs_d14(kr.new)
prob.Fh_d28<-rr.bs_d28(kr.new)
prob.Fh_w8<-rr.bs_w8(kr.new)
prob.Fh_w12<-rr.bs_w12(kr.new)
prob.Fh_w24<-rr.bs_w24(kr.new)


jumps.Fh_d1<-c(prob.Fh_d1[1],prob.Fh_d1[-1]-prob.Fh_d1[-length(prob.Fh_d1)])
jumps.Fh_d7<-c(prob.Fh_d7[1],prob.Fh_d7[-1]-prob.Fh_d7[-length(prob.Fh_d7)])
jumps.Fh_d14<-c(prob.Fh_d14[1],prob.Fh_d14[-1]-prob.Fh_d14[-length(prob.Fh_d14)])
jumps.Fh_d28<-c(prob.Fh_d28[1],prob.Fh_d28[-1]-prob.Fh_d28[-length(prob.Fh_d28)])
jumps.Fh_w8<-c(prob.Fh_w8[1],prob.Fh_w8[-1]-prob.Fh_w8[-length(prob.Fh_w8)])
jumps.Fh_w12<-c(prob.Fh_w12[1],prob.Fh_w12[-1]-prob.Fh_w12[-length(prob.Fh_w12)])
jumps.Fh_w24<-c(prob.Fh_w24[1],prob.Fh_w24[-1]-prob.Fh_w24[-length(prob.Fh_w24)])


jumps.Gh_d1<-c(prob.Gh_d1[1],prob.Gh_d1[-1]-prob.Gh_d1[-length(prob.Gh_d1)])
jumps.Gh_d7<-c(prob.Gh_d7[1],prob.Gh_d7[-1]-prob.Gh_d7[-length(prob.Gh_d7)])
jumps.Gh_d14<-c(prob.Gh_d14[1],prob.Gh_d14[-1]-prob.Gh_d14[-length(prob.Gh_d14)])
jumps.Gh_d28<-c(prob.Gh_d28[1],prob.Gh_d28[-1]-prob.Gh_d28[-length(prob.Gh_d28)])
jumps.Gh_w8<-c(prob.Gh_w8[1],prob.Gh_w8[-1]-prob.Gh_w8[-length(prob.Gh_w8)])
jumps.Gh_w12<-c(prob.Gh_w12[1],prob.Gh_w12[-1]-prob.Gh_w12[-length(prob.Gh_w12)])
jumps.Gh_w24<-c(prob.Gh_w24[1],prob.Gh_w24[-1]-prob.Gh_w24[-length(prob.Gh_w24)])


FG<-cbind(kr.new,prob.F_d1,prob.F_d7,prob.F_d14,prob.F_d28,prob.F_w8,prob.F_w12,prob.F_w24,
	jumps.F_d1,jumps.F_d7,jumps.F_d14,jumps.F_d28,jumps.F_w8,jumps.F_w12,jumps.F_w24,
	prob.Fh_d1,prob.Fh_d7,prob.Fh_d14,prob.Fh_d28,prob.Fh_w8,prob.Fh_w12,prob.Fh_w24,
	jumps.Fh_d1,jumps.Fh_d7,jumps.Fh_d14,jumps.Fh_d28,jumps.Fh_w8,jumps.Fh_w12,jumps.Fh_w24,
	prob.G_d1,prob.G_d7,prob.G_d14,prob.G_d28,prob.G_w8,prob.G_w12,prob.G_w24,
	jumps.G_d1,jumps.G_d7,jumps.G_d14,jumps.G_d28,jumps.G_w8,jumps.G_w12,jumps.G_w24,
	prob.Gh_d1,prob.Gh_d7,prob.Gh_d14,prob.Gh_d28,prob.Gh_w8,prob.Gh_w12,prob.Gh_w24,
	jumps.Gh_d1,jumps.Gh_d7,jumps.Gh_d14,jumps.Gh_d28,jumps.Gh_w8,jumps.Gh_w12,jumps.Gh_w24)










