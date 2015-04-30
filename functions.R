

est2stage<-function(data, beta.init, est.theta,theta.init=NULL, maxit=300,tol=10e-2){
  ##########################################################################
  ## beta= time varying coefficients
  ## theta= time independent coefficients
  ## beta.init= initial values for beta - a px1 vector for each time point
  ## theta.init= initial values for theta, of size qx1, default=NULL
  ## data= a data frame with the first column the id for each person,
  ##       the second column the time points, the third column the response,
  ##       the next p columns for covariates with time varying effect, and 
  ##       the next q columns for covariates with time independent effects.
  ##########################################################################
  
  if (is.null(theta.init)==F){
    iter=0
    converge=0
    theta.t<-theta.init
    beta.t<-beta.init
    while (TRUE) {
      iter=iter+1
      if (iter > maxit) {
        converge = 1
        break
      }
      beta.est<-est.beta(theta.t, data, beta.t)
      theta.est<-est.theta(beta.est,data, theta.t)
      delta1<-mean(abs(beta.est-beta.t))
      delta2<-mean(abs(theta.t-theta.est))
      if (delta1 < tol & delta2 < tol){
        converge = 0
        break
      }
      else {
        beta.t<-beta.est
        theta.t<-theta.est
      }
    }
  }
  else {
    beta.est<-est.beta(theta=NULL,data, beta)
    theta.est<-NULL
  }
  return(list(converge=converge, iteration=iter, beta=beta.est, theta=theta.est, delta1=delta1, delta2=delta2))
}
    
 
  
  
  
  est.beta<-function(theta,data,beta.init){
    
    id<-unique(data[,1])
    time<-unique(data[,2])  
    p<-dim(beta.init)[2]
    q<-dim(data)[2]-p-3
    y<-data[,3]
    X<-data[,4:(p+3),drop=F]
    Z<-data[,(p+4):(p+q+3),drop=F]
    beta<-matrix(0,ncol=ncol(beta.init),nrow=nrow(beta.init))
        
    for (t in 1:length(time)){
      X.t<-X[data[,2]==time[t],,drop=F]
      Z.t<-Z[data[,2]==time[t],,drop=F]
      y.t<-y[data[,2]==time[t]]
      solve.est.eq<-function(b){
        expr<-paste0("X.t[,",1:p,"]*b[",1:p,"]",collapse="+")
        D.t<-as.vector(exp(eval(parse(text=expr))+Z.t%*%theta)/((1+exp(eval(parse(text=expr))+Z.t%*%theta))^2))
        Y.bar.t<-as.vector(y.t-exp(eval(parse(text=expr))+Z.t%*%theta)/(1+exp(eval(parse(text=expr))+Z.t%*%theta)))
        V.inv.t<-as.vector((1+exp(eval(parse(text=expr))+Z.t%*%theta))^2/exp(eval(parse(text=expr))+Z.t%*%theta))
        est.eq1<-X.t*(D.t*V.inv.t*Y.bar.t)
        #est.eq2<-X.t*Y.bar.t
        return(colMeans(est.eq1))
      }
      beta[t,]<-multiroot(solve.est.eq, beta.init[t,])$root
    }
    return(beta)  
  }
  
  
  est.theta<-function(beta,data,theta.init){
    
    id<-unique(data[,1])
    time<-unique(data[,2])  
    q<-length(theta.init)
    p<-dim(data)[2]-q-3
    y<-data[,3]
    X<-data[,4:(p+3),drop=F]
    Z<-data[,(p+4):(p+q+3),drop=F]
    
    solve.est.eq2<-function(theta){
      expr<-paste0("Z.t[,",1:q,"]*theta[",1:q,"]",collapse="+")
      est.eq1<-array(0,c(length(time),length(id),length(theta)))
      est.eq2<-array(0,c(length(time),length(id),length(theta)))
      for (t in 1:length(time)){
        X.t<-X[data[,2]==time[t],,drop=F]
        Z.t<-Z[data[,2]==time[t],,drop=F]
        y.t<-y[data[,2]==time[t]]
        D.t<-as.vector(exp(X.t%*%beta[t,]+eval(parse(text=expr)))/((1+exp(X.t%*%beta[t,]+eval(parse(text=expr))))^2))
        Y.bar.t<-as.vector(y.t-exp(X.t%*%beta[t,]+eval(parse(text=expr)))/(1+exp(X.t%*%beta[t,]+eval(parse(text=expr)))))
        V.inv.t<-as.vector((1+exp(X.t%*%beta[t,]+eval(parse(text=expr))))^2/exp(X.t%*%beta[t,]+eval(parse(text=expr))))
        est.eq1[t,,]<-Z.t*(D.t*V.inv.t*Y.bar.t)
        #est.eq2[t,,]<-Z.t*Y.bar.t
      }
      t_diff<-time[2:length(time)]-time[1:(length(time)-1)]
      est.eq_sum<-(est.eq1[1:(length(time)-1),,]+est.eq1[2:length(time),,])/2
      int.est.eq<-apply(apply(est.eq_sum,c(2,3),function(x) x*t_diff),c(2,3),sum)
      return(colMeans(int.est.eq))
    }
    theta<-multiroot(solve.est.eq2, theta.init)$root
    return(theta)
  }
  
  
  
  

est.theta2<-function(beta,data,theta.init){
  
  id<-unique(data[,1])
  q<-length(theta.init)
  p<-dim(data)[2]-q-3
  y<-data[,3]
  X<-data[,4:(p+3),drop=F]
  Z<-data[,(p+4):(p+q+3),drop=F]
  
  solve.est.eq2<-function(theta){
    int.est.eq<-matrix(0,length(id),length(theta.init))
    for (i in 1:length(id)){
      data.i<-data[data[,1]==id[i],]
      X.i<-X[data[,1]==id[i],,drop=F]
      Z.i<-Z[data[,1]==id[i],,drop=F]
      y.i<-y[data[,1]==id[i]]
      time.i<-unique(data[,2][data[,1]==id[i]])  
      est.eq1<-array(0,c(length(time.i),length(theta.init)))
      est.eq2<-array(0,c(length(time.i),length(theta.init)))
      expr<-paste0("Z.t[,",1:q,"]*theta[",1:q,"]",collapse="+")
      for (t in 1:length(time.i)){
        X.t<-X.i[data.i[,2]==time.i[t],,drop=F]
        Z.t<-Z.i[data.i[,2]==time.i[t],,drop=F]
        y.t<-y.i[data.i[,2]==time.i[t]]
        D.t<-as.vector(exp(X.t%*%beta[t,]+eval(parse(text=expr)))/((1+exp(X.t%*%beta[t,]+eval(parse(text=expr))))^2))
        Y.bar.t<-as.vector(y.t-exp(X.t%*%beta[t,]+eval(parse(text=expr)))/(1+exp(X.t%*%beta[t,]+eval(parse(text=expr)))))
        V.inv.t<-as.vector((1+exp(X.t%*%beta[t,]+eval(parse(text=expr))))^2/exp(X.t%*%beta[t,]+eval(parse(text=expr))))
        est.eq1[t,]<-Z.t*(D.t*V.inv.t*Y.bar.t)
      }
      #if ()
      t_diff.i<-time.i[2:length(time.i)]-time.i[1:(length(time.i)-1)]
      est.eq_sum.i<-(est.eq1[1:(length(time.i)-1),]+est.eq1[2:length(time.i),])/2
      int.est.eq[i,]<-(apply(apply(est.eq_sum.i,2,function(x) x*t_diff.i),2,sum))/(max(time.i)-min(time.i))
    }
  
    return(colMeans(int.est.eq))
  }
  theta<-multiroot(solve.est.eq2, theta.init)$root
  return(theta)
}


cov.theta<-function(data,beta,theta){
  id<-unique(data[,1])
  time<-unique(data[,2])  
  q<-length(theta)
  p<-dim(data)[2]-q-3
  y<-data[,3]
  X<-data[,4:(p+3),drop=F]
  Z<-data[,(p+4):(p+q+3),drop=F]
  
  int.i<-array(0,c(length(id),dim(Z)[2],dim(Z)[2]))
  int.i2<-array(0,c(length(id),dim(Z)[2],dim(Z)[2]))
  for (i in 1:length(id)){
    data.i<-data[data[,1]==id[i],]
    X.i<-X[data[,1]==id[i],,drop=F]
    Z.i<-Z[data[,1]==id[i],,drop=F]
    y.i<-y[data[,1]==id[i]]
    time.i<-unique(data[,2][data[,1]==id[i]])  
    t_diff.i<-time.i[2:length(time.i)]-time.i[1:(length(time.i)-1)]
    integrand<-array(0,c(length(time.i),dim(Z.i.t)[2],dim(Z.i.t)[2]))
    integrand2<-matrix(0,length(time.i),dim(Z.i.t)[2])
    for (t in 1:length(time.i)){
      X.i.t<-X.i[data.i[,2]==time.i[t],,drop=F]
      Z.i.t<-Z.i[data.i[,2]==time.i[t],,drop=F]
      y.i.t<-y.i[data.i[,2]==time.i[t]]
      D.i.t<-as.vector(exp(X.i.t%*%beta[t,]+Z.i.t%*%theta)/((1+exp(X.i.t%*%beta[t,]+Z.i.t%*%theta))^2))
      Y.bar.i.t<-as.vector(y.i.t-exp(X.i.t%*%beta[t,]+Z.i.t%*%theta)/(1+exp(X.i.t%*%beta[t,]+Z.i.t%*%theta)))
      V.inv.i.t<-as.vector((1+exp(X.i.t%*%beta[t,]+Z.i.t%*%theta))^2/exp(X.i.t%*%beta[t,]+Z.i.t%*%theta))
      Z.bar.t<-zbar.t(t,data,beta,theta)
      M.i.t<-D.i.t*V.inv.i.t*Y.bar.i.t
      integrand[t,,]<-(D.i.t^2)*V.inv.i.t*t(Z.i.t)%*%(Z.i.t-t(X.i.t)%*%Z.bar.t)
      integrand2[t,]<-M.i.t*(Z.i.t-t(X.i.t)%*%Z.bar.t)
    }
    int.av<-(integrand[1:(length(time.i)-1),,]+integrand[2:length(time.i),,])/2
    int.av2<-(integrand2[1:(length(time.i)-1),]+integrand2[2:length(time.i),])/2
    int.sum<-array(0,c(length(t_diff.i),dim(Z.i.t)[2],dim(Z.i.t)[2]))
    int.sum2<-matrix(0,length(t_diff.i),dim(Z.i.t)[2])
    for (t in t_diff.i) {
      int.sum[t,,]<-int.av[t,,]*t_diff.i[t]
      int.sum2[t,]<-int.av2[t,]*t_diff.i[t]
    }
    int.i[i,,]<-apply(int.sum,c(2,3),sum)
    k.i<-apply(int.sum2,2,sum)
    int.i2[i,,]<-k.i%*%t(k.i)
  }
  return(list(A=apply(int.i,c(2,3),mean),B=apply(int.i2,c(2,3),mean)))  
  
}



zbar.t<-function(t,data,beta,theta){
  id<-unique(data[,1])
  time<-unique(data[,2])  
  q<-length(theta)
  p<-dim(data)[2]-q-3
  y<-data[,3]
  X<-data[,4:(p+3),drop=F]
  Z<-data[,(p+4):(p+q+3),drop=F]
  X.t<-X[data[,2]==time[t],,drop=F]
  Z.t<-Z[data[,2]==time[t],,drop=F]
  D.t<-as.vector(exp(X.t%*%beta[t,]++Z.t%*%theta)/((1+exp(X.t%*%beta[t,]++Z.t%*%theta))^2))
  V.inv.t<-as.vector((1+exp(X.t%*%beta[t,]++Z.t%*%theta))^2/exp(X.t%*%beta[t,]++Z.t%*%theta))
  XX.prime<-matrix(apply(X.t,1,function(x){x%*%t(x)}),dim(X.t)[1],dim(X.t)[2])
  XZ.den<-matrix(0,dim(X.t)[1],dim(Z.t)[2])
  for (i in 1:dim(X.t)[1]) {
    XZ.den[i,]<-(D.t[i]^2)*V.inv.t[i]*(X.t[i,]%*%t(Z.t[i,]))
    
  }
  zbar.t<-solve(colSums(matrix(mapply(function(x,y,z){(x^2)%*%y%*%z}, D.t,V.inv.t,XX.prime),dim(XX.prime)[1],dim(XX.prime)[2])))%*%colSums(XZ.den)
  return(zbar.t)
}
  

  