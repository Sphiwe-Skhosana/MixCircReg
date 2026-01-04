##Normalizing functions
minmaxscaler=function(x){
  (x-min(x))/(max(x)-min(x))
}

standardizer=function(x){
  z=(x-mean(x))/sd(x)
  return(z)
}

####Mixture of Circular Regressions####
##MLE of the mean direction and concentration parameter
dir_con_mle=function(theta,x,weig,initBeta){
  Beta0=as.matrix(initBeta)
  A=sum(weig*sin(theta-g(x%*%Beta0)))/sum(weig);B=sum(weig*cos(theta-g(x%*%Beta0)))/sum(weig)
  mu=NULL
  if(A>0 & B>0){
    mu=atan(A/B)
  }else{
    if(A<0 & B<0){
      mu=atan(A/B)+pi
    }else{
      if(A<0 & B>0){
        mu=atan(A/B)+2*pi
      }else{
        if(A>0 & B<0){
          mu=atan(A/B)+pi
        }
      }
    }
  }
  R=sum(weig*cos(theta-mu-g(x%*%Beta0)))/sum(weig)
  #R=sqrt(A^2+B^2)
  if(R<=0) R=0
  return(c(mu,A1inv(R)))
}

g=function(x) c(2*atan(x))
g1=function(x) c(2/(1+x^2))
g2=function(x) c(-4*x/(1+x^2)^2)

##mix_circ_regression
mix_f=function(theta,x,prop,Beta,mu,lmd){
  n=length(theta)
  k=length(prop)
  g=NULL
  for(j in 1:k){
    g=cbind(g,sapply(1:n,function(i) prop[j]*dvonmises(theta[i],circular(mu[j]+g(x%*%Beta[,j]))[i],lmd[j])))
  }
  return(g)
}
##MLE of Beta
Beta_IRWLS=function(theta,x,weig,mu,lmd,initBeta){
  n=length(theta)
  w=diag(weig)
  diff=1e6
  Beta_old=as.matrix(initBeta)
  count=0
  sing=0
  while(diff>1e-10){
    G1=diag(g1(x%*%Beta_old))
    G2=diag(g2(x%*%Beta_old))
    u=sin(theta-mu-g(x%*%Beta_old))
    v=cos(theta-mu-g(x%*%Beta_old))
    grad=t(x)%*%G1%*%w%*%u
    H=-t(x)%*%G1^2%*%w%*%diag(v)%*%x+t(x)%*%G2%*%w%*%diag(u)%*%x
    if(abs(det(H))<1e-100){diff=1e-100;sing=1}else{
      Beta_new=Beta_old-solve(H)%*%grad
      diff=max(abs(Beta_new-Beta_old))
      Beta_old=Beta_new
      count=count+1
      if(count==1e2) diff=1e-100 ##Set the max no. of iterations at 100
    }
  }
  return(list(Beta=Beta_old,sing=sing))
}

##
mix_circ_regression=function(theta,x1,x2,k=2,initBeta1=NULL,initBeta2=NULL,initmu=NULL,initlmd=NULL,initprop=NULL){
  ##map the circular response to the interval [0,2pi)
  if(any(theta<0)) theta=theta+pi
  ##
  n=length(theta)
  phi=as.matrix(x1)
  x1=NULL;for(j in 1:ncol(phi)) x1=cbind(x1,cbind(cos(phi[,j]),sin(phi[,j])))
  x2=as.matrix(x2)
  x=cbind(x1,x2)
  q=ncol(x1)
  p=ncol(x2)
  ###Initialize
  if(is.null(initBeta1)){Beta1_0=matrix(runif(q*k,0,1),q,k)}else{Beta1_0=initBeta1}
  if(is.null(initBeta2)){Beta2_0=matrix(runif(p*k,0,1),p,k)}else{Beta2_0=initBeta2}
  if(is.null(initmu)){mu0=runif(k,0.5*pi,1.5*pi)}else{mu0=initmu}
  if(is.null(initlmd)){lmd0=rgamma(k,1,1)}else{lmd0=initlmd}
  if(is.null(initprop)){prop0=rep(1/k,k)}else{prop0=initprop}
  ###
  Beta0=rbind(Beta1_0,Beta2_0)
  LogLik0=sum(log(rowSums(mix_f(theta,x,prop0,Beta0,mu0,lmd0))))
  diff1=1e6
  llk=NULL
  count=0
  while(diff1>1e-10){
    ###E-step
    resp=mix_f(theta,x,prop0,Beta0,mu0,lmd0)
    gn=resp/rowSums(resp)
    ###M-step
    prop1=colSums(gn)/n
    ##Estimating mu and kappa given Beta
    mu1=lmd1=Beta=temp1=temp2=NULL
    for(j in 1:k){
      Beta_old=Beta0[,j]
      diff2=1e10
      countNR=0
      #while(diff2>1e-10){
        res0=dir_con_mle(theta,x,gn[,j],Beta_old)
        mu0=res0[1];lmd0=res0[2];
        res1=Beta_IRWLS(theta,x,gn[,j],mu0,lmd0,Beta_old)
        Beta_new=res1$Beta;
        #if(res2$sing==1) stop("Execution stopped due to singularity...Restarting from new initial values!")
        diff2=max(abs(Beta_new-Beta_old))
        Beta_old=Beta_new
        countNR=countNR+1
        if(countNR==1e2) diff2=1e-100
      #}
      mu1=c(mu1,mu0);lmd1=c(lmd1,lmd0);Beta=cbind(Beta,Beta_old)
    }
    #if(max(temp2)==1) stop("Execution stopped due to bad initial values...Restarting from new initial values!")
    if(any(lmd1<0)) stop("Execution stopped due to negative kappa...Restarting from new initial values!")
    ##Checking for convergence
    LogLik1=sum(log(rowSums(mix_f(theta,x,prop1,Beta,mu1,lmd1))))
    diff1=abs(LogLik1-LogLik0)
    LogLik0=LogLik1
    mu0=mu1
    lmd0=lmd1
    prop0=prop1
    Beta0=Beta
    llk=c(llk,LogLik1)
    count=count+1
    if(count==1e2) diff1=1e-100 ##Set the max no. of iterations at 100
  }
  df=(3+p+q)*k-1
  BIC=-2*LogLik1+log(n)*df
  idx=order(mu1)
  prop1=prop1[idx];lmd1=lmd1[idx];mu1=mu1[idx];Beta1=Beta[1:q,idx];Beta2=matrix(Beta[(q+1):(q+p),idx],p,k)
  resp=mix_f(theta,x,prop1,Beta,mu1,lmd1)
  gn=resp/rowSums(resp)
  z=apply(gn,1,which.max)
  return(list(r=resp,z=z,mu=mu1,prop=prop1,lmd=lmd1,Beta1=Beta1,Beta2=Beta2,LL=LogLik0,BIC=BIC,xcov1=x1,xcov2=x2))
}

fit.mix_circ_regression=function(theta,x1,x2,k=2,initBeta1=NULL,initBeta2=NULL,initmu=NULL,initlmd=NULL,initprop=NULL){
  count=0
  R=10
  res=list()
  counter=0
  while(count<R){
    fit=NULL
    try({fit=mix_circ_regression(theta,x1,x2,k=k,initBeta1=initBeta1,initBeta2 = initBeta2,initmu=initmu,initlmd=initlmd,initprop=initprop)},silent=F)
    if(!is.null(fit) & !any(fit$lmd==0)) {count=count+1;res[[count]]=fit;print(count)}
    counter=counter+1
    if(counter==10) count=R
  }
  return(res)
}