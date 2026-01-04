##MLE of the mean direction and concentration parameter
dir_con_mle1=function(theta,x1,x2,weig,initBeta1,initBeta2){
  Beta1_0=as.matrix(initBeta1)
  Beta2_0=as.matrix(initBeta2)
  A=sum(weig*sin(theta-g(x1%*%Beta1_0)-g(x2%*%Beta2_0)))/sum(weig);B=sum(weig*cos(theta-g(x1%*%Beta1_0)-g(x2%*%Beta2_0)))/sum(weig)
  out=NULL
  if(A>0 & B>0){
    out=atan(A/B)
  }else{
    if(A<0 & B<0){
      out=atan(A/B)+pi
    }else{
      if(A<0 & B>0){
        out=atan(A/B)+2*pi
      }else{
        if(A>0 & B<0){
          out=atan(A/B)+pi
        }
      }
    }
  }
  R=sum(weig*cos(theta-out-g(x1%*%Beta1_0)-g(x2%*%Beta2_0)))/sum(weig)
  #R=sqrt(A^2+B^2)
  return(c(out,A1inv(R)))
}

LogLik=function(theta,x1,x2,Beta1,Beta2,mu,lmd){
  n=length(theta)
  out=sapply(1:n,function(i) dvonmises(theta[i],circular(mu+g(x1%*%Beta1)+g(x2%*%Beta2))[i],lmd))
}

g=function(x) c(2*atan(x))
g1=function(x) c(2/(1+x^2))
g2=function(x) c(-4*x/(1+x^2)^2)
##MLE for Beta
Beta2_IRWLS=function(theta,x1,x2,weig,mu,lmd,initBeta,Beta1){
  n=length(theta)
  w=diag(weig)
  diff=1e6
  Beta_old=as.matrix(initBeta)
  while(diff>1e-10){
    G1=diag(g1(x2%*%Beta_old))
    G2=diag(g2(x2%*%Beta_old))
    u=sin(theta-mu-g(x1%*%Beta1)-g(x2%*%Beta_old))
    v=cos(theta-mu-g(x1%*%Beta1)-g(x2%*%Beta_old))
    grad=t(x2)%*%G1%*%u
    H=-t(x2)%*%G1^2%*%diag(c(v))%*%x2+t(x2)%*%G2%*%diag(c(u))%*%x2
    Beta_new=Beta_old-solve(H)%*%grad
    diff=max(abs(Beta_new-Beta_old))
    Beta_old=Beta_new
  }
  return(Beta_new)
}

Beta1_IRWLS=function(theta,x1,x2,weig,mu,lmd,initBeta,Beta2){
  n=length(theta)
  w=diag(weig)
  diff=1e6
  Beta_old=as.matrix(initBeta)
  while(diff>1e-10){
    G1=diag(g1(x1%*%Beta_old))
    G2=diag(g2(x1%*%Beta_old))
    u=sin(theta-mu-g(x1%*%Beta_old)-g(x2%*%Beta2))
    v=cos(theta-mu-g(x1%*%Beta_old)-g(x2%*%Beta2))
    grad=t(x1)%*%G1%*%u
    H=-t(x1)%*%G1^2%*%diag(c(v))%*%x1+t(x1)%*%G2%*%diag(c(u))%*%x1
    Beta_new=Beta_old-solve(H)%*%grad
    diff=max(abs(Beta_new-Beta_old))
    Beta_old=Beta_new
  }
  return(Beta_new)
}

##
circ_regression=function(theta,x1,x2,initBeta1=NULL,initBeta2=NULL){
  n=length(theta)
  x1=as.matrix(x1)
  phi=NULL;for(j in 1:ncol(x1)) phi=cbind(phi,cbind(cos(x1[,j]),sin(x1[,j])))
  x2=as.matrix(x2)
  q=ncol(phi)
  p=ncol(x2)
  ###Initialize
  if(is.null(initBeta1)){Beta1_0=as.matrix(runif(q,-1,1))}else{Beta1_0=as.matrix(initBeta1)}
  if(is.null(initBeta2)){Beta2_0=as.matrix(runif(p,-1,1))}else{Beta2_0=as.matrix(initBeta2)}
  diff=1e6
  while(diff>1e-10){
  ##Estimating mu and kappa given Beta
  res0=dir_con_mle1(theta,phi,x2,rep(1/n,n),Beta1_0,Beta2_0)
  mu0=res0[1];lmd0=res0[2];
  res1=Beta1_IRWLS(theta,phi,x2,rep(1/n,n),mu0,lmd0,Beta1_0,Beta2_0)
  Beta1=res1;
  res2=Beta2_IRWLS(theta,phi,x2,rep(1/n,n),mu0,lmd0,Beta2_0,Beta1)
  Beta2=res2
  B1=c(Beta1,Beta2);B0=c(Beta1_0,Beta2_0)
  diff=max(abs(B1-B0))
  Beta1_0=Beta1
  Beta2_0=Beta2
  }
  mu=mu0;lmd=lmd0;Beta1=Beta1_0;Beta2=Beta2_0
  LogLik=sum(log(LogLik(theta,phi,x2,Beta1,Beta2,mu,lmd)))
  df=2+p+q
  BIC=-2*LogLik+df*log(n)
  return(list(mu=mu,BIC=BIC,lmd=lmd,Beta1=Beta1,Beta2=Beta2))
}

fit_circ_regression=function(theta,x1,x2,initBeta1=NULL,initBeta2=NULL,initmu=NULL,initlmd=NULL,initprop=NULL){
  count=0
  R=10
  res=list()
  counter=0
  while(count<R){
    fit=NULL
    try({fit=circ_regression(theta,x1,x2)},silent=F)
    if(!is.null(fit) & !any(fit$lmd==0)) {count=count+1;res[[count]]=fit;}
    counter=counter+1
    if(counter==10) count=R
  }
  return(res)
}

###Example: Simulating from a circular regression model using the von-Mises distribution
library(circular)
n=100
mu=pi/2
kp=10
Beta1=as.matrix(c(-1,0.5))
Beta2=as.matrix(c(0.5))
x2=sort(runif(n,-1,1))
x1=runif(n,1,4);phi=cbind(cos(x1),sin(x1))
res=NULL
count=0
while(count<500){
  fit=NULL
  theta=sapply(1:n,function(i) rvonmises(1,(mu+g(phi%*%Beta1)+g(x2%*%Beta2))[i],kp))
  print(any(theta>2*pi))
  plot(x2,theta)
  lines(x2,mu+g(x2%*%Beta2))
  try({fit=circ_regression(theta,x1,x2)},silent = T)
  if(!is.null(fit)){
  lines(x2,fit$mu+g(x2%*%fit$Beta2),col="red")
  res=rbind(res,c(fit$mu,fit$lmd,fit$Beta1,fit$Beta2))
  count=count+1
  }
}

