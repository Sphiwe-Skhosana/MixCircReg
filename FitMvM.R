##Mixture of von-Mises
##MLE of the mean direction and concentration parameter
dir_con_mle_mixVM=function(theta,weig){
  A=sum(weig*sin(theta))/sum(weig);B=sum(weig*cos(theta))/sum(weig)
  out=NULL
  #print(c(A,B))
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
  R=sum(weig*cos(theta-out))/sum(weig)
  return(c(out,A1inv(R)))
}

##Link function
g=function(x) c(2*atan(x))

mse_bias=function(est,true){
  x=est-matrix(true,nrow(est),ncol(est),byrow=T)
  return(apply(x,2,function(x) c(mean(x),var(x),sd(x))))
}

ColStat=function(x){
  apply(x,2,function(x) c(mean(x),sd(x)))
}

mix.vonMises=function(theta,k=2,initmu=NULL,initlmd=NULL,initprop=NULL){
  n=length(theta)
  ##Initialising
  if(!is.null(initmu)){mu0=initmu}else{mu0=runif(k,-pi,pi)}
  if(!is.null(initlmd)){lmd0=initlmd}else{lmd0=rgamma(k,1,1)}
  if(!is.null(initprop)){prop0=initprop}else{prop0=rep(1/k,k)}
  ####
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) prop0[j]*dvonmises(theta,circular(mu0[j]),lmd0[j])))))
  diff=1e6
  llk=NULL
  while(diff>1e-10){
    ###E-step
    g=sapply(1:k,function(j) prop0[j]*dvonmises(theta,circular(mu0[j]),lmd0[j]))
    gn=g/rowSums(g)
    ##M-step
    prop1=colSums(gn)/n
    mu1=lmd1=NULL
    for(j in 1:k){
      res=dir_con_mle_mixVM(theta,gn[,j])
      mu1=c(mu1,res[1])
      lmd1=c(lmd1,res[2])
    }
    ##Checking for convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) prop1[j]*dvonmises(theta,circular(mu1[j]),lmd1[j])))))
    diff=abs(LogLik1-LogLik0)
    LogLik0=LogLik1
    mu0=mu1
    lmd0=lmd1
    prop0=prop1
    llk=c(llk,LogLik1)
  }
  df=3*k-1
  BIC=-2*LogLik1+log(n)*df
  resp=sapply(1:k,function(j) prop1[j]*dvonmises(theta,circular(mu1[j]),lmd1[j]))
  gn=resp/rowSums(resp)
  z=apply(gn,1,which.max)
  return(list(z=z,mu=mu0,lmd=lmd0,prop=prop0,LL=LogLik1,BIC=BIC))
}

fit.mix_vonMises=function(theta,k=2,initmu=NULL,initlmd=NULL,initprop=NULL){
  count=0
  R=10
  res=list()
  counter=0
  while(count<R){
    fit=NULL
    try({fit=mix.vonMises(theta,k=k,initmu=initmu,initlmd=initlmd,initprop=initprop)},silent=F)
    if(!is.null(fit) & !any(fit$lmd==0)) {count=count+1;res[[count]]=fit;print(count)}
    counter=counter+1
    if(counter==10) count=R
  }
  return(res)
}