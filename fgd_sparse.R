x=c("RColorBrewer","lattice","plotly","gridExtra","colorspace","Matrix","purrr","SLOPE","sigmoid","Rcpp","grDevices","sfsmisc")
lapply(x, require, character.only = TRUE)
#set seed for reproducability
set.seed(1)
#set dimension, amount of training data, number of iterations and grid for rank(Sigma)
d=100
s=seq(10,100,by=10)
s=round(s)
n=50000
it=10
errsgd=rep(0,length(s)*it)
errsfgd=rep(0,length(s)*it)
errmfgd=rep(0,length(s)*it)
theta_star=runif(d,-10,10)
theta_star=theta_star
for(j in (1:length(s))){
  #generate the U matrix in order to generate the training data
U=runif(d*s[j],-1,1)
U=matrix(U,nrow=d)
#standardize U
U=U/norm(U,"2")
  for(m in (1:it)){
    #generate training data
    x=runif(s[j]*n,-sqrt(3),sqrt(3))
    x=matrix(x,nrow=n)
    x=x%*%t(U)
    ep=rnorm(n)
    y=x%*%theta_star+ep
    theta_sgd=rep(0,d)
    #compute SGD estimator
    for(i in (1:n)){
      theta_sgd=theta_sgd+(1/(s[j]+i))*(y[i]-as.numeric(x[i,]%*%theta_sgd))*x[i,]
    }
    theta_sfgd=rep(0,d)
    #compute FGD(1) estimator
    for(i in (1:n)){
      xi=rnorm(d)
      theta_sfgd=theta_sfgd+(1/(s[j]^2+i))*(y[i]-as.numeric(x[i,]%*%theta_sfgd))*
        as.numeric(x[i,]%*%xi)*xi
    }
    errsgd[(j-1)*it+m]=sum(((theta_sgd-theta_star)%*%U)^2)
    errsfgd[(j-1)*it+m]=sum(((theta_sfgd-theta_star)%*%U)^2)
    theta_mfgd=rep(0,d)
    #compute FGD(d) estimator
    for(i in (1:n)){
      for(k in (0:(s[j]-1))){
        xi=rnorm(d)
        theta_mfgd=theta_mfgd+(1/(s[j]*(1+i)))*(y[i]-as.numeric(x[i,]%*%theta_mfgd))*
          as.numeric(x[i,]%*%xi)*xi
      }
    }
    errmfgd[(j-1)*it+m]=sum(((theta_mfgd-theta_star)%*%U)^2)
  }
}
err=matrix(0,nrow=3,ncol=length(s))
sds=err
#compute average errors and standard deviation over the iterations
for(i in (1:length(s))){
  err[1,i]=mean(errsgd[(1+(i-1)*it):(i*it)])/s[i]
  err[2,i]=mean(errsfgd[(1+(i-1)*it):(i*it)])/s[i]
  err[3,i]=mean(errmfgd[(1+(i-1)*it):(i*it)])/s[i]
  sds[1,i]=sd(errsgd[(1+(i-1)*it):(i*it)])/s[i]
  sds[2,i]=sd(errsfgd[(1+(i-1)*it):(i*it)])/s[i]
  sds[3,i]=sd(errmfgd[(1+(i-1)*it):(i*it)])/s[i]
}
#generate plots
par(mar=c(5,5,3,3))
tred=adjustcolor("red",alpha.f=0.2)
tblue=adjustcolor("blue",alpha.f=0.2)
tblack=adjustcolor("black",alpha.f=0.2)
plot(s,err[1,],col="black",type="l",xaxs="i",ylim=c(min(err),max(err))
     ,lwd=1.5,xlab=expression(Rank(Sigma)),ylab=expression(MSPE/Rank(Sigma)),axes=FALSE)
lines(s,err[2,],col="blue",lwd=1.5)
lines(s,err[3,],col="red",lwd=1.5)
polygon(c(s,rev(s)),c(err[1,],rev((err[1,]+sds[1,]))),col=tblack,border= NA )
polygon(c(s,rev(s)),c(err[1,],rev((err[1,]-sds[1,]))),col=tblack,border= NA )
polygon(c(s,rev(s)),c(err[3,],rev((err[3,]+sds[3,]))),col=tred,border= NA )
polygon(c(s,rev(s)),c(err[3,],rev((err[3,]-sds[3,]))),col=tred,border= NA )
polygon(c(s,rev(s)),c(err[2,],rev((err[2,]+sds[2,]))),col=tblue,border= NA )
polygon(c(s,rev(s)),c(err[2,],rev((err[2,]-sds[2,]))),col=tblue,border= NA )
lines(s,err[1,],col="black",lwd=1.5)
lines(s,err[3,],col="red",lwd=1.5)
lines(s,err[2,],col="blue",lwd=1.5)
eaxis(side=1,at.small=FALSE,at=s,use.expr=FALSE)
eaxis(side=2,at.small=FALSE,at=c(-1,seq(0,2,by=0.5)))

