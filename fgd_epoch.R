x=c("RColorBrewer","lattice","plotly","gridExtra","colorspace","Matrix","purrr","SLOPE",
    "sigmoid","Rcpp","grDevices","sfsmisc")
lapply(x, require, character.only = TRUE)
#set seed for reproducability
set.seed(1)
#set dimension, amount of training data and number of iterations
d=30
n=10000
it=10
#generate logarithmic grid from which the number of epochs/samples is chosen
epoch=10^seq(0,2,by=0.2)
epoch=round(epoch)
errsgd=rep(0,length(epoch)*it)
errsgdep=rep(0,length(epoch)*it)
errsfgd=rep(0,length(epoch)*it)
errmfgd=rep(0,length(epoch)*it)
errmfgdep=rep(0,length(epoch)*it)
errnfgd=rep(0,n*it)
#generate theta_star
theta_star=runif(d,-10,10)
for(q in (1:length(epoch))){
  s=epoch[q]
  for(m in (1:it)){
    x=runif(d*n,-sqrt(3),sqrt(3))
    x=matrix(x,nrow=n)
    ep=rnorm(n)
    y=x%*%theta_star+ep
    theta_sgd=rep(0,d)
    #compute SGD estimator
    for(i in (1:n)){
      theta_sgd=theta_sgd+(1/(d+i))*(y[i]-as.numeric(x[i,]%*%theta_sgd))*x[i,]
      errsgd[(q-1)*it+m]=sum((theta_sgd-theta_star)^2)
    }
    theta_sgdep=rep(0,d)
    #compute SGD estimator with multiple epochs
    for(k in (0:(s-1))){
      for(i in (1:n)){
        theta_sgdep=theta_sgdep+(1/(d+i+k*n))*(y[i]-as.numeric(x[i,]%*%theta_sgdep))*x[i,]
        
      }
    }
    errsgdep[(q-1)*it+m]=sum((theta_sgdep-theta_star)^2)
    #compute FGD(1) estimator with epochs
    theta_sfgd=rep(0,d)
    for(k in (0:(s-1))){
      for(i in (1:n)){
        xi=rnorm(d)
        theta_sfgd=theta_sfgd+(1/(1*(d+i+k*n)))*(y[i]-as.numeric(x[i,]%*%theta_sfgd))*
          as.numeric(x[i,]%*%xi)*xi
      }
    }
    errsfgd[(q-1)*it+m]=sum((theta_mfgd-theta_star)^2)
    #compute FGD(s) estimator
    theta_mfgd=rep(0,d)
    for(i in (1:n)){
      for(k in (0:(s-1))){
        xi=rnorm(d)
        theta_mfgd=theta_mfgd+(1/(s*(d+i)))*(y[i]-as.numeric(x[i,]%*%theta_mfgd))*
          as.numeric(x[i,]%*%xi)*xi
      }
    }
    errmfgd[(q-1)*it+m]=sum((theta_mfgd-theta_star)^2)
    #compute FGD(d) estimator with multiple epochs
    theta_mfgdep=rep(0,d)
    for(u in (0:(s-1))){
    for(i in (1:n)){
      for(k in (0:(d-1))){
        xi=rnorm(d)
        theta_mfgdep=theta_mfgdep+(1/(d*(d+i+u*n)))*(y[i]-as.numeric(x[i,]%*%theta_mfgdep))*
          as.numeric(x[i,]%*%xi)*xi
      }
    } 
    }
    errmfgdep[(q-1)*it+m]=sum((theta_mfgdep-theta_star)^2)
  }
}
#compute errors and standard deviations over the iterations
err=matrix(0,nrow=5,ncol=length(epoch))
sds=matrix(0,nrow=5,ncol=length(epoch))
for(i in (1:length(epoch))){
  err[1,i]=mean(errmfgd[(1+(i-1)*it):(i*it)])
  err[2,i]=mean(errsfgd[(1+(i-1)*it):(i*it)])
  err[3,i]=mean(errsgd[(1+(i-1)*it):(i*it)])
  err[4,i]=mean(errsgdep[(1+(i-1)*it):(i*it)])
  err[5,i]=mean(errmfgdep[(1+(i-1)*it):(i*it)])
}
for(i in (1:length(epoch))){
  sds[1,i]=sd(errmfgd[(1+(i-1)*it):(i*it)])
  sds[2,i]=sd(errsfgd[(1+(i-1)*it):(i*it)])
  sds[3,i]=sd(errsgd[(1+(i-1)*it):(i*it)])
  sds[4,i]=sd(errsgdep[(1+(i-1)*it):(i*it)])
  sds[5,i]=sd(errmfgdep[(1+(i-1)*it):(i*it)])
}
#compute error bars for log log plot
logsds=(sds/err)/log(10)
#generate plots
tred=adjustcolor("red",alpha.f=0.2)
tblue=adjustcolor("blue",alpha.f=0.2)
tblack=adjustcolor("black",alpha.f=0.2)
torange=adjustcolor("orange",alpha.f=0.2)
tgreen=adjustcolor("darkgreen",alpha.f=0.2)
par(mar=c(5,5,3,3))
plot(epoch,err[2,],type="l",lwd=1.5,col="blue",log="xy",ylim=c(min(err),max(err)),
     ,xlab="Number of epochs/samples (log-scale)",xaxs="i",axes="false",ylab="")
lines(epoch,err[1,],lwd=1.5,col="red")
lines(epoch,err[3,],lwd=1.5)
lines(epoch,err[4,],lwd=1.5,col="darkgreen")
lines(epoch,err[5,],lwd=1.5,col="orange")
polygon(c(epoch,rev(epoch)),c(err[1,],rev(exp((log(err[1,])+logsds[1,])))),col=tred,border= NA )
polygon(c(epoch,rev(epoch)),c(err[1,],rev(exp(log(err[1,])-logsds[1,]))),col=tred,border= NA )
polygon(c(epoch,rev(epoch)),c(err[2,],rev(exp((log(err[2,])+logsds[2,])))),col=tblue,border= NA )
polygon(c(epoch,rev(epoch)),c(err[2,],rev(exp(log(err[2,])-logsds[2,]))),col=tblue,border= NA )
polygon(c(epoch,rev(epoch)),c(err[3,],rev(exp((log(err[3,])+logsds[3,])))),col=tblack,border= NA )
polygon(c(epoch,rev(epoch)),c(err[3,],rev(exp(log(err[3,])-logsds[3,]))),col=tblack,border= NA )
polygon(c(epoch,rev(epoch)),c(err[4,],rev(exp((log(err[4,])+logsds[4,])))),col=tgreen,border= NA )
polygon(c(epoch,rev(epoch)),c(err[4,],rev(exp(log(err[4,])-logsds[4,]))),col=tgreen,border= NA )
polygon(c(epoch,rev(epoch)),c(err[5,],rev(exp((log(err[5,])+logsds[5,])))),col=torange,border= NA )
polygon(c(epoch,rev(epoch)),c(err[5,],rev(exp(log(err[5,])-logsds[5,]))),col=torange,border= NA )
eaxis(side=1,at.small=FALSE,use.expr=FALSE)
eaxis(side=2,at.small=FALSE,use.expr=TRUE)
mtext(side=2,text="MSE (log-scale)",line=4)
legend(17,10^9,cex=0.8,legend=c("SGD","FGD(\u2113)","FGD(1) epochs","FGD(d) epochs","SGD epochs"),col=c("black","blue","red","orange","darkgreen"),lty=c(1,1,1,1,1), lwd=1.5)



